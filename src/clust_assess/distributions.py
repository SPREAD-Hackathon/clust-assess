"""Within-group distance distributions of a Pling run, for subcommunities and containment communities.

For every group (a Pling subcommunity, or a containment community) this collects the DCJ-Indel and
containment distances between all pairs of its members.

* A pair of members that is missing from a Pling distance file was filtered out before that distance
  was computed: for DCJ it failed the containment filter; for containment it was removed by the
  sourmash prefilter (only when Pling ran with --sourmash). Such pairs are counted as censored
  (+inf, "too different"), as in clust-assess, not dropped. Dropping them would make groups look tighter.
* Hub plasmids (hub_plasmids.csv) are removed by Pling before subcommunity detection, so they belong to
  no subcommunity. They are also left out of community distributions unless include_hubs is set.
* Quantiles use the "higher" rule, so the median equals statistics.median_high, as in clust-assess.
"""
from __future__ import annotations

import json
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd

from .accio_io import community_of

log = logging.getLogger("clust_assess.distributions")

LEVELS = ("subcommunity", "community")
METRICS = ("dcj", "containment")
QUANTILES = {"min": 0.0, "q25": 0.25, "median": 0.5, "q75": 0.75, "max": 1.0}

# Pling's own layout first, then a flat copy (e.g. /shared/team/pling_plsdb).
FILES = {
    "typing": ["dcj_thresh_*_graph/objects/typing.tsv", "typing.tsv"],
    "hubs": ["dcj_thresh_*_graph/objects/hub_plasmids.csv", "hub_plasmids.csv"],
    "communities": ["containment/containment_communities/objects/communities.tsv", "communities.tsv"],
    "dcj": ["all_plasmids_distances.tsv", "all_plasmid_distances.tsv"],
    "containment": ["containment/all_pairs_containment_distance.tsv", "all_pairs_containment_distance.tsv"],
}


def locate(pling_dir=None, **override) -> dict[str, Path | None]:
    """Paths of the Pling files; an explicit path in `override` wins over the search in `pling_dir`."""
    d = Path(pling_dir) if pling_dir else None
    out: dict[str, Path | None] = {}
    for key, patterns in FILES.items():
        if override.get(key):
            out[key] = Path(override[key])
            continue
        out[key] = None
        for pat in patterns if d else []:
            hits = sorted(d.glob(pat))
            if len(hits) > 1:
                raise ValueError(f"several {key} files in {d}: {[str(h) for h in hits]}; pass one with --{key}")
            if hits:
                out[key] = hits[0]
                break
    for key in ("typing", "dcj"):
        if out[key] is None or not out[key].exists():
            raise FileNotFoundError(f"Pling {key} file not found (looked in {d} for {FILES[key]})")
    return out


def membership(files: dict, include_hubs: bool = False) -> pd.DataFrame:
    """One row per plasmid: plasmid, subcommunity, community, hub (None where a plasmid has no group)."""
    t = pd.read_csv(files["typing"], sep="\t", dtype=str).iloc[:, :2]
    t.columns = ["plasmid", "subcommunity"]
    hubs: set[str] = set()
    if files.get("hubs"):
        hubs = set(pd.read_csv(files["hubs"], dtype=str).iloc[:, 0].dropna())
    if files.get("communities"):
        c = pd.read_csv(files["communities"], sep="\t", dtype=str).iloc[:, :2]
        c.columns = ["plasmid", "community"]
        m = c.merge(t, on="plasmid", how="outer")
    else:
        m = t.copy()
        m["community"] = None
    from_type = m["subcommunity"].map(community_of, na_action="ignore")
    clash = m["community"].notna() & from_type.notna() & (m["community"] != from_type)
    if clash.any():
        log.warning("%d plasmids have a community that differs from their subcommunity prefix", int(clash.sum()))
    m["community"] = m["community"].fillna(from_type)
    missing_hubs = sorted(hubs - set(m["plasmid"]))
    if missing_hubs:
        m = pd.concat([m, pd.DataFrame({"plasmid": missing_hubs})], ignore_index=True)
    m["hub"] = m["plasmid"].isin(hubs)
    m.loc[m["hub"], "subcommunity"] = None
    if not include_hubs:
        m.loc[m["hub"], "community"] = None
    m = m.astype(object).where(m.notna(), None)
    return m.drop_duplicates("plasmid").reset_index(drop=True)


def _read_within(path: Path, comm: dict, sub: dict, chunksize: int) -> pd.DataFrame:
    """Rows of a Pling distance file whose two plasmids share a community or a subcommunity."""
    keep = []
    for ch in pd.read_csv(path, sep="\t", usecols=[0, 1, 2], dtype={0: str, 1: str}, chunksize=chunksize):
        ch.columns = ["a", "b", "distance"]
        ch = ch[ch["a"] != ch["b"]]
        c1, c2 = ch["a"].map(comm), ch["b"].map(comm)
        s1, s2 = ch["a"].map(sub), ch["b"].map(sub)
        keep.append(ch[(c1.notna() & (c1 == c2)) | (s1.notna() & (s1 == s2))])
    df = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame(columns=["a", "b", "distance"])
    swap = df["a"] > df["b"]
    df.loc[swap, ["a", "b"]] = df.loc[swap, ["b", "a"]].to_numpy()
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce").replace(np.inf, np.nan)
    dup = df.duplicated(["a", "b"], keep=False)
    if dup.any():
        n_conflict = int((df[dup].groupby(["a", "b"])["distance"].nunique() > 1).sum())
        log.warning("%s: %d duplicated pairs (%d with different values); keeping the first",
                    path.name, int(df.duplicated(["a", "b"]).sum()), n_conflict)
    return df.drop_duplicates(["a", "b"]).reset_index(drop=True)


def within_pairs(files: dict, members: pd.DataFrame, chunksize: int = 2_000_000) -> pd.DataFrame:
    """All observed within-group pairs, one row each: plasmid_1, plasmid_2, community, subcommunity, dcj,
    containment. community/subcommunity is None when the two plasmids are not in the same one; a distance is
    NaN when Pling did not compute it for that pair."""
    comm = {p: c for p, c in zip(members["plasmid"], members["community"]) if c is not None}
    sub = {p: s for p, s in zip(members["plasmid"], members["subcommunity"]) if s is not None}
    out = None
    for metric in METRICS:
        if not files.get(metric):
            log.warning("no %s distance file: %s distributions are skipped", metric, metric)
            continue
        df = _read_within(Path(files[metric]), comm, sub, chunksize).rename(columns={"distance": metric})
        out = df if out is None else out.merge(df, on=["a", "b"], how="outer")
    for metric in METRICS:
        if metric not in out:
            out[metric] = np.nan
    c1, c2, s1, s2 = (out[x].map(d) for x, d in (("a", comm), ("b", comm), ("a", sub), ("b", sub)))
    out["community"] = c1.where(c1 == c2)
    out["subcommunity"] = s1.where(s1 == s2)
    out = out.rename(columns={"a": "plasmid_1", "b": "plasmid_2"})
    out = out[["plasmid_1", "plasmid_2", "community", "subcommunity", "dcj", "containment"]]
    return out.sort_values(["community", "subcommunity", "plasmid_1", "plasmid_2"], na_position="last") \
              .reset_index(drop=True)


def quantile_higher(observed: np.ndarray, n_total: int, p: float) -> float:
    """p-quantile of the sorted `observed` values padded with (n_total - len(observed)) +inf values, using the
    'higher' rule: element ceil(p * (n_total - 1)). p = 0.5 gives statistics.median_high."""
    if n_total <= 0:
        return math.nan
    i = math.ceil(p * (n_total - 1))
    return float(observed[i]) if i < len(observed) else math.inf


def summarise(pairs: pd.DataFrame, members: pd.DataFrame, level: str, thresholds: dict,
              available: tuple = METRICS, min_members: int = 2) -> pd.DataFrame:
    """One row per group of `level` with at least `min_members` members."""
    sizes = members.dropna(subset=[level]).groupby(level).size()
    sizes = sizes[sizes >= min_members]
    by_group = {g: d for g, d in pairs.dropna(subset=[level]).groupby(level)}
    rows = []
    for g, n in sizes.items():
        n_pairs = n * (n - 1) // 2
        row = {level: g, "n_members": int(n), "n_pairs": n_pairs}
        if level == "subcommunity":
            row["community"] = community_of(g)
        d = by_group.get(g)
        for metric in METRICS:
            obs = np.sort(d[metric].dropna().to_numpy(float)) if d is not None else np.array([])
            if len(obs) > n_pairs:
                raise ValueError(f"{g}: {len(obs)} {metric} pairs observed but only {n_pairs} possible")
            if metric not in available:
                row.update({f"{metric}_n": np.nan, f"{metric}_censored": np.nan})
                row.update({f"{metric}_{q}": np.nan for q in QUANTILES})
                row[f"{metric}_frac_le_threshold"] = np.nan
                continue
            row[f"{metric}_n"] = len(obs)
            row[f"{metric}_censored"] = n_pairs - len(obs)
            for q, p in QUANTILES.items():
                row[f"{metric}_{q}"] = quantile_higher(obs, n_pairs, p)
            row[f"{metric}_frac_le_threshold"] = float((obs <= thresholds[metric]).sum()) / n_pairs
        rows.append(row)
    cols = [level, "n_members", "n_pairs"] + (["community"] if level == "subcommunity" else [])
    for metric in METRICS:
        cols += [f"{metric}_n", f"{metric}_censored"] + [f"{metric}_{q}" for q in QUANTILES] \
                + [f"{metric}_frac_le_threshold"]
    out = pd.DataFrame(rows, columns=cols)
    return out.sort_values(["n_members", level], ascending=[False, True]).reset_index(drop=True)


# ----------------------------------------------------------------------------------------------- plots
INK, INK_2, MUTED, SURFACE, SERIES, CENSORED = "#0b0b0b", "#52514e", "#8a8984", "#fcfcfb", "#2a78d6", "#b9b8b2"
DCJ_CAP = 50  # DCJ values above this share the last histogram bar


def plot(pairs: pd.DataFrame, summary: pd.DataFrame, level: str, thresholds: dict, out_png, top: int = 30,
         title: str = "") -> bool:
    """One column per available metric (DCJ, containment): a within-group histogram (censored pairs as a separate
    grey bar; log scale when they dwarf the computed pairs) above box plots of the `top` largest groups
    (computed pairs only). Returns False, writing nothing, when there are no distances to plot."""
    metrics = [m for m in METRICS if len(summary) and summary[f"{m}_n"].notna().any()]
    if not metrics:
        log.warning("no within-%s distances: %s not written", level, out_png)
        return False
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "axes.titlecolor": INK,
                         "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})
    big = summary.head(top)
    box_h = 0.6 + 0.24 * max(len(big), 3)
    fig, axes = plt.subplots(2, len(metrics), figsize=(15 if len(metrics) > 1 else 10, 3.6 + box_h),
                             layout="constrained", squeeze=False, gridspec_kw={"height_ratios": [3.2, box_h]})
    within = pairs.dropna(subset=[level])
    for col, metric in enumerate(metrics):
        _histogram(axes[0, col], within[metric].dropna().to_numpy(float), int(summary[f"{metric}_censored"].sum()),
                   metric, thresholds[metric], level)
        _boxes(axes[1, col], within, big, level, metric, thresholds[metric])
    missing = [m for m in METRICS if m not in metrics]
    note = f"  (no {' or '.join(missing)} distances in the input)" if missing else ""
    fig.suptitle((title or f"Within-{level} distances") + note, color=INK, fontsize=11, x=0.01, ha="left")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return True


def _style(ax, grid_axis):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid_axis, color="#e4e3df", linewidth=0.6)
    ax.set_axisbelow(True)


def _histogram(ax, obs, n_cens, metric, threshold, level):
    if metric == "dcj":
        vals = np.minimum(obs, DCJ_CAP)
        edges = np.arange(0, DCJ_CAP + 2) - 0.5
        cens_x, width = DCJ_CAP + 8, 1.0
        ticks = list(range(0, DCJ_CAP + 1, 10))
        labels = [str(t) if t < DCJ_CAP else f"{DCJ_CAP}+" for t in ticks]
        xlabel = "DCJ-Indel distance"
    else:
        vals = obs
        edges = np.linspace(0, 1, 51)
        cens_x, width = 1.15, 0.02
        ticks = [0, 0.25, 0.5, 0.75, 1.0]
        labels = [f"{t:g}" for t in ticks]
        xlabel = "containment distance"
    from matplotlib.ticker import FuncFormatter

    counts, _, _ = ax.hist(vals, bins=edges, color=SERIES, edgecolor=SURFACE, linewidth=0.8)
    if n_cens:
        ax.bar([cens_x], [n_cens], width=width * 1.6, color=CENSORED, edgecolor=SURFACE, linewidth=0.8)
        ticks, labels = ticks + [cens_x], labels + ["not\ncomputed"]
    # A linear axis would flatten the computed pairs to nothing when the not-computed bar is far taller
    log_y = n_cens > 5 * max(counts.max() if len(counts) else 0, 1)
    if log_y:
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.axvline(threshold, color=INK_2, linestyle=(0, (4, 3)), linewidth=1)
    ax.text(threshold, 1.0, f" Pling threshold {threshold:g}", transform=ax.get_xaxis_transform(),
            color=INK_2, va="top", ha="left", fontsize=8)
    ax.set_xticks(ticks, labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("pairs (log scale)" if log_y else "pairs")
    n = len(obs) + n_cens
    pct = 100 * n_cens / n if n else 0
    ax.set_title(f"{metric.upper() if metric == 'dcj' else 'Containment'}: {n:,} within-{level} pairs, "
                 f"{pct:.0f}% not computed", loc="left", fontsize=9.5)
    _style(ax, "y")


def _boxes(ax, within, big, level, metric, threshold):
    data, labels = [], []
    for r in big.itertuples(index=False):
        g = getattr(r, level)
        v = within.loc[within[level] == g, metric].dropna().to_numpy(float)
        data.append(v if len(v) else np.array([np.nan]))
        cens = getattr(r, f"{metric}_censored")
        labels.append(f"{g}  (n={r.n_members}, {100 * cens / r.n_pairs:.0f}% n.c.)")
    if not data:
        ax.set_axis_off()
        return
    pos = np.arange(len(data))[::-1]
    style = dict(positions=pos, widths=0.6, patch_artist=True, showfliers=True,
                 medianprops={"color": INK, "linewidth": 1.2},
                 boxprops={"facecolor": "#b7d3f6", "edgecolor": SERIES, "linewidth": 1},
                 whiskerprops={"color": SERIES}, capprops={"color": SERIES},
                 flierprops={"marker": "o", "markersize": 3, "markerfacecolor": SERIES,
                             "markeredgecolor": SURFACE, "alpha": 0.6})
    try:
        ax.boxplot(data, orientation="horizontal", **style)   # matplotlib >= 3.10
    except TypeError:
        ax.boxplot(data, vert=False, **style)
    ax.set_yticks(pos, labels, fontsize=7.5)
    ax.axvline(threshold, color=INK_2, linestyle=(0, (4, 3)), linewidth=1)
    ax.set_xlabel(("DCJ-Indel" if metric == "dcj" else "containment") + " distance (computed pairs only)")
    if metric == "containment":
        ax.set_xlim(-0.02, 1.02)
    plural = "subcommunities" if level == "subcommunity" else "communities"
    ax.set_title(f"Largest {len(data)} {plural} (n.c. = pairs not computed)", loc="left", fontsize=9.5)
    _style(ax, "x")


# ----------------------------------------------------------------------------------------------- driver
def run(outdir, pling_dir=None, include_hubs=False, dcj_threshold=4, containment_threshold=0.5,
        min_members=2, plots=False, top=30, chunksize=2_000_000, **override) -> dict:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = locate(pling_dir, **override)
    members = membership(files, include_hubs=include_hubs)
    pairs = within_pairs(files, members, chunksize=chunksize)
    thresholds = {"dcj": dcj_threshold, "containment": containment_threshold}
    available = tuple(m for m in METRICS if files.get(m))
    pairs.to_csv(out / "within_pairs.tsv.gz", sep="\t", index=False, na_rep="NA")
    info = {"files": {k: str(Path(v).resolve()) if v else None for k, v in files.items()},
            "include_hubs": include_hubs, "thresholds": thresholds, "min_members": min_members,
            "n_plasmids": int(len(members)), "n_hubs": int(members["hub"].sum()),
            "n_within_pairs_observed": int(len(pairs)), "levels": {}}
    for level in LEVELS:
        s = summarise(pairs, members, level, thresholds, available, min_members)
        s.to_csv(out / f"{level}_summary.tsv", sep="\t", index=False, na_rep="NA")
        info["levels"][level] = {
            "n_groups": int(members[level].dropna().nunique()), f"n_groups_ge_{min_members}_members": int(len(s)),
            **{f"{m}_pairs_censored": int(s[f"{m}_censored"].sum()) for m in available}}
        if plots and plot(pairs, s, level, thresholds, out / f"{level}_distribution.png", top=top):
            info["levels"][level]["figure"] = f"{level}_distribution.png"
    (out / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    return info
