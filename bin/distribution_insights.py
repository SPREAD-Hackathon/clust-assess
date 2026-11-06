#!/usr/bin/env python3
"""More views of a distance_distribution.py output folder. Standalone: needs only numpy, pandas and matplotlib.

Reads run_info.json, subcommunity_summary.tsv, community_summary.tsv and within_pairs.tsv.gz from --dist_dir and
writes four figures (default: into the same folder):

  insight_group_sizes.png         How many plasmids sit in groups of each size. A plasmid alone in its subcommunity
                                  has no within-group distances, so there is no baseline to compare an input with.
  insight_same_vs_other.png       Cumulative share of pairs within distance x: pairs in the same subcommunity vs
                                  pairs that only share a community. A curve ends below 100% by the share of pairs
                                  that were never computed.
  insight_dcj_vs_containment.png  Every within-community pair, DCJ-Indel against containment (only when the Pling run
                                  had containment distances).
  insight_heatmaps.png            Distance matrices of the largest communities: DCJ-Indel above the diagonal,
                                  containment below, members grouped by subcommunity (boxed).

Pairs missing from the Pling output were never computed (for DCJ: the pair failed the containment filter). They are
counted from the group sizes in the summaries, so this also works when within_pairs.tsv.gz holds only computed pairs
(Pling runs without a containment file, such as the published PLSDB run).
Each plasmid's subcommunity comes from the typing file recorded in run_info.json (or --typing). Without it, it is
taken from the pairs, and a plasmid alone in its subcommunity is drawn without a subcommunity box.

Example:
  python bin/distribution_insights.py --dist_dir database/danish_distribution
"""
import argparse
import json
import logging
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("distribution_insights")

LEVELS = ("subcommunity", "community")
METRICS = ("dcj", "containment")
NAME = {"dcj": "DCJ-Indel", "containment": "containment"}
SAME, OTHER = "same subcommunity", "same community, other subcommunity"
SERIES = {SAME: "#2a78d6", OTHER: "#eb6834"}          # categorical slots 1 and 2, in that order
NOT_COMPUTED = "#b9b8b2"
INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]     # sequential, light = 0
ORANGES = ["#fde6da", "#f8c3a8", "#f29a74", "#eb6834", "#c4501f", "#963c17", "#682a10"]
DCJ_CAP = 50
SIZE_BINS = [(1, 1), (2, 2), (3, 4), (5, 9), (10, 19), (20, 49), (50, 99), (100, 499), (500, 999), (1000, 4999),
             (5000, 9999), (10000, math.inf)]
FIGURES = {"sizes": "insight_group_sizes.png", "same_vs_other": "insight_same_vs_other.png",
           "scatter": "insight_dcj_vs_containment.png", "heatmaps": "insight_heatmaps.png"}


# ----------------------------------------------------------------------------------------------- data
def load(dist_dir: Path):
    info = json.loads((dist_dir / "run_info.json").read_text())
    summaries = {lvl: pd.read_csv(dist_dir / f"{lvl}_summary.tsv", sep="\t") for lvl in LEVELS}
    pairs = pd.read_csv(dist_dir / "within_pairs.tsv.gz", sep="\t",
                        dtype={"plasmid_1": str, "plasmid_2": str, "community": "category",
                               "subcommunity": "category", "dcj": float, "containment": float})
    pairs["relation"] = pd.Categorical(np.where(pairs["subcommunity"].notna(), SAME, OTHER), categories=[SAME, OTHER])
    return info, summaries, pairs


def pair_totals(summaries: dict, pairs: pd.DataFrame) -> dict:
    """All within-group pairs per relation, computed or not (from the group sizes in the summaries)."""
    same = int(summaries["subcommunity"]["n_pairs"].sum())
    other = int(summaries["community"]["n_pairs"].sum()) - same
    seen = pairs["relation"].value_counts()
    return {SAME: max(same, int(seen.get(SAME, 0))), OTHER: max(other, int(seen.get(OTHER, 0)))}


def membership(pairs: pd.DataFrame, typing) -> pd.DataFrame:
    """plasmid, community, subcommunity (NaN when unknown) for every plasmid in a community."""
    long = pd.concat([pairs[[c, "community", "subcommunity"]].set_axis(["plasmid", "community", "subcommunity"], axis=1)
                      for c in ("plasmid_1", "plasmid_2")], ignore_index=True)
    comm = long.dropna(subset=["community"]).drop_duplicates("plasmid").set_index("plasmid")["community"].astype(str)
    if typing is not None:
        t = pd.read_csv(typing, sep="\t", dtype=str).iloc[:, :2]
        t.columns = ["plasmid", "subcommunity"]
        named = t["subcommunity"].str.extract(r"^(community_\d+)_subcommunity_\d+$", expand=False)
        t["community"] = t["plasmid"].map(comm).fillna(named)   # same rule as within_pairs.tsv.gz
        extra = comm[~comm.index.isin(t["plasmid"])]           # hubs placed in a community (--include_hubs)
        m = pd.concat([t, pd.DataFrame({"plasmid": extra.index, "community": extra.to_numpy(),
                                        "subcommunity": np.nan})], ignore_index=True)
    else:
        sub = long.dropna(subset=["subcommunity"]).drop_duplicates("plasmid").set_index("plasmid")["subcommunity"]
        m = pd.DataFrame({"plasmid": comm.index, "community": comm.to_numpy()})
        m["subcommunity"] = m["plasmid"].map(sub.astype(str))
    return m[["plasmid", "community", "subcommunity"]]


def group_sizes(level, summaries, info, members, exact):
    """Member count of every group at `level`, singletons included; plus a note when that is not possible."""
    if exact:
        return members[level].dropna().value_counts().to_numpy(), ""
    sizes = summaries[level]["n_members"].to_numpy(int)
    n_small = int(info["levels"][level]["n_groups"]) - len(sizes)
    min_members = int(info.get("min_members", 2))
    if min_members <= 2:
        return np.concatenate([sizes, np.ones(n_small, int)]), ""
    return sizes, f"{n_small:,} groups with fewer than {min_members} members are not shown"


# ----------------------------------------------------------------------------------------------- figures
def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2, "axes.titlecolor": INK,
                         "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "hatch.linewidth": 0.6,
                         "legend.frameon": False})
    return plt


def _style(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def _threshold(ax, x, label, vertical=True):
    if vertical:
        ax.axvline(x, color=INK_2, linestyle=(0, (4, 3)), linewidth=1)
        ax.text(x, 1.0, f" {label}", transform=ax.get_xaxis_transform(), color=INK_2, va="top", ha="left",
                fontsize=8)
    else:
        ax.axhline(x, color=INK_2, linestyle=(0, (4, 3)), linewidth=1)
        ax.text(1.0, x, f"{label} ", transform=ax.get_yaxis_transform(), color=INK_2, va="bottom", ha="right",
                fontsize=8)


def plot_sizes(sizes: dict, info: dict, out) -> str:
    plt = _setup()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), layout="constrained")
    lines = []
    for ax, level in zip(axes, LEVELS):
        s, note = sizes[level]
        rows = [(lo, hi, int(s[(s >= lo) & (s <= hi)].sum()), int(((s >= lo) & (s <= hi)).sum())) for lo, hi in SIZE_BINS]
        while len(rows) > 1 and rows[-1][3] == 0:
            rows.pop()
        x = np.arange(len(rows))
        plasmids = [r[2] for r in rows]
        ax.bar(x, plasmids, width=0.72, color=[NOT_COMPUTED if r[0] == 1 else SERIES[SAME] for r in rows],
               edgecolor=SURFACE, linewidth=0.8)
        for xi, v in zip(x, plasmids):
            if v:
                ax.annotate(f"{v:,}", (xi, v), xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
                            fontsize=8, color=INK_2)
        labels = []
        for lo, hi, _, n_groups in rows:
            span = f"{lo:,}" if lo == hi else (f"{lo:,}+" if math.isinf(hi) else f"{lo:,}–{hi:,}")
            labels.append(f"{span}\n{n_groups:,} group{'s' if n_groups != 1 else ''}")
        ax.set_xticks(x, labels, fontsize=8)
        ax.set_xlabel(f"members in the plasmid's {level}")
        ax.set_ylabel("plasmids")
        total, alone = int(s.sum()), int(s[s == 1].sum())
        pct = alone / total if total else 0
        ax.set_title(f"{level.capitalize()}: {alone:,} of {total:,} plasmids ({pct:.0%}) are alone"
                     + (f"\n{note}" if note else ""), loc="left", fontsize=10)
        ax.margins(y=0.12)
        _style(ax)
        lines.append(f"{level}: {alone:,} of {total:,} plasmids ({pct:.0%}) alone; "
                     f"{int((s >= 2).sum()):,} groups with >= 2 members, largest {int(s.max()) if len(s) else 0:,}")
    hubs = int(info.get("n_hubs", 0))
    hub_note = f"; {hubs:,} hub plasmids belong to no subcommunity" if hubs and not info.get("include_hubs") else ""
    fig.suptitle("Group sizes. Grey: plasmids alone in their group, so there are no within-group distances to "
                 f"compare an input with{hub_note}", color=INK, fontsize=11, x=0.01, ha="left")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return "\n".join(lines)


def plot_same_vs_other(pairs, totals, metrics, thresholds, out) -> str:
    plt = _setup()
    from matplotlib.ticker import PercentFormatter

    fig, axes = plt.subplots(1, len(metrics), figsize=(6.6 * len(metrics), 4.8), layout="constrained", squeeze=False)
    lines, headline = [], {}
    for ax, metric in zip(axes[0], metrics):
        thr = thresholds[metric]
        observed = pairs[metric].dropna()
        if metric == "dcj":
            hi = float(min(DCJ_CAP, max(observed.max() if len(observed) else 0, 2 * thr)))
            clipped = bool(len(observed) and observed.max() > hi)
        else:
            hi, clipped = 1.0, False
        ends = []
        for rel, color in SERIES.items():
            n = totals[rel]
            if not n:
                continue
            v = np.sort(np.minimum(pairs.loc[pairs["relation"] == rel, metric].dropna().to_numpy(float), hi))
            share = len(v) / n
            ax.step(np.r_[0, v, hi], np.r_[0, np.arange(1, len(v) + 1) / n, share], where="post", color=color,
                    linewidth=2, label=f"{rel} ({n:,} pairs)")
            within = float((v <= thr).sum() / n)
            ax.plot([thr], [within], "o", markersize=8, color=color, markeredgecolor=SURFACE, markeredgewidth=1.5,
                    zorder=4)
            ends.append((rel, within, share))
            lines.append(f"{NAME[metric]} <= {thr:g}: {within:.0%} of {rel} pairs; not computed {1 - share:.0%} "
                         f"(of {n:,})")
        by_value = sorted(range(len(ends)), key=lambda k: -ends[k][1])
        for rank, k in enumerate(by_value):   # higher share: above-left of its marker; lower: below-right (or above)
            within = ends[k][1]
            upper = rank == 0 or within < 0.08
            ax.annotate(f"{within:.0%}", (thr, within), xytext=(-7, 7) if rank == 0 else (7, 7 if upper else -7),
                        textcoords="offset points", ha="right" if rank == 0 else "left",
                        va="bottom" if upper else "top", fontsize=9, color=INK)
        order = sorted(range(len(ends)), key=lambda k: ends[k][2])
        for rank, k in enumerate(order):
            rel, _, share = ends[k]
            if share < 0.9995:
                upper = rank == len(order) - 1
                ax.annotate(f"{1 - share:.0%} not computed", (hi, share), xytext=(-2, 3 if upper else -3),
                            textcoords="offset points", ha="right", va="bottom" if upper else "top", fontsize=8,
                            color=INK_2)
        headline[metric] = [e[1] for e in ends]
        _threshold(ax, thr, f"Pling threshold {thr:g}")
        ax.set_xlim(-0.02 * hi, hi * 1.02)
        ax.set_ylim(0, 1.08)
        ax.yaxis.set_major_formatter(PercentFormatter(1.0))
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        if clipped:
            ticks = [t for t in ax.get_xticks() if 0 <= t < hi] + [hi]
            ax.set_xticks(ticks, [f"{t:g}" for t in ticks[:-1]] + [f"{hi:g}+"])
        ax.set_xlabel(f"{NAME[metric]} distance")
        ax.set_ylabel("share of pairs with distance ≤ x")
        ax.set_title(f"{NAME[metric]} distance", loc="left", fontsize=10)
        ax.legend(loc="lower right", fontsize=8.5)
        _style(ax, "both")
    d = headline.get("dcj", [])
    head = (f"DCJ-Indel ≤ {thresholds['dcj']:g}: {d[0]:.0%} of pairs in the same subcommunity, {d[1]:.0%} of pairs that "
            "only share a community" if len(d) == 2 else "Pairs in the same subcommunity vs pairs that only share a community")
    fig.suptitle(head + ". Curves end below 100% by the share of pairs never computed", color=INK, fontsize=11, x=0.01,
                 ha="left")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return "\n".join(lines)


def plot_scatter(pairs, thresholds, out, max_points, seed) -> str:
    plt = _setup()
    from matplotlib.ticker import MaxNLocator

    both = pairs.dropna(subset=["containment"])
    n_all = len(both)
    d = both.sample(max_points, random_state=seed) if n_all > max_points else both
    dcj = d["dcj"].to_numpy(float)
    finite = dcj[~np.isnan(dcj)]
    top = float(max(min(DCJ_CAP, finite.max() if len(finite) else 0), 2 * thresholds["dcj"]))
    strip = top + max(3.0, 0.12 * top)
    jitter = np.random.default_rng(seed).uniform(-0.25, 0.25, len(d))
    y = np.where(np.isnan(dcj), strip, np.minimum(dcj, top)) + jitter
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    for rel in (OTHER, SAME):                                   # same subcommunity drawn on top
        k = (d["relation"] == rel).to_numpy()
        if k.any():
            ax.scatter(d["containment"].to_numpy(float)[k], y[k], s=24, color=SERIES[rel], alpha=0.75,
                       edgecolors=SURFACE, linewidths=0.6, label=f"{rel} ({int(k.sum()):,} pairs)", zorder=3)
    ax.axhline((top + strip) / 2, color=GRID, linewidth=1)
    _threshold(ax, thresholds["containment"], f"containment threshold {thresholds['containment']:g}")
    _threshold(ax, thresholds["dcj"], f"DCJ threshold {thresholds['dcj']:g}", vertical=False)
    ticks = [t for t in MaxNLocator(6, integer=True).tick_values(0, top) if 0 <= t <= top]
    labels = [f"{t:g}" for t in ticks]
    if len(finite) and finite.max() > top:
        ticks, labels = ticks + [top], labels + [f"{top:g}+"]
    ax.set_yticks(ticks + [strip], labels + ["not\ncomputed"])
    ax.set_ylim(-1, strip + 1.5)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("containment distance")
    ax.set_ylabel("DCJ-Indel distance")
    rho = both[["dcj", "containment"]].dropna().corr(method="spearman").iloc[0, 1] if n_all else np.nan
    shown = f", {len(d):,} shown at random" if len(d) < n_all else ""
    ax.set_title(f"{n_all:,} within-community pairs{shown}: DCJ-Indel against containment, Spearman ρ = {rho:.2f}\n"
                 "DCJ-Indel is only computed left of the containment threshold; DCJ values jittered by ±0.25",
                 loc="left", fontsize=10)
    ax.legend(loc="center right", fontsize=8.5)
    _style(ax, "both")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return f"scatter: {n_all:,} pairs with containment, Spearman rho {rho:.2f}"


def plot_heatmaps(pairs, members, summaries, level, metrics, thresholds, out, n_max, max_members) -> str:
    sizes = members[level].value_counts().sort_index().sort_values(ascending=False, kind="stable")
    pick = sizes[(sizes >= 3) & (sizes <= max_members)].head(n_max)
    if pick.empty:
        log.warning("no %s with 3-%d members: %s not written", level, max_members, out)
        return ""
    plt = _setup()
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import Rectangle

    has_cont = "containment" in metrics
    cmaps = {"dcj": LinearSegmentedColormap.from_list("blues", BLUES),
             "containment": LinearSegmentedColormap.from_list("oranges", ORANGES)}
    norms = {"dcj": Normalize(0, 5 * thresholds["dcj"]),
             "containment": Normalize(0, min(1.0, 2 * thresholds["containment"]))}
    flagged = set(summaries["subcommunity"].loc[summaries["subcommunity"]["dcj_censored"] > 0, "subcommunity"])
    ncols = min(4, len(pick))
    nrows = math.ceil(len(pick) / ncols)
    gcols = max(ncols, 2)
    fig = plt.figure(figsize=(3.5 * gcols + 0.4, 3.75 * nrows + 1.3), layout="constrained")
    gs = fig.add_gridspec(nrows + 1, gcols, height_ratios=[1] * nrows + [0.06])
    lines = []
    for k, (g, n_expected) in enumerate(pick.items()):
        r, c = divmod(k, ncols)
        ax = fig.add_subplot(gs[r, :] if ncols == 1 else gs[r, c])
        mem = members[members[level] == g].copy()
        mem["_size"] = mem["subcommunity"].map(mem["subcommunity"].value_counts()).fillna(1)
        mem["_sub"] = mem["subcommunity"].fillna("~")
        mem = mem.sort_values(["_size", "_sub", "plasmid"], ascending=[False, True, True])
        names = mem["plasmid"].tolist()
        n = len(names)
        idx = {p: i for i, p in enumerate(names)}
        sel = pairs[pairs[level] == g]
        i, j = sel["plasmid_1"].map(idx), sel["plasmid_2"].map(idx)
        ok = (i.notna() & j.notna()).to_numpy()
        i, j = i.to_numpy()[ok].astype(int), j.to_numpy()[ok].astype(int)
        a, b = np.minimum(i, j), np.maximum(i, j)
        upper, lower = np.full((n, n), np.nan), np.full((n, n), np.nan)
        upper[a, b] = sel["dcj"].to_numpy(float)[ok]
        lower[b, a] = sel["containment" if has_cont else "dcj"].to_numpy(float)[ok]
        ax.pcolormesh(np.ma.masked_invalid(upper), cmap=cmaps["dcj"], norm=norms["dcj"], edgecolors=SURFACE,
                      linewidth=0.5)
        low = "containment" if has_cont else "dcj"
        ax.pcolormesh(np.ma.masked_invalid(lower), cmap=cmaps[low], norm=norms[low], edgecolors=SURFACE, linewidth=0.5)
        off = ~np.eye(n, dtype=bool)
        missing = np.argwhere(off & np.isnan(np.where(np.triu(off), upper, lower)))
        if len(missing):
            ax.add_collection(PatchCollection([Rectangle((cc, rr), 1, 1) for rr, cc in missing],
                                              facecolor=NOT_COMPUTED, edgecolor=SURFACE, linewidth=0.5, hatch="////"))
        subs = mem["subcommunity"].tolist()
        start, notes = 0, []
        for pos in range(1, n + 1):
            if pos == n or subs[pos] != subs[start] or pd.isna(subs[start]):
                length, sub = pos - start, subs[start]
                if level == "community" and isinstance(sub, str):   # a known subcommunity, even of one plasmid
                    ax.add_patch(Rectangle((start, start), length, length, fill=False, edgecolor=INK, linewidth=1.3))
                    if length >= 2 or n <= 12:                  # one-cell labels only where cells are wide enough
                        ax.text(start + length / 2, n + 0.15, "s" + sub.rsplit("_", 1)[-1], ha="center", va="top",
                                fontsize=7, color=INK_2, clip_on=False)   # below the matrix, clear of the title
                if isinstance(sub, str) and sub in flagged:
                    notes.append("s" + sub.rsplit("_", 1)[-1] if level == "community" else "")
                start = pos
        ax.set_xlim(0, n)
        ax.set_ylim(n, 0)
        ax.set_aspect("equal")
        ax.set_xticks([])
        if n <= 24:
            ax.set_yticks(np.arange(n) + 0.5, names, fontsize=5.5)
        else:
            ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        title = f"{g}  (n={n})"
        if n < n_expected:
            title += f", {n_expected - n} without distances"
        if notes:
            title += "\n" + (", ".join(notes) + ": pairs without DCJ" if level == "community"
                             else "pairs without DCJ")
        ax.set_title(title, loc="left", fontsize=8.5)
        lines.append(f"{g}: {n} members, subcommunities with >= 2 members: "
                     f"{sum(1 for s in set(subs) if isinstance(s, str) and subs.count(s) > 1)}, "
                     f"not-computed cells: {len(missing) // 2 if not has_cont else len(missing)}")
    half = gcols // 2
    bars = [("dcj", gs[-1, :half] if has_cont else gs[-1, :],
             "DCJ-Indel distance" + (" (above the diagonal)" if has_cont else ""))]
    if has_cont:
        bars.append(("containment", gs[-1, half:], "containment distance (below the diagonal)"))
    for metric, spec, label in bars:
        cax = fig.add_subplot(spec)
        cb = fig.colorbar(ScalarMappable(norm=norms[metric], cmap=cmaps[metric]), cax=cax, orientation="horizontal",
                          extend="max")
        cb.set_label(label, color=INK_2)
        cb.outline.set_visible(False)
        cax.axvline(thresholds[metric], color=INK, linestyle=(0, (3, 2)), linewidth=1.2)
    plural = "subcommunities" if level == "subcommunity" else "communities"
    what = ("DCJ-Indel above the diagonal, containment below" if has_cont else "DCJ-Indel")
    boxes = ("; boxes: subcommunities (s = subcommunity number; unboxed: subcommunity unknown)"
             if level == "community" else "")
    fig.suptitle(f"Largest {len(pick)} {plural} with 3–{max_members} members: {what}{boxes}.\n"
                 f"Grey hatched: not computed (DCJ: containment above {thresholds['containment']:g}); "
                 "dashed line on a colour bar: Pling threshold", color=INK, fontsize=10.5, x=0.01, ha="left")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return "\n".join(lines)


# ----------------------------------------------------------------------------------------------- driver
def run(dist_dir, outdir=None, typing=None, heatmap_level="community", n_heatmaps=12, max_heatmap_members=40,
        max_scatter=20000, seed=2026) -> dict:
    dist_dir = Path(dist_dir)
    out = Path(outdir) if outdir else dist_dir
    out.mkdir(parents=True, exist_ok=True)
    info, summaries, pairs = load(dist_dir)
    thresholds = info["thresholds"]
    if typing is None:
        recorded = (info.get("files") or {}).get("typing")
        typing = recorded if recorded and Path(recorded).exists() else None
        if typing is None:
            log.info("typing file %s not found: subcommunities are taken from the pairs", recorded)
    members = membership(pairs, typing)
    totals = pair_totals(summaries, pairs)
    metrics = [m for m in METRICS if pairs[m].notna().any()]
    written, report = {}, []

    sizes = {lvl: group_sizes(lvl, summaries, info, members, typing is not None) for lvl in LEVELS}
    report.append(plot_sizes(sizes, info, out / FIGURES["sizes"]))
    written["sizes"] = FIGURES["sizes"]
    if metrics:
        report.append(plot_same_vs_other(pairs, totals, metrics, thresholds, out / FIGURES["same_vs_other"]))
        written["same_vs_other"] = FIGURES["same_vs_other"]
    if set(METRICS) <= set(metrics):
        report.append(plot_scatter(pairs, thresholds, out / FIGURES["scatter"], max_scatter, seed))
        written["scatter"] = FIGURES["scatter"]
    else:
        log.info("no containment distances: %s not written", FIGURES["scatter"])
    if metrics:
        text = plot_heatmaps(pairs, members, summaries, heatmap_level, metrics, thresholds, out / FIGURES["heatmaps"],
                             n_heatmaps, max_heatmap_members)
        if text:
            report.append(text)
            written["heatmaps"] = FIGURES["heatmaps"]
    for block in report:
        for line in block.splitlines():
            log.info(line)
    return {"figures": written, "typing": str(typing) if typing else None, "pair_totals": totals}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dist_dir", required=True, help="output folder of distance_distribution.py")
    p.add_argument("--outdir", help="where to write the figures (default: --dist_dir)")
    p.add_argument("--typing", help="Pling typing.tsv (default: the one recorded in run_info.json, if it exists)")
    p.add_argument("--heatmap_level", choices=LEVELS, default="community")
    p.add_argument("--n_heatmaps", type=int, default=12, help="largest groups to draw (default 12)")
    p.add_argument("--max_heatmap_members", type=int, default=40,
                   help="skip groups larger than this in the heatmaps (default 40)")
    p.add_argument("--max_scatter", type=int, default=20000,
                   help="draw at most this many pairs in the scatter, chosen at random (default 20000)")
    p.add_argument("--seed", type=int, default=2026)
    a = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    res = run(a.dist_dir, a.outdir, a.typing, a.heatmap_level, a.n_heatmaps, a.max_heatmap_members, a.max_scatter,
              a.seed)
    log.info("written: %s", ", ".join(res["figures"].values()))
    return res


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
