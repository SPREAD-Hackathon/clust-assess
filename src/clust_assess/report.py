
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd

from . import __version__
from .seqio import input_name
from .stats import fmt, median_high, pair_lookup, resample, status, within

FINAL_COLS = ["Plasmid A", "Plasmid B", "DCJ", "community_DCJ", "containment", "community_containment"]
VERDICT_TEXT = {
    "tighter_than_db": "inputs are more closely related than the closest database plasmids of their subcommunity "
                       "(consistent with a recent shared origin / outbreak cluster)",
    "looser_than_db": "inputs are more variable than the closest database plasmids of their subcommunity "
                      "(not closely related)",
    "similar": "inputs vary as much as the closest database plasmids of their subcommunity",
    "single_input": "only one input plasmid in this subcommunity",
    "no_baseline": "too few database plasmids in this subcommunity for a baseline",
    "inputs_split": "inputs fall into different Pling subcommunities: they are not one related cluster",
    "no_assignment": "no input plasmid could be assigned to a database subcommunity",
}


def _read_pairs(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["a", "b", "containment", "dcj", "dcj_censored"])
    df["dcj"] = df["dcj"].replace("inf", math.inf).astype(float)
    return df


def build(queries_tsv, call_files, input_pairs_tsv, baseline_files, params: dict, outdir) -> dict:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    queries = pd.read_csv(queries_tsv, sep="\t", keep_default_na=False)
    calls = pd.concat([pd.read_csv(f, sep="\t", keep_default_na=False) for f in call_files], ignore_index=True) \
        if call_files else pd.DataFrame(columns=["query_id", "status", "subcommunity"])
    inp = queries.merge(calls, on="query_id", how="left", suffixes=("_input", ""))
    inp["status"] = inp["status"].where(inp["status_input"] == "ok", inp["status_input"])
    inp["status"] = inp["status"].replace("", "not_processed").fillna("not_processed")
    name = dict(zip(inp["query_id"], inp["input_id"]))
    sub_of = {r.query_id: r.subcommunity for r in inp.itertuples() if r.status == "called"}

    ipairs = _read_pairs(input_pairs_tsv)
    ilook = pair_lookup(ipairs[["a", "b", "containment", "dcj"]].itertuples(index=False, name=None))

    baselines = {}
    for f in baseline_files:
        f = Path(f)
        if f.name.endswith(".member_pairs.tsv"):
            sc = f.name[: -len(".member_pairs.tsv")]
            members = pd.read_csv(f.with_name(f"{sc}.members.tsv"), sep="\t")
            baselines[sc] = (members, _read_pairs(f))
    sub_rows = []
    for sc in sorted(set(sub_of.values())):
        q = sorted(k for k, v in sub_of.items() if v == sc)
        members, mpairs = baselines.get(sc, (pd.DataFrame(columns=["plasmid_id", "selected"]), _read_pairs_empty()))
        used = list(members.loc[members["selected"].astype(bool), "plasmid_id"].astype(str)) if len(members) else []
        mlook = pair_lookup(mpairs[["a", "b", "containment", "dcj"]].itertuples(index=False, name=None))
        obs_dcj = median_high(within(q, ilook, 1))
        obs_cont = median_high(within(q, ilook, 0))
        base_dcj = median_high(within(used, mlook, 1))
        base_cont = median_high(within(used, mlook, 0))
        rs = resample(used, len(q), mlook, obs_dcj, params.get("n_resample", 1000), params.get("seed", 2026))
        st = status(len(q), obs_dcj, base_dcj, rs["p_tighter"], rs["p_looser"], params.get("alpha", 0.05))
        sub_rows.append({"subcommunity": sc, "community": inp.loc[inp["subcommunity"] == sc, "community"].iloc[0],
                         "n_inputs": len(q), "inputs": ";".join(name[x] for x in q),
                         "n_members_in_db": len(members), "n_members_used": len(used),
                         "median_dcj_inputs": obs_dcj, "median_containment_inputs": obs_cont,
                         "community_DCJ": base_dcj, "community_containment": base_cont,
                         "p_tighter": rs["p_tighter"], "p_looser": rs["p_looser"], "n_subsets": rs["n_subsets"],
                         "exact_test": rs["exact"], "status": st})
    subs = pd.DataFrame(sub_rows, columns=["subcommunity", "community", "n_inputs", "inputs", "n_members_in_db",
                                           "n_members_used", "median_dcj_inputs", "median_containment_inputs",
                                           "community_DCJ", "community_containment", "p_tighter", "p_looser",
                                           "n_subsets", "exact_test", "status"])
    sub_info = subs.set_index("subcommunity").to_dict("index") if len(subs) else {}

    rows = []
    for r in ipairs.itertuples():
        sa, sb = sub_of.get(r.a), sub_of.get(r.b)
        same = sa is not None and sa == sb
        base = sub_info.get(sa, {}) if same else {}
        rows.append({"Plasmid A": name.get(r.a, input_name(r.a)), "Plasmid B": name.get(r.b, input_name(r.b)),
                     "DCJ": r.dcj, "community_DCJ": base.get("community_DCJ", math.nan),
                     "containment": r.containment, "community_containment": base.get("community_containment", math.nan),
                     "subcommunity_A": sa or "NA", "subcommunity_B": sb or "NA", "same_subcommunity": same,
                     "dcj_censored": bool(r.dcj_censored), "subcommunity_status": base.get("status", "NA"),
                     "p_tighter": base.get("p_tighter", math.nan)})
    final = pd.DataFrame(rows, columns=FINAL_COLS + ["subcommunity_A", "subcommunity_B", "same_subcommunity",
                                                     "dcj_censored", "subcommunity_status", "p_tighter"])
    _write_semicolon(final[FINAL_COLS], out / "final_output.csv")
    _write_semicolon(final, out / "final_output_extended.csv")

    in_rows = []
    for r in inp.itertuples():
        sc = r.subcommunity if r.status == "called" else None
        others = [x for x, s in sub_of.items() if s == sc and x != r.query_id] if sc else []
        d_others = [ilook.get(tuple(sorted((r.query_id, o))), (math.nan, math.nan))[1] for o in others]
        info = sub_info.get(sc, {}) if sc else {}
        in_rows.append({"input_id": r.input_id, "query_id": r.query_id, "length": r.length, "status": r.status,
                        "reason": getattr(r, "reason", "") or getattr(r, "note", ""),
                        "ref_plasmid": r.ref_plasmid if sc else "", "subcommunity": sc or "",
                        "community": r.community if sc else "", "accio_score": getattr(r, "final_score", ""),
                        "mash_score": getattr(r, "mash_score", ""), "n_inputs_in_subcommunity": len(others) + 1 if sc else 0,
                        "median_dcj_to_other_inputs": median_high(d_others),
                        "community_DCJ": info.get("community_DCJ", math.nan),
                        "subcommunity_status": info.get("status", "")})
    inputs = pd.DataFrame(in_rows)
    _write_tsv(inputs, out / "inputs.tsv")
    _write_tsv(subs, out / "subcommunities.tsv")

    called = [r for r in in_rows if r["status"] == "called"]
    n_sub = len(sub_info)
    if not called:
        verdict = "no_assignment"
    elif n_sub > 1:
        verdict = "inputs_split"
    else:
        verdict = next(iter(sub_info.values()))["status"]
    summary = {
        "clust_assess_version": __version__, "verdict": verdict, "verdict_text": VERDICT_TEXT.get(verdict, verdict),
        "n_inputs": int(len(inp)), "n_assigned": len(called), "n_subcommunities": n_sub,
        "params": params, "subcommunities": _jsonable(subs), "inputs": _jsonable(inputs),
        "pairs": _jsonable(final),
    }
    (out / "clust_assess.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def _read_pairs_empty() -> pd.DataFrame:
    return pd.DataFrame(columns=["a", "b", "containment", "dcj", "dcj_censored"])


def _write_semicolon(df: pd.DataFrame, path: Path) -> None:
    lines = [";".join(df.columns)]
    for row in df.itertuples(index=False):
        lines.append(";".join(fmt(v) if isinstance(v, float) else str(v) for v in row))
    path.write_text("\n".join(lines) + "\n")


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    lines = ["\t".join(df.columns)]
    for row in df.itertuples(index=False):
        lines.append("\t".join(fmt(v) if isinstance(v, float) else str(v) for v in row))
    path.write_text("\n".join(lines) + "\n")


def _jsonable(df: pd.DataFrame) -> list[dict]:
    recs = df.to_dict("records")
    for r in recs:
        for k, v in r.items():
            if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                r[k] = None
            elif hasattr(v, "item"):
                r[k] = v.item()
    return recs
