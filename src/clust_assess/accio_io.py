
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

CALL_COLS = ["query_id", "status", "reason", "ref_plasmid", "subcommunity", "community", "is_hub", "n_calls",
             "final_score", "mash_score", "blast_identity", "blast_coverage", "synteny_score", "length_score",
             "ref_length", "ref_rep_type", "all_calls"]
NO_CALL_REASONS = [
    ("Filtered to 0 plasmids with matching replicon types", "no_replicon_match"),
    ("Found 0 plasmids with BLAST matches", "no_blast_match"),
]


def community_of(subcommunity: str) -> str:
    m = re.match(r"^(.*?)_(subcommunity_\d+|hub\d+)$", subcommunity)
    return m.group(1) if m else subcommunity


def _f(x):
    try:
        return round(float(x), 6)
    except (TypeError, ValueError):
        return None


def parse_call(query_id: str, accio_dir: str | Path, exit_status: int = 0) -> dict:
    d = Path(accio_dir)
    row = {c: None for c in CALL_COLS}
    row.update(query_id=query_id, n_calls=0, is_hub=False)
    chosen = sorted(d.glob("*_plasmids_chosen.csv"))
    if exit_status != 0 or not chosen:
        row.update(status="error_accio", reason=f"accio find exit {exit_status}" if exit_status else "no output")
        return row
    df = pd.read_csv(chosen[0])
    calls = df[(df["ACCIO_Plasmids"] == "chosen_plasmid") & df["plasmid_id"].notna()].copy()
    if calls.empty:
        log = "".join(p.read_text(errors="ignore") for p in d.glob("*.log"))
        reason = next((r for pat, r in NO_CALL_REASONS if pat in log), "below_score_thresholds")
        row.update(status="no_assignment", reason=reason)
        return row
    calls["final_score"] = pd.to_numeric(calls.get("final_score"), errors="coerce")
    calls = calls.sort_values("final_score", ascending=False)
    top = calls.iloc[0]
    sub = str(top["pling_type"])
    row.update(status="called", reason="", ref_plasmid=str(top["plasmid_id"]), subcommunity=sub,
               community=community_of(sub), is_hub="_hub" in sub, n_calls=len(calls),
               final_score=_f(top.get("final_score")), mash_score=_f(top.get("mash_score")),
               blast_identity=_f(top.get("blast_identity")), blast_coverage=_f(top.get("blast_coverage")),
               synteny_score=_f(top.get("synteny_score")), length_score=_f(top.get("length_score")),
               ref_length=_f(top.get("length")), ref_rep_type=top.get("rep_type"),
               all_calls=";".join(f"{r.plasmid_id}|{r.pling_type}|{_f(r.final_score)}" for r in calls.itertuples()))
    if pd.isna(top["pling_type"]) or not sub or sub == "nan":
        row.update(status="no_assignment", reason="called plasmid has no Pling subcommunity")
    return row
