
from __future__ import annotations

import math
from itertools import combinations
from pathlib import Path

import pandas as pd

INF = math.inf
DCJ_FILES = ["all_plasmids_distances.tsv", "all_plasmid_distances.tsv"]
CONTAINMENT_FILE = "containment/all_pairs_containment_distance.tsv"


def _key(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def _read(path: Path) -> dict[tuple[str, str], float]:
    df = pd.read_csv(path, sep="\t")
    return {_key(str(r.plasmid_1), str(r.plasmid_2)): float(r.distance) for r in df.itertuples()}


def locate(pling_dir: str | Path) -> tuple[Path, Path | None]:
    d = Path(pling_dir)
    cont = d / CONTAINMENT_FILE
    if not cont.exists():
        raise FileNotFoundError(f"no containment table in {d}")
    dcj = next((d / f for f in DCJ_FILES if (d / f).exists()), None)
    return cont, dcj


QUERY_DB_COLS = ["query_id", "db_plasmid", "containment", "dcj", "dcj_censored"]


def split_query_pairs(df: pd.DataFrame, is_query) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split the pairs of one Pling run on database members + inputs into (member-member pairs,
    input-member pairs with the input first). Input-input pairs are dropped: PLING_INPUTS has them."""
    qa, qb = df["a"].map(is_query).astype(bool), df["b"].map(is_query).astype(bool)
    mixed = df[qa != qb]
    flip = qb[qa != qb]
    qd = pd.DataFrame({"query_id": mixed["b"].where(flip, mixed["a"]), "db_plasmid": mixed["a"].where(flip, mixed["b"]),
                       "containment": mixed["containment"], "dcj": mixed["dcj"], "dcj_censored": mixed["dcj_censored"]},
                      columns=QUERY_DB_COLS)
    return df[~qa & ~qb].reset_index(drop=True), qd.sort_values(["query_id", "db_plasmid"]).reset_index(drop=True)


def pairs(pling_dir: str | Path, ids: list[str]) -> pd.DataFrame:
    cont_path, dcj_path = locate(pling_dir)
    cont = _read(cont_path)
    dcj = _read(dcj_path) if dcj_path else {}
    rows = []
    for a, b in combinations(sorted(set(ids)), 2):
        k = _key(a, b)
        c = cont.get(k, 1.0)
        d = dcj.get(k)
        rows.append({"a": a, "b": b, "containment": c, "dcj": INF if d is None else d, "dcj_censored": d is None})
    return pd.DataFrame(rows, columns=["a", "b", "containment", "dcj", "dcj_censored"])
