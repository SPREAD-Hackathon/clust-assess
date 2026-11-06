
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .seqio import input_name, read_fasta, write_fasta

REQUIRED = ["plasmidDB.fasta", "plasmidDB_info.csv", "plasmidDB.fasta.msh"]
MANIFEST_FILES = REQUIRED + ["pling/dcj_thresh_4_graph/objects/typing.tsv",
                             "pling/containment/containment_communities/objects/communities.tsv"]


def load_info(db: str | Path) -> pd.DataFrame:
    df = pd.read_csv(Path(db) / "plasmidDB_info.csv", low_memory=False)
    df = df.loc[:, ~df.columns.duplicated()]
    keep = [c for c in ["plasmid_id", "length", "pling_type", "community", "rep_types"] if c in df.columns]
    df = df[keep].copy()
    df["plasmid_id"] = df["plasmid_id"].astype(str)
    return df


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(db: str | Path) -> dict:
    d = Path(db)
    lines = [f"{sha256(d / f)}  {f}" for f in MANIFEST_FILES if (d / f).exists()]
    (d / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    info = load_info(d)
    summary = {"db_name": d.name, "n_plasmids": int(len(info)),
               "n_subcommunities": int(info["pling_type"].nunique()),
               "n_communities": int(info["community"].nunique()) if "community" in info else None,
               "manifest_sha256": sha256(d / "MANIFEST.sha256")}
    (d / "db_info.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def verify(db: str | Path) -> dict:
    d = Path(db)
    missing = [f for f in REQUIRED if not (d / f).exists()]
    if missing:
        raise SystemExit(f"database {d} is missing {missing}")
    status = "no_manifest"
    man = d / "MANIFEST.sha256"
    if man.exists():
        for line in man.read_text().splitlines():
            digest, name = line.split(None, 1)
            if sha256(d / name.strip()) != digest:
                raise SystemExit(f"database file changed since build: {name}")
        status = "verified"
    info = load_info(d)
    return {"db_name": d.name, "n_plasmids": int(len(info)), "n_subcommunities": int(info["pling_type"].nunique()),
            "manifest": status, "manifest_sha256": sha256(man) if man.exists() else None}


def candidates(db: str | Path, subcommunity: str, mode: str, out_fasta: str | Path) -> pd.DataFrame:
    """Database plasmids of `subcommunity`. ACCIO de-duplicates near-identical plasmids (skani) and keeps
    one representative per cluster; mode 'all' adds the de-duplicated copies back from
    plasmids_before_clustering.fasta so the baseline reflects the database, not only representatives."""
    d = Path(db)
    info = load_info(d)
    reps = info[info["pling_type"] == subcommunity][["plasmid_id", "length"]].copy()
    reps["representative"] = reps["plasmid_id"]
    rep_map = d / "representative_plasmids.csv"
    source = d / "plasmidDB.fasta"
    if mode == "all" and rep_map.exists() and (d / "plasmids_before_clustering.fasta").exists():
        m = pd.read_csv(rep_map)
        m = m[m["Representative"].astype(str).isin(reps["plasmid_id"])]
        cands = pd.DataFrame({"plasmid_id": m["query"].astype(str), "length": m["len"],
                              "representative": m["Representative"].astype(str)})
        cands = pd.concat([reps, cands]).drop_duplicates("plasmid_id")
        source = d / "plasmids_before_clustering.fasta"
    else:
        cands = reps
    cands = cands.sort_values("plasmid_id").reset_index(drop=True)
    write_fasta(out_fasta, [(n, s) for n, _, s in read_fasta(source) if n in set(cands["plasmid_id"])])
    return cands


def _fasta_stem(name: str) -> str:
    n = Path(name).name
    for ext in (".fasta", ".fa", ".fna"):
        if n.endswith(ext):
            return n[: -len(ext)]
    return n


def _seq_key(seq: str) -> str:
    return hashlib.sha256(seq.upper().encode()).hexdigest()


def input_copies(cand_fasta: str | Path, query_fastas: list[str | Path]) -> dict[str, tuple[str, str]]:
    """Database candidates that are the input plasmids themselves, as {candidate id: (query id, match)}.
    match is 'sequence' (identical sequence, case-insensitive), 'name' (candidate id == input name) or both.
    A rotated or reverse-complemented copy is not recognised; near-identical plasmids are kept on purpose:
    they are genuine relatives, not the input."""
    by_seq, by_name = {}, {}
    for f in query_fastas:
        for qid, _, seq in read_fasta(f):
            by_seq.setdefault(_seq_key(seq), qid)
            by_name[input_name(qid)] = qid
    out = {}
    for cid, _, seq in read_fasta(cand_fasta):
        q_seq, q_name = by_seq.get(_seq_key(seq)), by_name.get(cid)
        if q_seq or q_name:
            match = "name+sequence" if q_seq and q_name and q_seq == q_name else "sequence" if q_seq else "name"
            out[cid] = (q_seq or q_name, match)
    return out


def select_members(cands: pd.DataFrame, mash_tsv: str | Path, query_ids: list[str], k: int,
                   copies: dict[str, tuple[str, str]] | None = None) -> pd.DataFrame:
    """The k candidates closest (min Mash distance) to any input; k <= 0 keeps all. Candidates in `copies`
    (the inputs themselves, see input_copies) are listed but never selected and take no slot of the k."""
    copies = copies or {}
    mash = pd.read_csv(mash_tsv, sep="\t", header=None, names=["ref", "query", "dist", "p", "shared"],
                       dtype={"ref": str, "query": str})
    mash["query"] = mash["query"].map(_fasta_stem)     # mash names a query FASTA by its file: Q__x.fasta -> Q__x
    mash = mash[mash["query"].isin(query_ids)]
    closest = mash.groupby("ref")["dist"].min()
    members = cands.copy()
    members["min_mash_to_inputs"] = members["plasmid_id"].map(closest).fillna(1.0)
    members["copy_of_input"] = members["plasmid_id"].map(lambda p: copies.get(p, ("", ""))[0])
    members["copy_match"] = members["plasmid_id"].map(lambda p: copies.get(p, ("", ""))[1])
    members = members.sort_values(["min_mash_to_inputs", "plasmid_id"]).reset_index(drop=True)
    eligible = members["copy_of_input"] == ""
    members["rank"] = pd.Series(pd.NA, index=members.index, dtype="Int64")
    members.loc[eligible, "rank"] = range(1, int(eligible.sum()) + 1)
    members["selected"] = (eligible & ((members["rank"] <= k).fillna(False) if k > 0 else eligible)).astype(bool)
    return members


def extract(fasta: str | Path, ids: list[str], outdir: str | Path) -> list[Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    want, paths = set(ids), []
    for name, _, seq in read_fasta(fasta):
        if name in want:
            p = out / f"{name}.fasta"
            write_fasta(p, [(name, seq)])
            paths.append(p)
    missing = want - {p.stem for p in paths}
    if missing:
        raise SystemExit(f"members not found in {fasta}: {sorted(missing)[:5]}")
    return paths
