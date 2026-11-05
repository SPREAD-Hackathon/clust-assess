
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .seqio import header, query_id, read_fasta, sanitise, write_fasta

ACCIO_MIN_LENGTH = 10_000  


def prep_queries(inputs: list[str], outdir: str, circular: bool = True) -> pd.DataFrame:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    rows, seen = [], {}
    for src in inputs:
        for name, desc, seq in read_fasta(src):
            qid = query_id(name)
            if qid in seen:
                seen[qid] += 1
                qid = f"{qid}_{seen[qid]}"
            else:
                seen[qid] = 1
            status, note = "ok", ""
            if not seq:
                status, note = "empty", "no sequence"
            elif len(seq) < ACCIO_MIN_LENGTH:
                note = f"shorter than {ACCIO_MIN_LENGTH} bp: ACCIO databases exclude such plasmids by default"
            if status == "ok":
                write_fasta(out / f"{qid}.fasta", [(header(qid, circular), seq)])
            rows.append({"query_id": qid, "input_id": name, "description": desc, "source_file": Path(src).name,
                         "length": len(seq), "status": status, "note": note})
    df = pd.DataFrame(rows, columns=["query_id", "input_id", "description", "source_file", "length", "status", "note"])
    if not (df["status"] == "ok").any():
        raise SystemExit("no usable input plasmids")
    return df


def prep_refs(inputs: list[str], out_fasta: str, circular: bool = True) -> pd.DataFrame:
    rows, recs, seen = [], [], set()
    for src in inputs:
        for name, desc, seq in read_fasta(src):
            rid = sanitise(name)
            if rid in seen:
                raise SystemExit(f"duplicate reference id after sanitising: {rid}")
            if not seq:
                continue
            seen.add(rid)
            recs.append((header(rid, circular), seq))
            rows.append({"ref_id": rid, "input_id": name, "source_file": Path(src).name, "length": len(seq)})
    if not recs:
        raise SystemExit("no reference sequences")
    write_fasta(out_fasta, recs)
    return pd.DataFrame(rows)
