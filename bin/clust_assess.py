#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import pandas as pd  

from clust_assess import accio_io, db, pling_io, prep, report  # noqa: E402


def cmd_prep_queries(a):
    t = prep.prep_queries(a.inputs, a.outdir, circular=not a.no_circular)
    t.to_csv(a.table, sep="\t", index=False)


def cmd_prep_refs(a):
    t = prep.prep_refs(a.inputs, a.out, circular=not a.no_circular)
    t.to_csv(a.table, sep="\t", index=False)


def cmd_db_manifest(a):
    print(json.dumps(db.write_manifest(a.db)))


def cmd_verify_db(a):
    Path(a.out).write_text(json.dumps(db.verify(a.db), indent=2) + "\n")


def cmd_parse_accio(a):
    row = accio_io.parse_call(a.query_id, a.dir, a.exit_status)
    pd.DataFrame([row], columns=accio_io.CALL_COLS).to_csv(a.out, sep="\t", index=False)


def cmd_candidates(a):
    db.candidates(a.db, a.subcommunity, a.mode, a.out).to_csv(a.table, sep="\t", index=False)


def cmd_select_members(a):
    cands = pd.read_csv(a.candidates, sep="\t", dtype={"plasmid_id": str, "representative": str})
    m = db.select_members(cands, a.mash, a.queries, a.k)
    m.to_csv(a.table, sep="\t", index=False)
    db.extract(a.candidates_fasta, list(m.loc[m["selected"], "plasmid_id"]), a.outdir)


def cmd_pairs(a):
    ids = [Path(p).name.rsplit(".fasta", 1)[0] for p in Path(a.genomes).read_text().split()]
    if a.pling_dir and Path(a.pling_dir, pling_io.CONTAINMENT_FILE).exists():
        df = pling_io.pairs(a.pling_dir, ids)
    else:
        df = pd.DataFrame(columns=["a", "b", "containment", "dcj", "dcj_censored"])
    df.to_csv(a.out, sep="\t", index=False)


def cmd_report(a):
    params = json.loads(Path(a.params).read_text()) if a.params else {}
    s = report.build(a.queries, a.calls, a.input_pairs, a.baselines or [], params, a.outdir)
    print(f"verdict: {s['verdict']} ({s['n_assigned']}/{s['n_inputs']} inputs assigned, "
          f"{s['n_subcommunities']} subcommunities)")


def main(argv=None):
    p = argparse.ArgumentParser(prog="clust_assess.py")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("prep-queries")
    s.add_argument("--in", dest="inputs", nargs="+", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--table", required=True)
    s.add_argument("--no-circular", action="store_true")
    s.set_defaults(fn=cmd_prep_queries)

    s = sub.add_parser("prep-refs")
    s.add_argument("--in", dest="inputs", nargs="+", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--table", required=True)
    s.add_argument("--no-circular", action="store_true")
    s.set_defaults(fn=cmd_prep_refs)

    s = sub.add_parser("db-manifest")
    s.add_argument("--db", required=True)
    s.set_defaults(fn=cmd_db_manifest)

    s = sub.add_parser("verify-db")
    s.add_argument("--db", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_verify_db)

    s = sub.add_parser("parse-accio")
    s.add_argument("--query-id", required=True)
    s.add_argument("--dir", required=True)
    s.add_argument("--exit-status", type=int, default=0)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_parse_accio)

    s = sub.add_parser("candidates")
    s.add_argument("--db", required=True)
    s.add_argument("--subcommunity", required=True)
    s.add_argument("--mode", choices=["all", "representatives"], default="all")
    s.add_argument("--out", required=True)
    s.add_argument("--table", required=True)
    s.set_defaults(fn=cmd_candidates)

    s = sub.add_parser("select-members")
    s.add_argument("--candidates", required=True)
    s.add_argument("--candidates-fasta", required=True)
    s.add_argument("--mash", required=True)
    s.add_argument("--queries", nargs="+", required=True)
    s.add_argument("--k", type=int, default=20)
    s.add_argument("--outdir", required=True)
    s.add_argument("--table", required=True)
    s.set_defaults(fn=cmd_select_members)

    s = sub.add_parser("pairs")
    s.add_argument("--pling-dir")
    s.add_argument("--genomes", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_pairs)

    s = sub.add_parser("report")
    s.add_argument("--queries", required=True)
    s.add_argument("--calls", nargs="*", default=[])
    s.add_argument("--input-pairs", required=True)
    s.add_argument("--baselines", nargs="*", default=[])
    s.add_argument("--params")
    s.add_argument("--outdir", required=True)
    s.set_defaults(fn=cmd_report)

    a = p.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
