#!/usr/bin/env python3

import argparse
import itertools
import random
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


def read(path, want):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif name in want:
                seqs.setdefault(name, []).append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def write(path, recs):
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def mash_dist(seqs, mash):
    with tempfile.TemporaryDirectory() as d:
        fa = Path(d) / "x.fasta"
        write(fa, seqs.items())
        subprocess.run([mash, "sketch", "-i", "-o", f"{d}/x", str(fa)], check=True, capture_output=True)
        out = subprocess.run([mash, "dist", f"{d}/x.msh", f"{d}/x.msh"], check=True, capture_output=True, text=True)
    return {(a, b): float(dd) for a, b, dd, *_ in (l.split("\t") for l in out.stdout.splitlines())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool-fasta", required=True)
    ap.add_argument("--pool-meta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mash", default="mash")
    ap.add_argument("--f1", default="AA038")
    ap.add_argument("--f2", default="AA002")
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    meta = pd.read_csv(a.pool_meta, sep="\t", dtype=str, keep_default_na=False)
    meta["length"] = meta["length"].astype(int)
    ok = meta[(meta["length"] >= 10_000) & (meta["length"] <= 120_000) & (meta["replicons"] != "")]

    f1 = sorted(ok.loc[ok["primary_cluster"] == a.f1, "accession"])
    f2 = sorted(ok.loc[ok["primary_cluster"] == a.f2, "accession"])
    seqs = read(a.pool_fasta, set(f1) | set(f2) | set(meta["accession"]))
    d = mash_dist({k: seqs[k] for k in f1}, a.mash)
    trio = min(itertools.combinations(f1, 3),
               key=lambda t: (max(d[(x, y)] for x, y in itertools.combinations(t, 2)), t))
    near = sorted((x for x in f1 if x not in trio), key=lambda x: (min(d[(x, t)] for t in trio), x))
    f1_refs = near[:15]
    f2_query = rng.choice(f2)
    f2_refs = rng.sample([x for x in f2 if x != f2_query], min(10, len(f2) - 1))
    others = rng.sample(sorted(set(ok["accession"]) - set(f1) - set(f2)), 8)
    small = sorted(meta.loc[(meta["length"] < 10_000) & (meta["replicons"] != ""), "accession"])[0]
    refs = f1_refs + f2_refs + others

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    write(out / "refs.fasta", [(r, seqs[r]) for r in refs])
    inputs = [(f"outbreak_{t}", seqs[t]) for t in trio] + [(f"other_{f2_query}", seqs[f2_query]),
                                                          (f"dbmember_{f1_refs[0]}", seqs[f1_refs[0]]),
                                                          (f"small_{small}", seqs[small])]
    write(out / "inputs.fasta", inputs)
    write(out / "inputs_outbreak.fasta", inputs[:3])
    trio_max = max(d[(x, y)] for x, y in itertools.combinations(trio, 2))
    print(f"refs={len(refs)} (F1 {len(f1_refs)}, F2 {len(f2_refs)}, other {len(others)}); "
          f"trio={trio} max mash {trio_max:.5f}; inputs={len(inputs)}")


if __name__ == "__main__":
    main()
