#!/usr/bin/env python3
"""Within-subcommunity and within-community (containment) distance distributions of a Pling run.

Example (PLSDB Pling run on CLIMB):
    bin/distance_distribution.py --pling_dir /shared/team/pling_plsdb --outdir dist_plsdb --plots

Outputs in --outdir:
    subcommunity_summary.tsv   one row per subcommunity: size, pairs, DCJ and containment quantiles, censored pairs
    community_summary.tsv      the same per containment community
    within_pairs.tsv.gz        every computed within-group pair (plasmid_1, plasmid_2, community, subcommunity,
                               dcj, containment); pairs Pling did not compute are not listed, they are counted
                               as *_censored in the summaries
    *_distribution.png         histograms + box plots of the largest groups (with --plots)
    run_info.json              files used, settings and counts
"""
import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from clust_assess import distributions  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pling_dir", help="Pling output folder (Pling's own layout, or a flat copy of its files)")
    p.add_argument("--outdir", required=True)
    for key in distributions.FILES:
        p.add_argument(f"--{key}", help=f"path of the {key} file (overrides the search in --pling_dir)")
    p.add_argument("--include_hubs", action="store_true",
                   help="keep hub plasmids in their containment community (they never have a subcommunity)")
    p.add_argument("--dcj_threshold", type=float, default=4, help="for the fraction of pairs at or below it [4]")
    p.add_argument("--containment_threshold", type=float, default=0.5, help="as above, for containment [0.5]")
    p.add_argument("--min_members", type=int, default=2, help="smallest group size reported [2]")
    p.add_argument("--plots", action="store_true", help="write PNG figures (needs matplotlib)")
    p.add_argument("--top", type=int, default=30, help="largest groups shown as box plots [30]")
    p.add_argument("--chunksize", type=int, default=2_000_000, help="rows read at a time from distance files")
    a = p.parse_args(argv)
    if not a.pling_dir and not (a.typing and a.dcj):
        p.error("give --pling_dir, or at least --typing and --dcj")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    override = {k: getattr(a, k) for k in distributions.FILES}
    info = distributions.run(a.outdir, a.pling_dir, include_hubs=a.include_hubs, dcj_threshold=a.dcj_threshold,
                             containment_threshold=a.containment_threshold, min_members=a.min_members,
                             plots=a.plots, top=a.top, chunksize=a.chunksize, **override)
    for level, d in info["levels"].items():
        logging.info("%s: %s", level, d)
    logging.info("written to %s", Path(a.outdir).resolve())


if __name__ == "__main__":
    main()
