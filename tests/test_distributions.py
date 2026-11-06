import math
import random
import statistics
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))   # not bin/: bin/clust_assess.py would shadow the clust_assess package

from clust_assess import distributions as dist  # noqa: E402


def _tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


@pytest.fixture
def pling_dir(tmp_path):
    """community_0: subcommunity_0 = p1 p2 p3, subcommunity_1 = p4 p5, hub h1; community_1: p6 alone."""
    d = tmp_path / "pling_out"
    _tsv(d / "dcj_thresh_4_graph/objects/typing.tsv", ["plasmid", "type"],
         [("p1", "community_0_subcommunity_0"), ("p2", "community_0_subcommunity_0"),
          ("p3", "community_0_subcommunity_0"), ("p4", "community_0_subcommunity_1"),
          ("p5", "community_0_subcommunity_1"), ("p6", "community_1_subcommunity_0")])
    (d / "dcj_thresh_4_graph/objects/hub_plasmids.csv").write_text("hub_plasmids\nh1\n")
    _tsv(d / "containment/containment_communities/objects/communities.tsv", ["plasmid", "community"],
         [(p, "community_0") for p in ("p1", "p2", "p3", "p4", "p5", "h1")] + [("p6", "community_1")])
    # p2-p3 failed containment (no DCJ); p3-p1 is written in reverse order; p1-p6 crosses communities
    _tsv(d / "all_plasmids_distances.tsv", ["plasmid_1", "plasmid_2", "distance"],
         [("p1", "p2", 1), ("p3", "p1", 3), ("p4", "p5", 2), ("p1", "p4", 7), ("h1", "p1", 0), ("p1", "p6", 5)])
    _tsv(d / "containment/all_pairs_containment_distance.tsv", ["plasmid_1", "plasmid_2", "distance"],
         [("p1", "p2", 0.1), ("p1", "p3", 0.2), ("p2", "p3", 0.6), ("p4", "p5", 0.0), ("p1", "p4", 0.3),
          ("h1", "p1", 0.0), ("p1", "p6", 0.9)])
    return d


def test_quantile_higher_matches_median_high():
    rnd = random.Random(1)
    for _ in range(200):
        obs = sorted(rnd.randint(0, 9) for _ in range(rnd.randint(0, 8)))
        n = len(obs) + rnd.randint(0, 5)
        full = obs + [math.inf] * (n - len(obs))
        expect = statistics.median_high(full) if full else math.nan
        got = dist.quantile_higher(np.array(obs, float), n, 0.5)
        assert (math.isnan(got) and math.isnan(expect)) or got == expect


def test_membership_and_pairs(pling_dir):
    files = dist.locate(pling_dir)
    m = dist.membership(files).set_index("plasmid")
    assert m.loc["h1", "hub"] and m.loc["h1", "subcommunity"] is None and m.loc["h1", "community"] is None
    assert m.loc["p4", "community"] == "community_0"
    pairs = dist.within_pairs(files, m.reset_index())
    keys = set(zip(pairs["plasmid_1"], pairs["plasmid_2"]))
    assert ("p1", "p3") in keys and ("p3", "p1") not in keys          # put in canonical order
    assert ("p1", "p6") not in keys and ("h1", "p1") not in keys      # other community / hub
    row = pairs.set_index(["plasmid_1", "plasmid_2"]).loc[("p2", "p3")]
    assert math.isnan(row["dcj"]) and row["containment"] == 0.6 and row["subcommunity"] == "community_0_subcommunity_0"
    assert pd.isna(pairs.set_index(["plasmid_1", "plasmid_2"]).loc[("p1", "p4"), "subcommunity"])


def test_summaries_count_missing_pairs_as_censored(pling_dir, tmp_path):
    info = dist.run(tmp_path / "out", pling_dir)
    sub = pd.read_csv(tmp_path / "out" / "subcommunity_summary.tsv", sep="\t").set_index("subcommunity")
    s0 = sub.loc["community_0_subcommunity_0"]
    assert s0["n_members"] == 3 and s0["n_pairs"] == 3
    assert s0["dcj_n"] == 2 and s0["dcj_censored"] == 1                 # DCJ for p2-p3 was never computed
    assert s0["dcj_min"] == 1 and s0["dcj_median"] == 3 and s0["dcj_max"] == math.inf
    assert s0["containment_censored"] == 0 and s0["containment_median"] == 0.2
    assert s0["dcj_frac_le_threshold"] == pytest.approx(2 / 3)
    assert sub.loc["community_0_subcommunity_1", "dcj_median"] == 2
    assert "community_1_subcommunity_0" not in sub.index                # one member: no pairs
    com = pd.read_csv(tmp_path / "out" / "community_summary.tsv", sep="\t").set_index("community")
    c0 = com.loc["community_0"]
    assert c0["n_members"] == 5 and c0["n_pairs"] == 10 and c0["dcj_n"] == 4 and c0["dcj_censored"] == 6
    assert c0["dcj_median"] == math.inf                                 # most pairs fail containment
    assert info["n_hubs"] == 1 and (tmp_path / "out" / "within_pairs.tsv.gz").exists()


def test_include_hubs_and_flat_layout(pling_dir, tmp_path):
    flat = tmp_path / "flat"
    flat.mkdir()
    for f in ["dcj_thresh_4_graph/objects/typing.tsv", "dcj_thresh_4_graph/objects/hub_plasmids.csv",
              "all_plasmids_distances.tsv"]:
        (flat / Path(f).name).write_text((pling_dir / f).read_text())
    dist.run(tmp_path / "o", flat, include_hubs=True)                   # no communities / containment files
    com = pd.read_csv(tmp_path / "o" / "community_summary.tsv", sep="\t").set_index("community")
    assert com.loc["community_0", "n_members"] == 5                     # hub has no community without communities.tsv
    assert com["containment_n"].isna().all()
    dist.run(tmp_path / "h", pling_dir, include_hubs=True)
    com = pd.read_csv(tmp_path / "h" / "community_summary.tsv", sep="\t").set_index("community")
    assert com.loc["community_0", "n_members"] == 6 and com.loc["community_0", "dcj_n"] == 5


def test_plots_without_containment(pling_dir, tmp_path):
    pytest.importorskip("matplotlib")
    flat = tmp_path / "flat"
    flat.mkdir()
    for f in ["dcj_thresh_4_graph/objects/typing.tsv", "all_plasmids_distances.tsv"]:
        (flat / Path(f).name).write_text((pling_dir / f).read_text())
    info = dist.run(tmp_path / "o", flat, plots=True)            # one-column figures, DCJ only
    assert info["levels"]["subcommunity"]["figure"] == "subcommunity_distribution.png"
    assert (tmp_path / "o" / "community_distribution.png").stat().st_size > 0
    empty = dist.summarise(pd.DataFrame(columns=["plasmid_1", "plasmid_2", "community", "subcommunity", "dcj",
                                                 "containment"]),
                           pd.DataFrame(columns=["plasmid", "subcommunity", "community", "hub"]),
                           "subcommunity", {"dcj": 4, "containment": 0.5})
    assert not dist.plot(pd.DataFrame(columns=["subcommunity", "dcj", "containment"]), empty, "subcommunity",
                         {"dcj": 4, "containment": 0.5}, tmp_path / "none.png")
    assert not (tmp_path / "none.png").exists()


def _load_cli(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_cli_with_plots(pling_dir, tmp_path):
    pytest.importorskip("matplotlib")
    distance_distribution = _load_cli("distance_distribution")
    distance_distribution.main(["--pling_dir", str(pling_dir), "--outdir", str(tmp_path / "cli"), "--plots"])
    for f in ["subcommunity_distribution.png", "community_distribution.png", "run_info.json"]:
        assert (tmp_path / "cli" / f).stat().st_size > 0


def test_insights_script(pling_dir, tmp_path):
    pytest.importorskip("matplotlib")
    ins = _load_cli("distribution_insights")
    dist.run(tmp_path / "d", pling_dir)
    res = ins.main(["--dist_dir", str(tmp_path / "d")])          # typing.tsv recorded in run_info.json
    assert res["typing"] and set(res["figures"]) == {"sizes", "same_vs_other", "scatter", "heatmaps"}
    assert res["pair_totals"] == {ins.SAME: 4, ins.OTHER: 6}      # community_0: 5 members, 10 pairs, 4 within subs
    for f in res["figures"].values():
        assert (tmp_path / "d" / f).stat().st_size > 0
    flat = tmp_path / "flat"                                     # DCJ only, then no typing file either
    flat.mkdir()
    for f in ["dcj_thresh_4_graph/objects/typing.tsv", "all_plasmids_distances.tsv"]:
        (flat / Path(f).name).write_text((pling_dir / f).read_text())
    dist.run(tmp_path / "f", flat)
    (flat / "typing.tsv").unlink()
    res = ins.main(["--dist_dir", str(tmp_path / "f"), "--outdir", str(tmp_path / "f_out")])
    assert res["typing"] is None and set(res["figures"]) == {"sizes", "same_vs_other", "heatmaps"}
    assert res["pair_totals"] == {ins.SAME: 4, ins.OTHER: 6}      # missing pairs counted from the summaries
