import math
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from clust_assess import accio_io, pling_io, prep, report, seqio, stats  # noqa: E402

FIX = ROOT / "tests" / "fixtures"


def test_ids():
    assert seqio.sanitise("NZ_CP044032.1 plasmid pX") == "NZ_CP044032.1_plasmid_pX"
    assert seqio.sanitise("a|b:c") == "a_b_c"
    assert seqio.query_id("NZ_CP1.1") == "Q__NZ_CP1.1"
    assert seqio.input_name("Q__NZ_CP1.1") == "NZ_CP1.1" and seqio.input_name("x") == "x"
    assert seqio.header("Q__a", True) == "Q__a circular=true" and seqio.header("Q__a", False) == "Q__a"


def test_prep_queries_splits_and_dedups(tmp_path):
    fa = tmp_path / "in.fasta"
    fa.write_text(">pA desc\n" + "A" * 12000 + "\n>pA\n" + "C" * 12000 + "\n>small\nACGT\n>empty\n")
    t = prep.prep_queries([str(fa)], tmp_path / "q")
    assert list(t["query_id"]) == ["Q__pA", "Q__pA_2", "Q__small", "Q__empty"]
    assert list(t["status"]) == ["ok", "ok", "ok", "empty"]
    assert "10000" in t.loc[2, "note"]
    assert sorted(p.name for p in (tmp_path / "q").iterdir()) == ["Q__pA.fasta", "Q__pA_2.fasta", "Q__small.fasta"]
    assert (tmp_path / "q" / "Q__pA.fasta").read_text().startswith(">Q__pA circular=true\n")


def test_parse_accio_called():
    r = accio_io.parse_call("Q__q2", FIX / "accio" / "called")
    assert r["status"] == "called" and r["ref_plasmid"] == "NZ_CP093504.1"
    assert r["subcommunity"] == "community_1_subcommunity_1" and r["community"] == "community_1"
    assert r["final_score"] > 0.9 and r["n_calls"] == 1 and not r["is_hub"]


def test_parse_accio_no_call_reason():
    r = accio_io.parse_call("Q__q1", FIX / "accio" / "no_call")
    assert r["status"] == "no_assignment" and r["reason"] == "no_replicon_match"


def test_parse_accio_error(tmp_path):
    assert accio_io.parse_call("Q__x", tmp_path, exit_status=1)["status"] == "error_accio"


def test_community_of():
    assert accio_io.community_of("community_12_subcommunity_3") == "community_12"
    assert accio_io.community_of("community_2_hub1") == "community_2"


def test_pling_pairs_censoring():
    ids = ["NZ_CP044032.1", "NZ_CP133885.1", "CP117218.1"]
    df = pling_io.pairs(FIX / "pling", ids).set_index(["a", "b"])
    assert len(df) == 3
    ok = df.loc[("NZ_CP044032.1", "NZ_CP133885.1")]
    assert ok["dcj"] == 2 and not ok["dcj_censored"]
    cens = df.loc[("CP117218.1", "NZ_CP044032.1")]
    assert math.isinf(cens["dcj"]) and cens["dcj_censored"] and cens["containment"] > 0.5


def test_median_high_and_censoring():
    assert stats.median_high([1, 2, 3, 4]) == 3
    assert stats.median_high([0, math.inf, math.inf]) == math.inf
    assert math.isnan(stats.median_high([]))


def test_resample_exact_and_status():
    members = ["m1", "m2", "m3", "m4"]
    look = stats.pair_lookup([("m1", "m2", .1, 2), ("m1", "m3", .1, 6), ("m1", "m4", .1, 8),
                              ("m2", "m3", .1, 7), ("m2", "m4", .1, 9), ("m3", "m4", .1, 10)])
    rs = stats.resample(members, 2, look, observed=2, n_draws=1000)
    assert rs["exact"] and rs["n_subsets"] == 6 and rs["p_tighter"] == pytest.approx(1 / 6)
    rs3 = stats.resample(members, 3, look, observed=0)
    assert rs3["p_tighter"] == 0 and rs3["n_subsets"] == 4
    assert stats.status(3, 0, 8, 0.0, 1.0) == "tighter_than_db"
    assert stats.status(3, 12, 8, 1.0, 0.0) == "looser_than_db"
    assert stats.status(3, 8, 8, 0.5, 0.5) == "similar"
    assert stats.status(1, math.nan, 8, math.nan, math.nan) == "single_input"
    assert stats.status(2, 1, math.nan, math.nan, math.nan) == "no_baseline"


def test_fmt():
    assert stats.fmt(3.0) == "3" and stats.fmt(0.25) == "0.25" and stats.fmt(math.inf) == "NA"


def _write(path, df):
    df.to_csv(path, sep="\t", index=False)
    return path


def test_report_final_output_format(tmp_path):
    q = _write(tmp_path / "queries.tsv", pd.DataFrame({
        "query_id": ["Q__pA", "Q__pB", "Q__pC", "Q__pD"], "input_id": ["pA", "pB", "pC", "pD"],
        "description": "", "source_file": "in.fasta", "length": 50000, "status": "ok", "note": ""}))
    calls = []
    for qid, sc in [("Q__pA", "community_0_subcommunity_0"), ("Q__pB", "community_0_subcommunity_0"),
                    ("Q__pC", "community_0_subcommunity_0"), ("Q__pD", "")]:
        row = {c: "" for c in accio_io.CALL_COLS}
        row.update(query_id=qid, status="called" if sc else "no_assignment", subcommunity=sc,
                   community="community_0" if sc else "", ref_plasmid="m1" if sc else "")
        calls.append(_write(tmp_path / f"{qid}.call.tsv", pd.DataFrame([row])))
    ip = _write(tmp_path / "input_pairs.tsv", pd.DataFrame(
        [("Q__pA", "Q__pB", 0.05, 1, False), ("Q__pA", "Q__pC", 0.5, 5, False), ("Q__pB", "Q__pC", 0.5, 4, False),
         ("Q__pA", "Q__pD", 0.9, math.inf, True), ("Q__pB", "Q__pD", 0.9, math.inf, True),
         ("Q__pC", "Q__pD", 0.9, math.inf, True)], columns=["a", "b", "containment", "dcj", "dcj_censored"]))
    sc = "community_0_subcommunity_0"
    _write(tmp_path / f"{sc}.members.tsv", pd.DataFrame({"plasmid_id": ["m1", "m2", "m3"], "selected": True}))
    mp = _write(tmp_path / f"{sc}.member_pairs.tsv", pd.DataFrame(
        [("m1", "m2", 0.7, 3, False), ("m1", "m3", 0.7, 3, False), ("m2", "m3", 0.6, 5, False)],
        columns=["a", "b", "containment", "dcj", "dcj_censored"]))
    s = report.build(q, calls, ip, [mp, tmp_path / f"{sc}.members.tsv"], {"n_resample": 100}, tmp_path / "out")
    lines = (tmp_path / "out" / "final_output.csv").read_text().splitlines()
    assert lines[0] == "Plasmid A;Plasmid B;DCJ;community_DCJ;containment;community_containment"
    assert "pA;pB;1;3;0.05;0.7" in lines
    assert "pA;pD;NA;NA;0.9;NA" in lines
    assert len(lines) == 7
    sub = pd.read_csv(tmp_path / "out" / "subcommunities.tsv", sep="\t").iloc[0]
    assert sub["median_dcj_inputs"] == 4 and sub["community_DCJ"] == 3 and sub["status"] in {"similar", "looser_than_db"}
    assert s["verdict"] == sub["status"] and s["n_assigned"] == 3


def test_report_json_is_strict_and_uses_null(tmp_path):
    import json

    test_report_final_output_format(tmp_path)
    text = (tmp_path / "out" / "clust_assess.json").read_text()

    def reject(token):
        raise ValueError(f"non-standard JSON token {token}")

    d = json.loads(text, parse_constant=reject)          # NaN / Infinity would break JSON.parse in a browser
    rows = d["inputs"] + d["pairs"] + d["subcommunities"]
    assert not [(k, v) for r in rows for k, v in r.items() if v in ("", "NA")]
    pd_row = next(r for r in d["inputs"] if r["input_id"] == "pD")
    assert pd_row["status"] == "no_assignment" and pd_row["ref_plasmid"] is None and pd_row["subcommunity"] is None
    pair = next(p for p in d["pairs"] if {p["Plasmid A"], p["Plasmid B"]} == {"pA", "pD"})
    assert pair["DCJ"] is None and pair["dcj_censored"] is True and pair["containment"] == 0.9


def test_report_input_to_db_and_copies(tmp_path):
    test_report_final_output_format(tmp_path)                 # queries, calls, input pairs and a baseline in tmp_path
    sc = "community_0_subcommunity_0"
    _write(tmp_path / f"{sc}.members.tsv", pd.DataFrame({
        "plasmid_id": ["m1", "m2", "m3", "pA_copy"], "selected": [True, True, True, False],
        "copy_of_input": ["", "", "", "Q__pA"], "copy_match": ["", "", "", "sequence"]}))
    qd = _write(tmp_path / f"{sc}.input_db_pairs.tsv", pd.DataFrame(
        [("Q__pA", "m1", 0.1, 2, False), ("Q__pA", "m2", 0.05, 2, False), ("Q__pA", "m3", 0.6, math.inf, True),
         ("Q__pB", "m1", 0.7, math.inf, True), ("Q__pB", "m2", 0.4, math.inf, True)], columns=pling_io.QUERY_DB_COLS))
    s = report.build(tmp_path / "queries.tsv", sorted(tmp_path.glob("Q__*.call.tsv")), tmp_path / "input_pairs.tsv",
                     [tmp_path / f"{sc}.member_pairs.tsv", tmp_path / f"{sc}.members.tsv", qd], {"n_resample": 100},
                     tmp_path / "out2")
    inp = pd.read_csv(tmp_path / "out2" / "inputs.tsv", sep="\t").set_index("input_id")
    a = inp.loc["pA"]                                         # DCJ tie (2, 2) broken by containment: m2
    assert a["nearest_db_plasmid"] == "m2" and a["nearest_db_dcj"] == 2 and a["nearest_db_containment"] == 0.05
    assert a["median_dcj_to_db"] == 2 and a["n_db_compared"] == 3 and a["in_db_as"] == "pA_copy"
    b = inp.loc["pB"]                                         # no DCJ to any member: nearest by containment
    assert b["nearest_db_plasmid"] == "m2" and pd.isna(b["nearest_db_dcj"]) and b["nearest_db_containment"] == 0.4
    assert pd.isna(inp.loc["pD", "nearest_db_plasmid"]) and inp.loc["pD", "n_db_compared"] == 0
    sub = pd.read_csv(tmp_path / "out2" / "subcommunities.tsv", sep="\t").iloc[0]
    assert sub["n_members_in_db"] == 4 and sub["n_members_used"] == 3 and sub["n_members_excluded"] == 1
    assert pd.isna(sub["median_dcj_inputs_to_db"])             # [2, 2, inf, inf, inf]: more than half not computed
    pb = next(r for r in s["inputs"] if r["input_id"] == "pB")
    assert pb["nearest_db_dcj"] is None and pb["nearest_db_plasmid"] == "m2"


def test_split_query_pairs():
    df = pd.DataFrame([("Q__a", "m1", 0.1, 2.0, False), ("m1", "m2", 0.2, 3.0, False), ("Q__a", "Q__b", 0.0, 0.0, False),
                       ("m2", "Q__b", 0.6, math.inf, True)], columns=["a", "b", "containment", "dcj", "dcj_censored"])
    members, qd = pling_io.split_query_pairs(df, seqio.is_query)
    assert list(zip(members["a"], members["b"])) == [("m1", "m2")]                   # input-input pair dropped
    assert list(zip(qd["query_id"], qd["db_plasmid"])) == [("Q__a", "m1"), ("Q__b", "m2")]   # input always first
    assert bool(qd.iloc[1]["dcj_censored"]) and math.isinf(qd.iloc[1]["dcj"])
    m0, q0 = pling_io.split_query_pairs(df.iloc[:0], seqio.is_query)
    assert m0.empty and q0.empty and list(q0.columns) == pling_io.QUERY_DB_COLS


def test_input_copies_left_out_of_baseline(tmp_path):
    from clust_assess import db
    cand = tmp_path / "cands.fasta"
    cand.write_text(">R1 circular=true\nACGTACGT\n>NZ_X.1\nTTTT\n>R2\nGGGG\n>R3\nCCCC\n>Q2\nACGTACGT\n")
    qa, qx = tmp_path / "Q__a.fasta", tmp_path / "Q__NZ_X.1.fasta"
    qa.write_text(">Q__a circular=true\nacgtacgt\n")                                 # R1 and Q2: same sequence
    qx.write_text(">Q__NZ_X.1 circular=true\nAAAA\n")                                # NZ_X.1: same name only
    copies = db.input_copies(cand, [qa, qx])
    assert copies == {"R1": ("Q__a", "sequence"), "Q2": ("Q__a", "sequence"), "NZ_X.1": ("Q__NZ_X.1", "name")}
    cands = pd.DataFrame({"plasmid_id": ["R1", "NZ_X.1", "R2", "R3", "Q2"], "length": 8})
    mash = tmp_path / "m.tsv"
    mash.write_text("R1\tQ__a\t0.0\t0\t1/1\nNZ_X.1\tQ__NZ_X.1\t0.0\t0\t1/1\nR2\tQ__a\t0.01\t0\t1/1\n"
                    "R3\tQ__a\t0.02\t0\t1/1\nQ2\tQ__a\t0.0\t0\t1/1\n")
    m = db.select_members(cands, mash, ["Q__a", "Q__NZ_X.1"], k=1, copies=copies).set_index("plasmid_id")
    assert list(m.index[m["selected"]]) == ["R2"]                                      # copies take no slot of k
    assert m.loc["R1", "copy_of_input"] == "Q__a" and pd.isna(m.loc["R1", "rank"]) and m.loc["R2", "rank"] == 1
    every = db.select_members(cands, mash, ["Q__a"], k=0, copies=copies).set_index("plasmid_id")["selected"]
    assert every.to_dict() == {"R1": False, "Q2": False, "R2": True, "R3": True, "NZ_X.1": False}   # k=0: all but copies
    kept = db.select_members(cands, mash, ["Q__a"], k=1)                              # --keep-input-copies
    assert kept["copy_of_input"].eq("").all() and kept.loc[kept["selected"], "plasmid_id"].tolist() == ["Q2"]


def test_select_members_reads_mash_file_names(tmp_path):
    """`mash dist ref.msh Q__a.fasta` names the query by its file; the k closest must still be by distance."""
    from clust_assess import db
    cands = pd.DataFrame({"plasmid_id": ["A1", "B1", "C1"], "length": 8})
    mash = tmp_path / "m.tsv"
    mash.write_text("A1\tQ__a.1.fasta\t0.09\t0\t1/1\nB1\tQ__a.1.fasta\t0.001\t0\t1/1\n"
                    "C1\tdir/Q__b.fasta\t0.02\t0\t1/1\nA1\tQ__b.fasta\t0.05\t0\t1/1\n")
    m = db.select_members(cands, mash, ["Q__a.1", "Q__b"], k=2).set_index("plasmid_id")
    assert m["min_mash_to_inputs"].to_dict() == {"B1": 0.001, "C1": 0.02, "A1": 0.05}
    assert list(m.index[m["selected"]]) == ["B1", "C1"]                              # not alphabetical A1, B1


def test_candidates_include_deduplicated_copies(tmp_path):
    from clust_assess import db
    d = tmp_path / "db"
    d.mkdir()
    (d / "plasmidDB_info.csv").write_text(
        "plasmid_id,length,pling_type,community\nR1,40000,community_1_subcommunity_0,community_1\n"
        "R2,45000,community_1_subcommunity_0,community_1\nX1,60000,community_0_subcommunity_0,community_0\n")
    (d / "representative_plasmids.csv").write_text(
        "query,cluster,len,Plasmid,Representative\nR1,1,40000,R1,R1\nc1,1,39990,c1,R1\nc2,1,39995,c2,R1\n"
        "R2,2,45000,R2,R2\nX1,3,60000,X1,X1\n")
    seqs = "".join(f">{n} circular=true\nACGT\n" for n in ["R1", "c1", "c2", "R2", "X1"])
    (d / "plasmids_before_clustering.fasta").write_text(seqs)
    (d / "plasmidDB.fasta").write_text("".join(f">{n}\nACGT\n" for n in ["R1", "R2", "X1"]))
    sc = "community_1_subcommunity_0"
    all_ = db.candidates(d, sc, "all", tmp_path / "a.fasta")
    assert sorted(all_["plasmid_id"]) == ["R1", "R2", "c1", "c2"]
    assert dict(zip(all_["plasmid_id"], all_["representative"]))["c2"] == "R1"
    reps = db.candidates(d, sc, "representatives", tmp_path / "r.fasta")
    assert sorted(reps["plasmid_id"]) == ["R1", "R2"]
    assert (tmp_path / "a.fasta").read_text().count(">") == 4
    mash = tmp_path / "m.tsv"
    mash.write_text("R1\tQ__a\t0.01\t0\t1/1\nc1\tQ__a\t0.0\t0\t1/1\nc2\tQ__a\t0.02\t0\t1/1\nR2\tQ__a\t0.05\t0\t1/1\n")
    m = db.select_members(all_, mash, ["Q__a"], k=2)
    assert list(m.loc[m["selected"], "plasmid_id"]) == ["c1", "R1"]
    assert db.select_members(all_, mash, ["Q__a"], k=0)["selected"].all()
