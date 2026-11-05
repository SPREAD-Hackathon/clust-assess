
from __future__ import annotations

import math
import random
import statistics
from itertools import combinations
from math import comb

INF = math.inf


def median_high(values) -> float:
    v = list(values)
    return statistics.median_high(v) if v else math.nan


def pair_lookup(rows) -> dict[tuple[str, str], tuple[float, float]]:
    out = {}
    for a, b, c, d in rows:
        out[(a, b) if a <= b else (b, a)] = (float(c), float(d))
    return out


def within(ids: list[str], lookup: dict, field: int) -> list[float]:
    vals = []
    for a, b in combinations(sorted(ids), 2):
        v = lookup.get((a, b))
        if v is not None:
            vals.append(v[field])
    return vals


def resample(members: list[str], n: int, lookup: dict, observed: float, n_draws: int = 1000,
             seed: int = 2026) -> dict:
    out = {"p_tighter": math.nan, "p_looser": math.nan, "n_subsets": 0, "exact": False}
    if n < 2 or len(members) < n or math.isnan(observed):
        return out
    members = sorted(members)
    total = comb(len(members), n)
    if total <= n_draws:
        subsets, exact = list(combinations(members, n)), True
    else:
        rng = random.Random(seed)
        subsets, exact = [tuple(rng.sample(members, n)) for _ in range(n_draws)], False
    meds = [median_high(within(list(s), lookup, 1)) for s in subsets]
    out.update(p_tighter=sum(m <= observed for m in meds) / len(meds),
               p_looser=sum(m >= observed for m in meds) / len(meds), n_subsets=len(meds), exact=exact)
    return out


def status(n_inputs: int, observed: float, baseline: float, p_tighter: float, p_looser: float,
           alpha: float = 0.05) -> str:
    if n_inputs < 2:
        return "single_input"
    if math.isnan(baseline):
        return "no_baseline"
    if not math.isnan(p_tighter):
        if p_tighter <= alpha and observed < baseline:
            return "tighter_than_db"
        if p_looser <= alpha and observed > baseline:
            return "looser_than_db"
        return "similar"
    return "tighter_than_db" if observed < baseline else "looser_than_db" if observed > baseline else "similar"


def fmt(x) -> str:
    if x is None:
        return "NA"
    if isinstance(x, float):
        if math.isnan(x) or math.isinf(x):
            return "NA"
        return str(int(x)) if x.is_integer() else f"{x:.4g}"
    return str(x)
