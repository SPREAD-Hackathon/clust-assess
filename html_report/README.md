# Interactive HTML report

Turns the pipeline's `final_output.csv` and the reference pair table of the input's subcommunity into one self-contained `report.html`. It needs no Python, no internet and no server, and opens with a double-click (it can also be emailed). The reference table is embedded compressed, so a 17 MB table adds about 2 MB to the report. Opening it needs a current browser (Chrome, Edge, Firefox or Safari 16.4+).

An example report is in [`../example_output/report.html`](../example_output/report.html).

## Build it

Run from the repository root, after the pipeline has written `results/final_output.csv`:

```
sh html_report/make_report.sh                  # Mac / Linux
html_report\make_report.bat                    # Windows (double-clicking it also works)

# explicit paths, e.g. inside a workflow step: <results> <reference> <report>
sh html_report/make_report.sh results/final_output.csv results/reference_db.tsv results/report.html
```

Defaults are `results/final_output.csv`, `results/reference_db.tsv` and `results/report.html`. Exit code `0` = report written, `1` = failed (the reason is printed). With no arguments the report opens in the browser; with arguments it does not. Set `CA_NOOPEN=1` to never open it.

## Inputs

**`final_output.csv`** (the pipeline's main result; `;`, `,` or tab separated)

| Column | Required | Meaning |
|---|---|---|
| `Plasmid A`, `Plasmid B` | yes | the two plasmids compared |
| `DCJ` | yes | Pling DCJ-Indel distance, `NA` rows are skipped |
| `community` | no | subcommunity of the pair (`7129` or the full id). Otherwise taken from known reference plasmids in the file, or from the only subcommunity in the reference |
| `containment` | no | drawn as its own pairwise heatmap (not part of the judgement) |
| `mash` | no | Mash distance; drawn as its own pairwise heatmap when present (not part of the judgement) |
| `Plasmid A country` / `Plasmid B country`, `... host`, `... date` | no | metadata; each one unlocks its own section |
| `community_DCJ`, `community_containment` | no | ignored: the reference pairs are the baseline |

A pair with exactly one plasmid in the reference is treated as a measured input-to-reference pair: it is shown in the network and the plasmid tables and is not part of the judgement.

**`reference_db.tsv`** (tab separated): the Pling pair table of the subcommunity (or subcommunities) of the input, one row per pair of known plasmids.

| Column | Required | Meaning |
|---|---|---|
| `plasmid_1`, `plasmid_2` | yes | the two plasmids of the pair |
| `subcommunity` | yes (or `community`) | e.g. `community_0_subcommunity_7129` |
| `community` | no | the Pling community of the subcommunity |
| `dcj` | yes | DCJ distance, >= 0 |
| `containment` | no | read but not used (NA in the current table) |

## How the report judges the input

Within each subcommunity, the DCJ distribution of its known pairs is the yardstick. Each input pair is placed on it and given a "closer than X% of reference pairs" value; each subcommunity gets the average with a 95% interval (a plasmid-level bootstrap), where 50% means no difference from the reference. The classification mode (reference mean by default, a percentile, or fixed thresholds) is secondary: about half of the pairs are "closer than the mean" even when the input is no different from the reference.

A subcommunity that is missing from the reference, or whose reference has fewer than 10 pairs, is not assessed and is listed in the Summary. A reference whose 95th-percentile DCJ is above 25 is flagged as broad.

## What is in the report

Summary (checks, key figures), Distance distributions by subcommunity (with a Details view of which plasmid sits where), Pairwise heatmaps (DCJ, plus containment and Mash when present), Position of the input in the reference (rings and ranked pairs, with drill-down), Plasmid network (drill-down), and the optional Geography / Host / Dates, Plasmid statistics, Pair explorer and Method & keys sections. Every chart has an ⓘ explaining what it shows and how to read it.

## Files

```
index.html        page template (all charts and logic)
echarts.min.js    chart library (ECharts 5.5), embedded into the report
make_report.sh    builder for Mac / Linux (POSIX shell)
make_report.bat   builder for Windows (PowerShell, built into Windows)
```

## Notes and follow-ups

- The reference table is the one named in the Input section; the pipeline's own `baseline/*.member_pairs.tsv` files have the columns `a`, `b`, `containment`, `dcj`, `dcj_censored` and carry the subcommunity in the file name. Feeding them straight to the report needs a small step that combines them into one table with a `subcommunity` column (or an extension of the parser).
- Containment and Mash are shown as supplied; the reference currently holds DCJ only, so there is no baseline for them.
