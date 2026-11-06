# Interactive HTML report

Turns the pipeline's `final_output_extended.csv` into one self-contained `report.html`, comparing your plasmids with a reference database that **ships with this folder**. It needs no Python, no internet and no server, and opens with a double-click (it can also be emailed). The reference is embedded compressed, so it adds only about 2 MB to the report. Opening the report needs a current browser (Chrome, Edge, Firefox or Safari 16.4+).

A sample input and the report built from it are in [`../example_output/html_report/`](../example_output/html_report/).

## Build it

Run from the repository root, after the pipeline has written `results/final_output_extended.csv`:

```
sh html_report/make_report.sh                  # Mac / Linux
html_report\make_report.bat                    # Windows (double-clicking it also works)

# explicit paths, e.g. inside a workflow step: <final_output_extended.csv> <report.html>
sh html_report/make_report.sh results/final_output_extended.csv results/report.html
```

Defaults are `results/final_output_extended.csv` and `results/report.html`. An optional fourth argument is the pipeline's `inputs.tsv`; by default `inputs.tsv` in the folder of the first file is used when it exists. Exit code `0` = report written, `1` = failed (the reason is printed). With no arguments the report opens in the browser; with arguments it does not. Set `CA_NOOPEN=1` to never open it.

## Input: `final_output_extended.csv`

The pipeline's main result (`;`, `,` or tab separated).

| Column | Required | Meaning |
|---|---|---|
| `Plasmid A`, `Plasmid B` | yes | the two plasmids compared |
| `DCJ` | yes | Pling DCJ-Indel distance, `NA` rows are skipped |
| `subcommunity_A`, `subcommunity_B` | recommended | subcommunity of each plasmid (`7129` or the full id, e.g. `community_0_subcommunity_7129`). A pair whose plasmids are in different subcommunities is not assessed |
| `dcj_censored` | no | `True` = Pling did not compute the DCJ (containment above its threshold); these pairs are left out of the judgement, counted, and drawn grey in the DCJ heatmap |
| `community` | no | older alternative to `subcommunity_A/B` |
| `containment` | no | Pling containment **distance** (0 = fully contained, 1 = nothing shared); drawn as its own pairwise heatmap (not part of the judgement) |
| `mash` | no | pairwise Mash distance (0 = identical); drawn as its own heatmap when the column is present. The pipeline does not write it yet: its outputs only carry the ACCIO Mash score of each input against its reference plasmid |
| `Plasmid A country` / `Plasmid B country`, `... host`, `... date` | no | metadata; each one unlocks its own section |
| `community_DCJ`, `community_containment`, `same_subcommunity`, `subcommunity_status`, `p_tighter` | no | ignored: the reference pairs are the baseline |

### Optional: `inputs.tsv`

The pipeline's table of inputs. Each input's nearest database plasmid (`nearest_db_plasmid`, `nearest_db_dcj`, `nearest_db_containment`) is added as an input-to-reference pair when that plasmid is in the reference, so the reference plasmid shows in the heatmaps and the network. Those known plasmids also place the input in their subcommunity, which takes priority over the subcommunity names in the results file. This is how inputs assigned to `community_0_subcommunity_0` by the pipeline's database are judged against subcommunity 7129 of the bundled reference (when the nearest plasmid is in both).

**Subcommunity ids belong to one database build.** If the file names subcommunities that are not in the bundled reference (for example `community_0_subcommunity_0` from another database), they are not assessed. When none can be assessed, the page says so and offers to judge them against the reference's only subcommunity anyway; use that only if the plasmids really belong to that group. A file with no subcommunity information at all is assigned to the reference's only subcommunity.

A pair with exactly one plasmid in the reference is treated as a measured input-to-reference pair: it is shown in the network and the plasmid tables and is not part of the judgement.

## The reference database (part of the program)

`reference/reference_db.tsv.gz` is the Pling pair table of the reference subcommunity: community 0, subcommunity 7129, 730 plasmids and 230,842 pairs, DCJ only. Users do not supply it.

- **Coverage:** a subcommunity that is not in the bundled table, or whose reference has fewer than 10 pairs, is not assessed and is listed in the Summary.
- **Format** (tab separated, plain or `.gz`): `plasmid_1`, `plasmid_2`, `subcommunity` (or `community`), `dcj` (>= 0); `containment` (distance) is optional and, when present, gives the background for the containment dots. The bundled table is DCJ only; the table the pipeline builds from `baseline/*.member_pairs.tsv` includes containment.
- **To update or extend it:** replace `reference/reference_db.tsv.gz` with a table in the same format (it can hold several subcommunities), or pass another table as a third argument: `sh html_report/make_report.sh results/final_output_extended.csv results/report.html my_reference.tsv`.

## How the report judges the input

Within each subcommunity, the DCJ distribution of its known pairs is the yardstick. Each input pair is placed on it and given a "closer than X% of reference pairs" value; each subcommunity gets the average with a 95% interval (a plasmid-level bootstrap), where 50% means no difference from the reference. The classification mode (reference mean by default, a percentile, or fixed thresholds) is secondary: about half of the pairs are "closer than the mean" even when the input is no different from the reference.

A reference whose 95th-percentile DCJ is above 25 is flagged as broad.

## What is in the report

Summary (checks, key figures), Distance distributions by subcommunity (with a Details view of which plasmid sits where), Pairwise heatmaps (DCJ, plus containment and Mash when present; darker = closer, scale 0 to twice the Pling threshold on a square-root scale, threshold marked on the colour bar, 1-3 dots for pairs closer than the background mean by 1, 1.5 and 2 SD (the analyst's z-score levels; background = the subcommunity's reference pairs when the reference has that measure, otherwise the input pairs themselves), grey cells for censored DCJ pairs, plasmids ordered by DCJ and containment together), Position of the input in the reference (rings and ranked pairs, with drill-down), Plasmid network (drill-down), and the optional Geography / Host / Dates, Plasmid statistics, Pair explorer and Method & keys sections. Every chart has an ⓘ explaining what it shows and how to read it.

## Files

```
index.html                        page template (all charts and logic)
echarts.min.js                    chart library (ECharts 5.5), embedded into the report
make_report.sh                    builder for Mac / Linux (POSIX shell)
make_report.bat                   builder for Windows (PowerShell, built into Windows)
reference/reference_db.tsv.gz     bundled reference database
```

## Follow-ups

- The pipeline builds the reference from its own `baseline/*.member_pairs.tsv` files in the `HTML_REPORT` step (see `README_workflow.md`); the bundled table is the fallback. The Pling thresholds used for the heatmap scale (DCJ 4, containment 0.5) are the pipeline defaults and are fixed in the page.
- Containment and Mash are shown as supplied; the reference currently holds DCJ only, so there is no baseline for them.
