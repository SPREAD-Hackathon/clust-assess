# Sample HTML report

- `final_output_extended.csv` and `inputs.tsv`: a sample pipeline output (15 OXA-48 plasmids, 105 pairs). The pipeline placed them in `community_0_subcommunity_0` of its own database; the report places them in subcommunity 7129 of the bundled reference because their nearest database plasmid, `NZ_CP117991.1` (listed in `inputs.tsv`), belongs to it.
- `report.html`: the report built from them with the bundled reference database. Open it with a double-click.

Rebuild it from the repository root (`inputs.tsv` is picked up from the folder of the first file):

```
sh html_report/make_report.sh example_output/html_report/final_output_extended.csv example_output/html_report/report.html
```

See [`../../html_report/README.md`](../../html_report/README.md) for the input format and how the report works.
