# Sample HTML report

- `final_output_extended.csv`: a sample pipeline output (15 OXA-48 plasmids, 105 pairs). Its subcommunity columns were relabelled from `community_0_subcommunity_0` to `community_0_subcommunity_7129` so that the pairs match the bundled reference database (the subcommunity ids of one database build are not the ids of another).
- `report.html`: the report built from it with the bundled reference database. Open it with a double-click.

Rebuild it from the repository root:

```
sh html_report/make_report.sh example_output/html_report/final_output_extended.csv example_output/html_report/report.html
```

See [`../../html_report/README.md`](../../html_report/README.md) for the input format and how the report works.
