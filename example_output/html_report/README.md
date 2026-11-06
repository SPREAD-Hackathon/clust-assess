# Sample HTML report

- `final_output.csv`: a sample pipeline output (15 OXA-48 plasmids plus one known reference plasmid, `NZ_CP117991.1`).
- `report.html`: the report built from it with the bundled reference database. Open it with a double-click.

Rebuild it from the repository root:

```
sh html_report/make_report.sh example_output/html_report/final_output.csv example_output/html_report/report.html
```

See [`../../html_report/README.md`](../../html_report/README.md) for the input format and how the report works.
