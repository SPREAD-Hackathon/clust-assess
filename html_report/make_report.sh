#!/bin/sh
# Builds results/report.html from results/final_output.csv and results/reference_db.tsv. Nothing to install (POSIX shell, tr, cat).
# Run it, or call:  sh html_report/make_report.sh [final_output.csv] [reference_db.tsv] [report.html]
HERE=$(cd "$(dirname "$0")" && pwd) || exit 1
ROOT=$(dirname "$HERE")
CSV=${1:-"$ROOT/results/final_output.csv"}
REF=${2:-"$ROOT/results/reference_db.tsv"}
OUT=${3:-"$ROOT/results/report.html"}

for f in "$CSV" "$REF" "$HERE/index.html" "$HERE/echarts.min.js"; do
  [ -f "$f" ] || { echo "ERROR: required file not found: $f" >&2; exit 1; }
done
[ -s "$CSV" ] || { echo "ERROR: $CSV is empty." >&2; exit 1; }
[ -s "$REF" ] || { echo "ERROR: $REF is empty." >&2; exit 1; }

# tr -d '\r' makes the page's line endings irrelevant (Windows or Unix)
tr -d '\r' < "$HERE/index.html" | while IFS= read -r line || [ -n "$line" ]; do
  case "$line" in
    '<script src="echarts.min.js"></script>')
      printf '<script>\n'; cat "$HERE/echarts.min.js"; printf '\n</script>\n' ;;
    '<script id="csv-data" type="text/csv"></script>')
      printf '<script id="csv-data" type="text/csv">\n'; cat "$CSV"; printf '\n</script>\n' ;;
    '<script id="ref-data" type="text/csv"></script>')
      printf '<script id="ref-data" type="text/plain" data-enc="gzip-base64">\n'; gzip -9 -c "$REF" | base64 | tr -d '\n\r'; printf '\n</script>\n' ;;
    *) printf '%s\n' "$line" ;;
  esac
done > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"

grep -q 'id="csv-data" type="text/csv">' "$OUT" && grep -q 'id="ref-data" type="text/plain"' "$OUT" || { echo "ERROR: index.html is missing its placeholder lines." >&2; exit 1; }
echo "Created $OUT"
# open the report only for interactive use (no arguments)
[ -n "$1" ] || [ -n "$CA_NOOPEN" ] && exit 0
if command -v open >/dev/null 2>&1; then open "$OUT"
elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$OUT" >/dev/null 2>&1; fi
