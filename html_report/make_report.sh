#!/bin/sh
# Builds results/report.html from results/final_output_extended.csv and the reference database that ships with this folder (reference/). Nothing to install (POSIX shell, tr, cat).
# Run it, or call:  sh html_report/make_report.sh [final_output_extended.csv] [report.html] [reference_db.tsv(.gz)] [inputs.tsv]
HERE=$(cd "$(dirname "$0")" && pwd) || exit 1
ROOT=$(dirname "$HERE")
CSV=${1:-"$ROOT/results/final_output_extended.csv"}
OUT=${2:-"$ROOT/results/report.html"}
REF=${3:-"$HERE/reference/reference_db.tsv.gz"}   # bundled; pass another table (plain or .gz) to override
INP=${4:-"$(dirname "$CSV")/inputs.tsv"}            # optional: the pipeline's inputs.tsv (nearest database plasmid of each input); used when present

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
    '<script id="inputs-data" type="text/csv"></script>')
      if [ -s "$INP" ]; then printf '<script id="inputs-data" type="text/csv">
'; cat "$INP"; printf '
</script>
'; else printf '%s
' "$line"; fi ;;
    '<script id="ref-data" type="text/csv"></script>')
      printf '<script id="ref-data" type="text/plain" data-enc="gzip-base64">\n'; case "$REF" in *.gz) base64 < "$REF" ;; *) gzip -9 -c "$REF" | base64 ;; esac | tr -d '\n\r'; printf '\n</script>\n' ;;
    *) printf '%s\n' "$line" ;;
  esac
done > "$OUT.tmp" && mv "$OUT.tmp" "$OUT"

grep -q 'id="csv-data" type="text/csv">' "$OUT" && grep -q 'id="ref-data" type="text/plain"' "$OUT" || { echo "ERROR: index.html is missing its placeholder lines." >&2; exit 1; }
echo "Created $OUT"
# open the report only for interactive use (no arguments)
[ -n "$1" ] || [ -n "$CA_NOOPEN" ] && exit 0
if command -v open >/dev/null 2>&1; then open "$OUT"
elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$OUT" >/dev/null 2>&1; fi
