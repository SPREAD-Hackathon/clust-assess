<# : batch portion
@echo off
setlocal
rem  Builds results\report.html from results\final_output.csv and the reference database that ships with this folder (reference\). Nothing to install (uses PowerShell, built into Windows).
rem  Double-click it, or call:  make_report.bat [final_output.csv] [report.html] [reference_db.tsv(.gz)]
set "CA_SELF=%~f0"
set "CA_CSV=%~1"
set "CA_OUT=%~2"
set "CA_REF=%~3"
powershell -NoProfile -ExecutionPolicy Bypass -Command "& ([scriptblock]::Create([IO.File]::ReadAllText($env:CA_SELF)))"
if errorlevel 1 (
  if "%~1"=="" pause
  exit /b 1
)
exit /b 0
: end batch / begin PowerShell #>

$ErrorActionPreference = 'Stop'
try {
  $here = Split-Path -Parent $env:CA_SELF
  $root = Split-Path -Parent $here
  $csv  = if ($env:CA_CSV) { [IO.Path]::GetFullPath($env:CA_CSV) } else { Join-Path (Join-Path $root 'results') 'final_output.csv' }
  $ref  = if ($env:CA_REF) { [IO.Path]::GetFullPath($env:CA_REF) } else { Join-Path (Join-Path $here 'reference') 'reference_db.tsv.gz' }
  $out  = if ($env:CA_OUT) { [IO.Path]::GetFullPath($env:CA_OUT) } else { Join-Path (Join-Path $root 'results') 'report.html' }
  $page = Join-Path $here 'index.html'
  $lib  = Join-Path $here 'echarts.min.js'
  foreach ($f in @($csv, $ref, $page, $lib)) {
    if (-not (Test-Path -LiteralPath $f)) { throw "Required file not found: $f" }
  }
  $utf8 = New-Object System.Text.UTF8Encoding($false)
  $csvText = [IO.File]::ReadAllText($csv, $utf8).TrimStart([char]0xFEFF)
  if (-not $csvText.Trim()) { throw "$csv is empty." }
  $refBytes = [IO.File]::ReadAllBytes($ref)
  if ($refBytes.Length -eq 0) { throw "$ref is empty." }
  if ($ref -like '*.gz') { $refB64 = [Convert]::ToBase64String($refBytes) } else {
    $ms = New-Object IO.MemoryStream
    $gz = New-Object IO.Compression.GZipStream($ms, [IO.Compression.CompressionLevel]::Optimal)
    $gz.Write($refBytes, 0, $refBytes.Length); $gz.Close()
    $refB64 = [Convert]::ToBase64String($ms.ToArray())
  }
  $html = [IO.File]::ReadAllText($page, $utf8)
  $js   = [IO.File]::ReadAllText($lib,  $utf8).Replace('</script', '<\/script')
  $libTag = '<script src="echarts.min.js"></script>'
  $csvTag = '<script id="csv-data" type="text/csv"></script>'
  $refTag = '<script id="ref-data" type="text/csv"></script>'
  if (-not $html.Contains($libTag) -or -not $html.Contains($csvTag) -or -not $html.Contains($refTag)) { throw 'index.html is missing its placeholder lines.' }
  $html = $html.Replace($libTag, "<script>`n$js`n</script>")
  $html = $html.Replace($csvTag, "<script id=`"csv-data`" type=`"text/csv`">`n$csvText`n</script>")
  $html = $html.Replace($refTag, "<script id=`"ref-data`" type=`"text/plain`" data-enc=`"gzip-base64`">`n$refB64`n</script>")
  [IO.File]::WriteAllText($out, $html, $utf8)
  Write-Host "Created $out"
  if (-not $env:CA_CSV -and -not $env:CA_NOOPEN) { Start-Process -FilePath $out }
} catch {
  Write-Host ''
  Write-Host "ERROR: $($_.Exception.Message)" -ForegroundColor Red
  exit 1
}
