#!/usr/bin/env bash
# Run the API first: Rscript api/run.R
# Then: bash api/test_api.sh
set -euo pipefail

BASE="${BASE:-http://localhost:8000}"

echo "== /health =="
curl -s "${BASE}/health" | python3 -m json.tool
echo

echo "== /forecast?horizon=3 =="
curl -s "${BASE}/forecast?horizon=3" | python3 -m json.tool
echo

echo "== /model-info (keys only) =="
curl -s "${BASE}/model-info" | python3 -c 'import sys,json; d=json.load(sys.stdin); print(list(d.keys()))'
echo

echo "== assert avg30_garch matches latest.rds =="
Rscript -e '
  latest <- readRDS("outputs/latest.rds")
  api    <- jsonlite::fromJSON(Sys.getenv("BASE","http://localhost:8000") %+% "/forecast")
  stopifnot(abs(api$avg30_garch - latest$forecast$avg30_garch) < 1e-6)
  cat("OK avg30_garch =", api$avg30_garch, "\n")
' 2>/dev/null || Rscript -e '
  latest <- readRDS("outputs/latest.rds")
  api    <- jsonlite::fromJSON(paste0(Sys.getenv("BASE","http://localhost:8000"), "/forecast"))
  stopifnot(abs(api$avg30_garch - latest$forecast$avg30_garch) < 1e-6)
  cat("OK avg30_garch =", api$avg30_garch, "\n")
'