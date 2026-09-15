#!/bin/bash
# Compactly transfer validated legacy ppi_screen scores into the V2 pair cache.

set -euo pipefail

: "${GBI_REPO_ROOT:?GBI_REPO_ROOT is required}"
: "${GBI_PPI_SCREEN_DIR:?GBI_PPI_SCREEN_DIR is required}"
: "${GBI_PAIR_CACHE:?GBI_PAIR_CACHE is required}"

cd "$GBI_REPO_ROOT"
if [[ -f .venv/bin/activate ]]; then
  source .venv/bin/activate
fi

python dynamic_search.py import-legacy-screen \
  --screen-dir "$GBI_PPI_SCREEN_DIR" \
  --cache "$GBI_PAIR_CACHE" \
  --audit-dir "$GBI_REPO_ROOT/runtime/legacy_ppi_screen_import"
