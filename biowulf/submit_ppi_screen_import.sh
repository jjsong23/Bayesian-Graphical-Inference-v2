#!/bin/bash
# Submit the compact legacy-screen transfer. No model/MSA files are copied.

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
screen_dir="${1:-/data/$USER/ppi_screen}"
cache_path="${2:-$repo_root/runtime/pair_cache.sqlite3}"

if [[ ! -d "$screen_dir/scores" ]]; then
  echo "Legacy scores directory is missing: $screen_dir/scores" >&2
  exit 2
fi
mkdir -p "$(dirname "$cache_path")" "$repo_root/runtime/legacy_ppi_screen_import"

job_id=$(sbatch --parsable \
  --job-name=gbi_import_ppi \
  --cpus-per-task=2 --mem=4g --time=01:00:00 \
  --output="$repo_root/runtime/legacy_ppi_screen_import/import_%j.out" \
  --error="$repo_root/runtime/legacy_ppi_screen_import/import_%j.err" \
  --export="ALL,GBI_REPO_ROOT=$repo_root,GBI_PPI_SCREEN_DIR=$screen_dir,GBI_PAIR_CACHE=$cache_path" \
  "$repo_root/biowulf/import_ppi_screen.sbatch")

echo "Submitted compact ppi_screen validation/import as Slurm job $job_id"
echo "Monitor with: squeue -j $job_id"
