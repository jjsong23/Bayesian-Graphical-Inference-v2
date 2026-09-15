#!/bin/bash
# Submit one prepared backward-search round to NIH Biowulf.
# Usage: bash biowulf/submit_round.sh runs/my_search/round_000

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: bash biowulf/submit_round.sh <round-directory>" >&2
  exit 2
fi

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
round_dir=$(cd "$1" && pwd)
# Features are deliberately run-local. That makes post-frontier cleanup safe:
# one run can never delete another run's reusable feature objects.
search_dir=$(dirname "$round_dir")
feature_dir="${GBI_FEATURE_DIR:-$search_dir/features}"
max_template_date="${GBI_MAX_TEMPLATE_DATE:-2026-09-11}"
max_parallel="${GBI_MAX_CONCURRENT_GPU_JOBS:-10}"

mkdir -p "$round_dir/logs" "$feature_dir"
sequence_count=$(grep -c '^>' "$round_dir/sequences.fasta")
pair_count=$(( $(wc -l < "$round_dir/pairs.tsv") - 1 ))
if [[ "$sequence_count" -lt 1 || "$pair_count" -lt 1 ]]; then
  echo "The round package has no sequences or pairs." >&2
  exit 2
fi

num_cycle=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["num_cycle"])' "$round_dir/round_config.json")
num_predictions=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["num_predictions_per_model"])' "$round_dir/round_config.json")
model_names=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["model_names"])' "$round_dir/round_config.json")

feature_job=$(sbatch --parsable \
  --job-name=gbi_features \
  --cpus-per-task=8 --mem=50g --time=10:00:00 \
  --array="1-${sequence_count}%10" \
  --output="$round_dir/logs/features_%A_%a.out" \
  --error="$round_dir/logs/features_%A_%a.err" \
  --export="ALL,GBI_SEQUENCE_FASTA=$round_dir/sequences.fasta,GBI_FEATURE_DIR=$feature_dir,GBI_MAX_TEMPLATE_DATE=$max_template_date,GBI_ALPHAFOLD_DATA_DIR=${GBI_ALPHAFOLD_DATA_DIR:-}" \
  "$repo_root/biowulf/run_features.sbatch")

prediction_job=$(sbatch --parsable \
  --dependency="afterok:$feature_job" \
  --job-name=gbi_pairs \
  --cpus-per-task=8 --mem=50g --time=3:00:00 \
  --gres=lscratch:50,gpu:a100:1 --partition=gpu \
  --array="1-${pair_count}%${max_parallel}" \
  --output="$round_dir/logs/pairs_%A_%a.out" \
  --error="$round_dir/logs/pairs_%A_%a.err" \
  --export="ALL,GBI_ROUND_DIR=$round_dir,GBI_FEATURE_DIR=$feature_dir,GBI_NUM_CYCLE=$num_cycle,GBI_NUM_PREDICTIONS=$num_predictions,GBI_MODEL_NAMES=$model_names" \
  "$repo_root/biowulf/run_pairs.sbatch")

printf 'stage\tjob_id\nfeatures\t%s\npredictions\t%s\n' \
  "$feature_job" "$prediction_job" > "$round_dir/submission.tsv"

echo "Feature job: $feature_job"
echo "Prediction job: $prediction_job"
echo "After prediction completes:"
echo "  python dynamic_search.py collect --run-dir $(dirname "$round_dir")"
echo "  python dynamic_search.py step --run-dir $(dirname "$round_dir")"
