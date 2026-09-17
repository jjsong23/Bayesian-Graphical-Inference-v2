#!/bin/bash
# One manifest pair per GPU-array task. This is intentionally non-interactive.

set -euo pipefail

module load alphapulldown

: "${GBI_ROUND_DIR:?GBI_ROUND_DIR is required}"
: "${GBI_FEATURE_DIR:?GBI_FEATURE_DIR is required}"
: "${GBI_NUM_CYCLE:=3}"
: "${GBI_NUM_PREDICTIONS:=1}"
: "${GBI_MODEL_NAMES:=model_1_multimer_v3}"
: "${GBI_DB_PRESET:=full_dbs}"

manifest="$GBI_ROUND_DIR/pairs.tsv"
record=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$manifest")
if [[ -z "$record" ]]; then
  echo "No pair for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

IFS=$'\t' read -r pair_id node_a node_b backward_parent backward_current <<< "$record"
work_root="${LSCRATCH:-${TMPDIR:-/tmp}}/gbi_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$work_root"
trap 'rm -rf -- "$work_root"' EXIT
printf '%s\n' "$node_a" > "$work_root/bait.txt"
printf '%s\n' "$node_b" > "$work_root/candidate.txt"

scratch_output="$work_root/models"
compact_output="$GBI_ROUND_DIR/models/$pair_id"
mkdir -p "$scratch_output" "$compact_output"

run_multimer_jobs.py \
  --mode=pulldown \
  --num_cycle="$GBI_NUM_CYCLE" \
  --num_predictions_per_model="$GBI_NUM_PREDICTIONS" \
  --model_names="$GBI_MODEL_NAMES" \
  --compress_result_pickles=True \
  --remove_result_pickles=True \
  --output_path="$scratch_output" \
  --protein_lists="$work_root/bait.txt,$work_root/candidate.txt" \
  --monomer_objects_dir="$GBI_FEATURE_DIR" \
  --job_index=1

# Analyze on node-local scratch. A task is successful only if the analysis
# creates a nonempty table containing at least one finite ipTM value in [0,1].
(
  cd "$scratch_output"
  run_get_good_pae.sh --output_dir .
)

score_file=$(find "$scratch_output" -type f -name predictions_with_good_interpae.csv -print -quit)
if [[ -z "$score_file" ]]; then
  echo "AlphaPulldown did not produce predictions_with_good_interpae.csv for $pair_id" >&2
  exit 3
fi
python3 -c 'import math,pandas as pd,sys; f=pd.read_csv(sys.argv[1]); c={str(x).strip().lower():x for x in f.columns}.get("iptm"); v=[] if c is None else pd.to_numeric(f[c],errors="coerce").dropna().tolist(); assert any(math.isfinite(float(x)) and 0.0 <= float(x) <= 1.0 for x in v), "no valid ipTM score"' "$score_file"
top_structure=$(find "$scratch_output" -type f -name ranked_0.pdb -size +0c -print -quit)
if [[ -z "$top_structure" ]]; then
  echo "AlphaPulldown did not produce a nonempty ranked_0.pdb for $pair_id" >&2
  exit 4
fi

# Return only the compact scientific record to /data. Large result pickles,
# PAE JSON files, and all other temporary model products disappear with
# $LSCRATCH at task exit.
cp -f "$score_file" "$compact_output/predictions_with_good_interpae.csv"
cp -f "$top_structure" "$compact_output/ranked_0.pdb"
for artifact_name in ranking_debug.json timings.json; do
  artifact=$(find "$scratch_output" -type f -name "$artifact_name" -print -quit)
  if [[ -n "$artifact" ]]; then
    cp -f "$artifact" "$compact_output/$artifact_name"
  fi
done
printf 'pair_id\tnode_a\tnode_b\tmodel_names\tnum_cycle\tnum_predictions_per_model\tfeature_db_preset\tscratch_job_id\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$pair_id" "$node_a" "$node_b" "$GBI_MODEL_NAMES" "$GBI_NUM_CYCLE" \
  "$GBI_NUM_PREDICTIONS" "$GBI_DB_PRESET" "$SLURM_JOB_ID" \
  > "$compact_output/prediction_manifest.tsv"
