#!/bin/bash
# Ordinary deterministic AlphaPulldown feature generation; no AI assistant is used.
# Resource requests are supplied by submit_round.sh so they remain easy to edit.

set -euo pipefail

module load alphapulldown

: "${GBI_SEQUENCE_FASTA:?GBI_SEQUENCE_FASTA is required}"
: "${GBI_FEATURE_DIR:?GBI_FEATURE_DIR is required}"
: "${GBI_MAX_TEMPLATE_DATE:=2026-09-11}"
: "${GBI_DB_PRESET:=full_dbs}"

case "$GBI_DB_PRESET" in
  full_dbs|reduced_dbs) ;;
  *)
    echo "GBI_DB_PRESET must be full_dbs or reduced_dbs (received: $GBI_DB_PRESET)" >&2
    exit 2
    ;;
esac

mkdir -p "$GBI_FEATURE_DIR"

arguments=(
  --fasta_paths="$GBI_SEQUENCE_FASTA"
  --output_dir="$GBI_FEATURE_DIR"
  --max_template_date="$GBI_MAX_TEMPLATE_DATE"
  --db_preset="$GBI_DB_PRESET"
  --save_msa_files=False
  --use_precomputed_msas=False
  --skip_existing=True
  --seq_index="$SLURM_ARRAY_TASK_ID"
)

if [[ -n "${GBI_ALPHAFOLD_DATA_DIR:-}" ]]; then
  arguments+=(--data_dir="$GBI_ALPHAFOLD_DATA_DIR")
fi

create_individual_features.py "${arguments[@]}"
