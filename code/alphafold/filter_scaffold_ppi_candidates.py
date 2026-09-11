"""Filter AlphaFold candidates to pairs involving adaptor/scaffold proteins.

The project taxonomy has no literal ``structural_protein`` class. This analysis
therefore uses ``adaptor_scaffold`` as an explicit, reproducible proxy and
reports both:

1. at least one scaffold endpoint (includes scaffold--scaffold), and
2. exactly one scaffold endpoint (scaffold--non-scaffold only).

Only adjacency-aware colocalization-supported pairs are considered. Pairs with
Tier A or Tier B experimental PPI support are removed; ambiguous-only prior
records are retained.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


PROJECT_ROOT = Path(
    r"C:\Users\songjj\Documents\Codex\2026-07-21\un\graphical_bayesian_inference"
)
UNIVERSE = (
    PROJECT_ROOT
    / "data"
    / "node_selection"
    / "alternate_no_kinase_activity"
    / "processed"
    / "node_universe_no_kinase_activity.tsv"
)
SUPPORTED_PAIRS = (
    PROJECT_ROOT
    / "results"
    / "colocalization_adjacency_aware_no_kinase_activity_subset"
    / "supported_colocalization_edges.tsv.gz"
)
EXPERIMENTAL_PAIRS = (
    PROJECT_ROOT
    / "results"
    / "experimental_ppi_tiers_621"
    / "experimental_ppi_pairs.tsv"
)
RUNTIME_ESTIMATES = (
    PROJECT_ROOT
    / "results"
    / "alphafold_runtime_estimate_621_proteins"
    / "supported_pair_runtime_estimates.tsv.gz"
)
OUTPUT_DIR = PROJECT_ROOT / "results" / "alphafold_scaffold_pair_filter"
PROTEIN_INDEX = (
    PROJECT_ROOT
    / "data"
    / "edge_characterization"
    / "kinase_predictor"
    / "phosphosite_database"
    / "protein_index.tsv"
)

STRUCTURAL_PROXY_CLASS = "adaptor_scaffold"
STARTUP_SECONDS = 73.0
DEFAULT_PAIRS_PER_JOB = 10
FOLDING_COEFFICIENT = 1.29
FOLDING_EXPONENT = 0.95
TOKENS_PER_RESIDUE = 1.2


def canonical_pair(node_a: str, node_b: str) -> tuple[str, str]:
    return tuple(sorted((node_a, node_b), key=str.casefold))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe-tsv", type=Path, default=UNIVERSE)
    parser.add_argument("--supported-pairs", type=Path, default=SUPPORTED_PAIRS)
    parser.add_argument("--experimental-pairs", type=Path, default=EXPERIMENTAL_PAIRS)
    parser.add_argument("--protein-index", type=Path, default=PROTEIN_INDEX)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.universe_tsv.open("r", encoding="utf-8", newline="") as handle:
        universe_rows = list(csv.DictReader(handle, delimiter="\t"))
    symbol_column = (
        "gene_symbol" if "gene_symbol" in universe_rows[0] else "symbol"
    )
    protein_rows = [
        row
        for row in universe_rows
        if not row.get("node_type") or row["node_type"] == "protein"
    ]
    all_symbols = {row[symbol_column] for row in protein_rows}
    scaffold_symbols = {
        row[symbol_column]
        for row in protein_rows
        if STRUCTURAL_PROXY_CLASS in row["classes"].split(";")
    }
    non_scaffold_symbols = all_symbols - scaffold_symbols

    known_tier_ab: set[tuple[str, str]] = set()
    ambiguous_only: set[tuple[str, str]] = set()
    with args.experimental_pairs.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pair = canonical_pair(row["node_a"], row["node_b"])
            if row["pair_tier"] in {"A", "B"}:
                known_tier_ab.add(pair)
            elif row["pair_tier"] == "ambiguous_only":
                ambiguous_only.add(pair)

    sequence_lengths: dict[str, int] = {}
    with args.protein_index.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["symbol"] in all_symbols and row["sequence_length"]:
                sequence_lengths[row["symbol"]] = int(row["sequence_length"])
    missing_sequence_lengths = sorted(
        all_symbols - set(sequence_lengths), key=str.casefold
    )
    if missing_sequence_lengths:
        raise ValueError(
            "Missing sequence lengths for: " + ", ".join(missing_sequence_lengths)
        )

    counts = {
        "colocalization_supported_total": 0,
        "colocalization_supported_after_tier_ab_exclusion": 0,
        "supported_at_least_one_scaffold_before_prior_exclusion": 0,
        "supported_exactly_one_scaffold_before_prior_exclusion": 0,
        "supported_scaffold_scaffold_before_prior_exclusion": 0,
        "at_least_one_scaffold_after_tier_ab_exclusion": 0,
        "exactly_one_scaffold_after_tier_ab_exclusion": 0,
        "scaffold_scaffold_after_tier_ab_exclusion": 0,
        "non_scaffold_non_scaffold_after_tier_ab_exclusion": 0,
        "tier_ab_excluded_from_at_least_one_scaffold": 0,
        "tier_ab_excluded_from_exactly_one_scaffold": 0,
        "tier_ab_excluded_from_scaffold_scaffold": 0,
        "ambiguous_only_retained_in_at_least_one_scaffold": 0,
    }
    at_least_one_rows: list[dict[str, object]] = []
    exactly_one_rows: list[dict[str, object]] = []

    with gzip.open(
        args.supported_pairs, "rt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_headers = list(reader.fieldnames or [])
        for row in reader:
            counts["colocalization_supported_total"] += 1
            node_a = row["node_a"]
            node_b = row["node_b"]
            pair = canonical_pair(node_a, node_b)
            a_scaffold = node_a in scaffold_symbols
            b_scaffold = node_b in scaffold_symbols
            at_least_one = a_scaffold or b_scaffold
            exactly_one = a_scaffold != b_scaffold
            both_scaffold = a_scaffold and b_scaffold

            if at_least_one:
                counts[
                    "supported_at_least_one_scaffold_before_prior_exclusion"
                ] += 1
            if exactly_one:
                counts[
                    "supported_exactly_one_scaffold_before_prior_exclusion"
                ] += 1
            if both_scaffold:
                counts[
                    "supported_scaffold_scaffold_before_prior_exclusion"
                ] += 1

            if pair in known_tier_ab:
                if at_least_one:
                    counts["tier_ab_excluded_from_at_least_one_scaffold"] += 1
                if exactly_one:
                    counts["tier_ab_excluded_from_exactly_one_scaffold"] += 1
                if both_scaffold:
                    counts["tier_ab_excluded_from_scaffold_scaffold"] += 1
                continue

            counts["colocalization_supported_after_tier_ab_exclusion"] += 1
            if at_least_one:
                counts["at_least_one_scaffold_after_tier_ab_exclusion"] += 1
                combined_length = sequence_lengths[node_a] + sequence_lengths[node_b]
                estimated_tokens = TOKENS_PER_RESIDUE * combined_length
                folding_seconds = FOLDING_COEFFICIENT * (
                    estimated_tokens**FOLDING_EXPONENT
                )
                output_row = {
                    **row,
                    "structural_proxy_class": STRUCTURAL_PROXY_CLASS,
                    "node_a_is_scaffold": a_scaffold,
                    "node_b_is_scaffold": b_scaffold,
                    "pair_category": (
                        "scaffold_scaffold" if both_scaffold else "scaffold_other"
                    ),
                    "ambiguous_prior_record_retained": pair in ambiguous_only,
                    "combined_length_aa": combined_length,
                    "estimated_tokens": estimated_tokens,
                    "estimated_folding_seconds_benchmark_configuration": (
                        folding_seconds
                    ),
                }
                at_least_one_rows.append(output_row)
                if exactly_one:
                    counts["exactly_one_scaffold_after_tier_ab_exclusion"] += 1
                    exactly_one_rows.append(output_row)
                else:
                    counts["scaffold_scaffold_after_tier_ab_exclusion"] += 1
                if pair in ambiguous_only:
                    counts[
                        "ambiguous_only_retained_in_at_least_one_scaffold"
                    ] += 1
            else:
                counts["non_scaffold_non_scaffold_after_tier_ab_exclusion"] += 1

    output_headers = [
        *source_headers,
        "structural_proxy_class",
        "node_a_is_scaffold",
        "node_b_is_scaffold",
        "pair_category",
        "ambiguous_prior_record_retained",
        "combined_length_aa",
        "estimated_tokens",
        "estimated_folding_seconds_benchmark_configuration",
    ]
    for path, rows in [
        (
            args.output_dir / "scaffold_involving_candidates.tsv.gz",
            at_least_one_rows,
        ),
        (
            args.output_dir / "scaffold_to_non_scaffold_candidates.tsv.gz",
            exactly_one_rows,
        ),
    ]:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=output_headers,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    def runtime_summary(rows: list[dict[str, object]]) -> dict[str, float]:
        folding_seconds = sum(
            float(row["estimated_folding_seconds_benchmark_configuration"])
            for row in rows
        )
        job_count = (len(rows) + DEFAULT_PAIRS_PER_JOB - 1) // DEFAULT_PAIRS_PER_JOB
        startup_seconds = STARTUP_SECONDS * job_count
        total_seconds = folding_seconds + startup_seconds
        one_prediction_approx_seconds = folding_seconds / 5.0 + startup_seconds
        return {
            "pair_count": len(rows),
            "folding_gpu_hours_benchmark_configuration": folding_seconds / 3600,
            "serial_days_benchmark_configuration_pairs_per_job_10": (
                total_seconds / 86400
            ),
            "serial_days_one_prediction_per_model_approximation": (
                one_prediction_approx_seconds / 86400
            ),
        }

    eligible_after_prior = counts[
        "colocalization_supported_after_tier_ab_exclusion"
    ]
    original_pair_count = len(all_symbols) * (len(all_symbols) - 1) // 2
    summary = {
        "structural_definition": {
            "requested_term": "structural proteins",
            "literal_structural_class_present": False,
            "proxy_used": STRUCTURAL_PROXY_CLASS,
            "proxy_node_count": len(scaffold_symbols),
            "non_proxy_node_count": len(non_scaffold_symbols),
            "interpretation_caveat": (
                "adaptor_scaffold is a signaling-role class and is not identical "
                "to a cytoskeletal/extracellular-matrix structural-protein set."
            ),
        },
        "counts": counts,
        "reductions": {
            "at_least_one_scaffold_vs_eligible_after_prior": (
                1.0
                - counts["at_least_one_scaffold_after_tier_ab_exclusion"]
                / eligible_after_prior
            ),
            "exactly_one_scaffold_vs_eligible_after_prior": (
                1.0
                - counts["exactly_one_scaffold_after_tier_ab_exclusion"]
                / eligible_after_prior
            ),
            "at_least_one_scaffold_vs_original_pair_space": (
                1.0
                - counts["at_least_one_scaffold_after_tier_ab_exclusion"]
                / original_pair_count
            ),
            "exactly_one_scaffold_vs_original_pair_space": (
                1.0
                - counts["exactly_one_scaffold_after_tier_ab_exclusion"]
                / original_pair_count
            ),
        },
        "runtime": {
            "at_least_one_scaffold": runtime_summary(at_least_one_rows),
            "exactly_one_scaffold": runtime_summary(exactly_one_rows),
        },
        "inputs": {
            "protein_count": len(all_symbols),
            "original_undirected_pair_count": original_pair_count,
            "universe_tsv": str(args.universe_tsv),
            "supported_pairs": str(args.supported_pairs),
            "experimental_pairs": str(args.experimental_pairs),
            "protein_index": str(args.protein_index),
        },
    }
    (args.output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "scaffold_symbols.txt").write_text(
        "\n".join(sorted(scaffold_symbols, key=str.casefold)) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
