"""Estimate AlphaFold PPI runtime for the alternate colocalization graph.

Runtime model supplied by the user:
    startup_seconds = 73 per submitted job
    folding_seconds = 1.29 * tokens ** 0.95

The two supplied calibration examples imply approximately 1.2 model tokens per
amino-acid residue in the combined pair:
    Calm1 + Nme2: 301 aa -> 361.2 tokens -> 347.1 s
    Camk2d + Pde4d: 1246 aa -> 1495.2 tokens -> 1338.3 s

The script evaluates three useful scopes:
  1. supported colocalization edges (> neutral baseline),
  2. conservative AlphaFold-retained pairs (not explicitly disjoint), and
  3. every possible pair in the 621-protein universe.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import statistics
from pathlib import Path


PROJECT_ROOT = Path(
    r"C:\Users\songjj\Documents\Codex\2026-07-21\un\graphical_bayesian_inference"
)
PROTEIN_INDEX = (
    PROJECT_ROOT
    / "data"
    / "edge_characterization"
    / "kinase_predictor"
    / "phosphosite_database"
    / "protein_index.tsv"
)
ALTERNATE_UNIVERSE = (
    PROJECT_ROOT
    / "data"
    / "node_selection"
    / "alternate_no_kinase_activity"
    / "processed"
    / "node_universe_no_kinase_activity.tsv"
)
PAIR_DIR = (
    PROJECT_ROOT
    / "results"
    / "colocalization_adjacency_aware_no_kinase_activity_subset"
)
OUTPUT_DIR = PROJECT_ROOT / "results" / "alphafold_runtime_estimate_621_proteins"

STARTUP_SECONDS = 73.0
FOLDING_COEFFICIENT = 1.29
FOLDING_EXPONENT = 0.95
TOKENS_PER_RESIDUE = 1.2
BATCH_SIZES = [1, 5, 10, 25, 50, 100, 250, 500, 1000]
PARALLEL_WORKERS = [1, 4, 8, 16, 32, 64, 128]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def percentile(sorted_values: list[float], probability: float) -> float:
    if not sorted_values:
        return math.nan
    position = (len(sorted_values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return (
        sorted_values[lower] * (1.0 - fraction)
        + sorted_values[upper] * fraction
    )


def format_duration(seconds: float) -> dict[str, float]:
    return {
        "seconds": seconds,
        "hours": seconds / 3600,
        "days": seconds / 86400,
        "years_365_25_days": seconds / (365.25 * 86400),
    }


def load_lengths() -> tuple[dict[str, int], dict[str, object]]:
    selected_rows = read_tsv(ALTERNATE_UNIVERSE)
    selected_symbols = {row["gene_symbol"] for row in selected_rows}
    index_rows = read_tsv(PROTEIN_INDEX)
    length_by_symbol = {
        row["symbol"]: int(row["sequence_length"])
        for row in index_rows
        if row["symbol"] in selected_symbols and row["sequence_length"]
    }
    missing = sorted(selected_symbols - set(length_by_symbol))
    if missing:
        raise ValueError(f"Missing sequence lengths for {len(missing)} proteins: {missing}")
    lengths = sorted(length_by_symbol.values())
    audit = {
        "selected_proteins": len(selected_symbols),
        "proteins_with_sequence_length": len(length_by_symbol),
        "missing_sequence_lengths": missing,
        "sequence_length_summary_aa": {
            "minimum": min(lengths),
            "median": statistics.median(lengths),
            "mean": statistics.fmean(lengths),
            "p90": percentile([float(value) for value in lengths], 0.90),
            "p95": percentile([float(value) for value in lengths], 0.95),
            "p99": percentile([float(value) for value in lengths], 0.99),
            "maximum": max(lengths),
        },
    }
    return length_by_symbol, audit


def evaluate_pair_file(
    path: Path,
    length_by_symbol: dict[str, int],
    per_pair_writer: csv.DictWriter | None = None,
) -> dict[str, object]:
    folding_times: list[float] = []
    combined_lengths: list[float] = []
    token_counts: list[float] = []
    maximum_record: dict[str, object] | None = None
    minimum_record: dict[str, object] | None = None

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            node_a = row["node_a"]
            node_b = row["node_b"]
            length_a = length_by_symbol[node_a]
            length_b = length_by_symbol[node_b]
            combined_length = length_a + length_b
            tokens = TOKENS_PER_RESIDUE * combined_length
            folding_seconds = FOLDING_COEFFICIENT * tokens**FOLDING_EXPONENT
            record = {
                "edge_id": row["edge_id"],
                "node_a": node_a,
                "node_b": node_b,
                "length_a_aa": length_a,
                "length_b_aa": length_b,
                "combined_length_aa": combined_length,
                "estimated_tokens": tokens,
                "estimated_folding_seconds": folding_seconds,
                "estimated_pair_seconds_if_one_pair_per_job": (
                    folding_seconds + STARTUP_SECONDS
                ),
            }
            if per_pair_writer is not None:
                per_pair_writer.writerow(record)
            folding_times.append(folding_seconds)
            combined_lengths.append(float(combined_length))
            token_counts.append(tokens)
            if (
                maximum_record is None
                or folding_seconds
                > float(maximum_record["estimated_folding_seconds"])
            ):
                maximum_record = record
            if (
                minimum_record is None
                or folding_seconds
                < float(minimum_record["estimated_folding_seconds"])
            ):
                minimum_record = record

    sorted_folding = sorted(folding_times)
    sorted_lengths = sorted(combined_lengths)
    sorted_tokens = sorted(token_counts)
    total_folding_seconds = sum(folding_times)
    return {
        "pair_count": len(folding_times),
        "combined_length_aa": {
            "minimum": min(combined_lengths),
            "median": statistics.median(combined_lengths),
            "mean": statistics.fmean(combined_lengths),
            "p90": percentile(sorted_lengths, 0.90),
            "p95": percentile(sorted_lengths, 0.95),
            "p99": percentile(sorted_lengths, 0.99),
            "maximum": max(combined_lengths),
        },
        "estimated_tokens": {
            "minimum": min(token_counts),
            "median": statistics.median(token_counts),
            "mean": statistics.fmean(token_counts),
            "p90": percentile(sorted_tokens, 0.90),
            "p95": percentile(sorted_tokens, 0.95),
            "p99": percentile(sorted_tokens, 0.99),
            "maximum": max(token_counts),
        },
        "folding_seconds_per_pair": {
            "minimum": min(folding_times),
            "median": statistics.median(folding_times),
            "mean": statistics.fmean(folding_times),
            "p90": percentile(sorted_folding, 0.90),
            "p95": percentile(sorted_folding, 0.95),
            "p99": percentile(sorted_folding, 0.99),
            "maximum": max(folding_times),
        },
        "total_folding_time": format_duration(total_folding_seconds),
        "minimum_runtime_pair": minimum_record,
        "maximum_runtime_pair": maximum_record,
    }


def add_batching_scenarios(scope: dict[str, object]) -> None:
    pair_count = int(scope["pair_count"])
    folding_seconds = float(scope["total_folding_time"]["seconds"])
    scenarios = []
    for batch_size in BATCH_SIZES:
        jobs = math.ceil(pair_count / batch_size)
        startup_seconds = STARTUP_SECONDS * jobs
        total_seconds = folding_seconds + startup_seconds
        scenarios.append(
            {
                "pairs_per_job": batch_size,
                "job_count": jobs,
                "startup_time": format_duration(startup_seconds),
                "total_serial_time": format_duration(total_seconds),
                "startup_fraction": startup_seconds / total_seconds,
            }
        )
    scope["batching_scenarios"] = scenarios

    reference_batch_size = 10
    reference = next(
        scenario
        for scenario in scenarios
        if scenario["pairs_per_job"] == reference_batch_size
    )
    reference_seconds = float(reference["total_serial_time"]["seconds"])
    scope["ideal_parallel_wall_time_pairs_per_job_10"] = [
        {
            "parallel_workers": workers,
            "ideal_wall_time": format_duration(reference_seconds / workers),
        }
        for workers in PARALLEL_WORKERS
    ]


def write_scenario_table(
    path: Path, scopes: dict[str, dict[str, object]]
) -> None:
    headers = [
        "scope",
        "pair_count",
        "pairs_per_job",
        "job_count",
        "startup_days",
        "folding_days",
        "total_serial_days",
        "total_serial_years",
        "startup_fraction",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=headers, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for scope_name, scope in scopes.items():
            for scenario in scope["batching_scenarios"]:
                writer.writerow(
                    {
                        "scope": scope_name,
                        "pair_count": scope["pair_count"],
                        "pairs_per_job": scenario["pairs_per_job"],
                        "job_count": scenario["job_count"],
                        "startup_days": scenario["startup_time"]["days"],
                        "folding_days": scope["total_folding_time"]["days"],
                        "total_serial_days": scenario["total_serial_time"]["days"],
                        "total_serial_years": scenario["total_serial_time"][
                            "years_365_25_days"
                        ],
                        "startup_fraction": scenario["startup_fraction"],
                    }
                )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    length_by_symbol, length_audit = load_lengths()

    supported_per_pair_path = OUTPUT_DIR / "supported_pair_runtime_estimates.tsv.gz"
    per_pair_headers = [
        "edge_id",
        "node_a",
        "node_b",
        "length_a_aa",
        "length_b_aa",
        "combined_length_aa",
        "estimated_tokens",
        "estimated_folding_seconds",
        "estimated_pair_seconds_if_one_pair_per_job",
    ]
    with gzip.open(
        supported_per_pair_path, "wt", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=per_pair_headers,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        supported = evaluate_pair_file(
            PAIR_DIR / "supported_colocalization_edges.tsv.gz",
            length_by_symbol,
            writer,
        )

    scopes = {
        "supported_colocalization_pairs": supported,
        "conservative_alphafold_retained_pairs": evaluate_pair_file(
            PAIR_DIR / "alphafold_candidate_pairs.tsv.gz", length_by_symbol
        ),
        "all_621_protein_pairs": evaluate_pair_file(
            PAIR_DIR / "colocalization_protein_pairs.tsv.gz", length_by_symbol
        ),
    }
    for scope in scopes.values():
        add_batching_scenarios(scope)

    calibration = {
        "Calm1_Nme2": {
            "combined_length_aa": 149 + 152,
            "estimated_tokens": TOKENS_PER_RESIDUE * (149 + 152),
            "model_folding_seconds": (
                FOLDING_COEFFICIENT
                * (TOKENS_PER_RESIDUE * (149 + 152)) ** FOLDING_EXPONENT
            ),
            "user_reported_folding_seconds": 348,
        },
        "Camk2d_Pde4d": {
            "combined_length_aa": 499 + 747,
            "estimated_tokens": TOKENS_PER_RESIDUE * (499 + 747),
            "model_folding_seconds": (
                FOLDING_COEFFICIENT
                * (TOKENS_PER_RESIDUE * (499 + 747)) ** FOLDING_EXPONENT
            ),
            "user_reported_folding_seconds": 1336,
        },
    }
    summary = {
        "runtime_model": {
            "startup_seconds_per_job": STARTUP_SECONDS,
            "folding_seconds_formula": "1.29 * tokens^0.95",
            "tokens_per_combined_amino_acid_residue": TOKENS_PER_RESIDUE,
            "token_conversion_basis": (
                "Inferred from the two user-supplied examples; it reproduces "
                "their folding times to within approximately 0.2%."
            ),
        },
        "calibration_check": calibration,
        "protein_length_coverage": length_audit,
        "scopes": scopes,
        "interpretation": {
            "primary_scope": "supported_colocalization_pairs",
            "serial_time_definition": (
                "Total folding work plus 73 seconds per submitted batch/job."
            ),
            "parallel_time_caveat": (
                "Parallel wall times are ideal lower-bound estimates assuming "
                "perfect scheduling, no queueing, no failures, and one job per worker."
            ),
            "extrapolation_caveat": (
                "The fit has r^2=0.89 and is extrapolated to very large proteins; "
                "actual memory/token limits may make some long pairs infeasible."
            ),
        },
    }

    summary_path = OUTPUT_DIR / "runtime_estimate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_scenario_table(OUTPUT_DIR / "batching_scenarios.tsv", scopes)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
