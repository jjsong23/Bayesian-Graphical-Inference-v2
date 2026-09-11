#!/usr/bin/env python3
"""Integrate compartment adjacency into the Bayesian colocalization graph.

For HPA and COMPARTMENTS, each source obtains two similarities:

1. the original exact-overlap cosine similarity; and
2. expected compatibility p_i.T @ C @ p_j using L1-normalized location
   profiles and the scientist-reviewable native-location C matrix.

The adjacency-aware similarity is their elementwise maximum, so adjacency can
add but cannot erase exact-overlap evidence.  T_q backgrounds and the
complement-of-minimum likelihood are recomputed.  To keep this positive-only
evidence stream monotone, the final source likelihood is the maximum of the
original likelihood and the recomputed adjacency-aware likelihood.

The primary mpkCCD C matrix is identity because its five columns are
centrifugation fractions, not named compartments.  Its technical-proximity C
matrix is reported as a sensitivity analysis.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from build_colocalization_graph import (
    EPSILON,
    MINIMUM_LIKELIHOOD,
    NO_COLOCALIZATION_LIKELIHOOD,
    Q,
    build_compartments_profiles,
    load_compartments,
    posterior_from_log_bf,
    sha256_file,
    write_matrix,
)


PROJECT_DEFAULT = Path(__file__).resolve().parents[2]
NEUTRAL_TOLERANCE = 1e-12
FRACTIONS = ["1K", "4K", "17K", "200Kp", "200Ks"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_DEFAULT)
    return parser.parse_args()


def read_aligned_matrix(path: Path, labels: list[str]) -> np.ndarray:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.index.tolist() != labels or frame.columns.tolist() != labels:
        raise ValueError(f"Matrix label/order mismatch: {path}")
    return frame.to_numpy(dtype=float)


def read_c_matrix(path: Path) -> tuple[list[str], np.ndarray]:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.index.tolist() != frame.columns.tolist():
        raise ValueError(f"C matrix row/column mismatch: {path}")
    matrix = frame.to_numpy(dtype=float)
    if (
        not np.allclose(matrix, matrix.T)
        or not np.allclose(np.diag(matrix), 1.0)
        or np.any(matrix < 0)
        or np.any(matrix > 1)
    ):
        raise ValueError(f"Invalid C matrix: {path}")
    return frame.index.tolist(), matrix


def boolean_series(frame: pd.DataFrame, column: str) -> np.ndarray:
    return (
        frame[column]
        .astype(str)
        .str.strip()
        .str.casefold()
        .eq("true")
        .to_numpy()
    )


def score_similarity(
    similarity_observed: np.ndarray,
    observed_indices: np.ndarray,
    node_count: int,
) -> dict[str, object]:
    """Apply the per-node q75 complement-of-minimum procedure."""

    similarity_observed = np.asarray(similarity_observed, dtype=float)
    similarity_observed = np.clip(similarity_observed, 0.0, None)
    threshold_background = similarity_observed.copy()
    np.fill_diagonal(threshold_background, np.nan)
    raw_thresholds = np.nanquantile(threshold_background, Q, axis=1)
    thresholds = raw_thresholds.copy()
    methods: list[str] = []
    positive_sizes: list[int] = []

    for row_index, raw_threshold in enumerate(raw_thresholds):
        row = threshold_background[row_index]
        positive = row[np.isfinite(row) & (row > EPSILON)]
        positive_sizes.append(int(len(positive)))
        if np.isfinite(raw_threshold) and raw_threshold > EPSILON:
            methods.append("all_observed_nodes_q75")
        elif positive.size:
            thresholds[row_index] = float(np.quantile(positive, Q))
            methods.append("positive_overlap_q75_fallback")
        else:
            thresholds[row_index] = np.nan
            methods.append("no_positive_overlap_neutral")

    directed = np.full_like(similarity_observed, MINIMUM_LIKELIHOOD)
    for row_index, threshold in enumerate(thresholds):
        if not np.isfinite(threshold) or threshold <= EPSILON:
            continue
        ratio = similarity_observed[row_index] / threshold
        with np.errstate(over="ignore", invalid="ignore"):
            evidence = 1.0 - np.exp(-0.5 * np.square(ratio))
        directed[row_index] = np.maximum(MINIMUM_LIKELIHOOD, evidence)
    symmetric = (directed + directed.T) / 2.0

    similarity = np.zeros((node_count, node_count), dtype=float)
    likelihood = np.full(
        (node_count, node_count), MINIMUM_LIKELIHOOD, dtype=float
    )
    similarity[np.ix_(observed_indices, observed_indices)] = (
        similarity_observed
    )
    likelihood[np.ix_(observed_indices, observed_indices)] = symmetric
    np.fill_diagonal(similarity, 0.0)
    np.fill_diagonal(likelihood, 0.0)

    threshold_full = np.full(node_count, np.nan, dtype=float)
    raw_threshold_full = np.full(node_count, np.nan, dtype=float)
    method_full = np.full(node_count, "profile_not_observed", dtype=object)
    positive_size_full = np.zeros(node_count, dtype=int)
    threshold_full[observed_indices] = thresholds
    raw_threshold_full[observed_indices] = raw_thresholds
    method_full[observed_indices] = methods
    positive_size_full[observed_indices] = positive_sizes
    return {
        "similarity": similarity,
        "likelihood": likelihood,
        "threshold": threshold_full,
        "raw_threshold": raw_threshold_full,
        "threshold_method": method_full,
        "positive_background_size": positive_size_full,
    }


def adjacency_aware_similarity(
    profiles: np.ndarray,
    c_matrix: np.ndarray,
) -> dict[str, object]:
    """Return exact, expected-C, and positive-only combined similarities."""

    profiles = np.asarray(profiles, dtype=float)
    observed_mask = profiles.sum(axis=1) > 0
    observed_indices = np.flatnonzero(observed_mask)
    observed = profiles[observed_mask]
    l2 = observed / np.linalg.norm(observed, axis=1, keepdims=True)
    l1 = observed / observed.sum(axis=1, keepdims=True)
    exact = np.clip(l2 @ l2.T, 0.0, 1.0)
    expected = np.clip((l1 @ c_matrix) @ l1.T, 0.0, 1.0)
    adjacency_aware = np.maximum(exact, expected)
    np.fill_diagonal(exact, 1.0)
    np.fill_diagonal(expected, 1.0)
    np.fill_diagonal(adjacency_aware, 1.0)
    return {
        "observed_mask": observed_mask,
        "observed_indices": observed_indices,
        "exact_observed": exact,
        "expected_observed": expected,
        "adjacency_observed": adjacency_aware,
    }


def full_similarity(
    observed_matrix: np.ndarray,
    observed_indices: np.ndarray,
    node_count: int,
) -> np.ndarray:
    result = np.zeros((node_count, node_count), dtype=float)
    result[np.ix_(observed_indices, observed_indices)] = observed_matrix
    np.fill_diagonal(result, 0.0)
    return result


def log_bf_from_likelihood(likelihood: np.ndarray) -> np.ndarray:
    adjusted = np.maximum(likelihood, NO_COLOCALIZATION_LIKELIHOOD)
    log_bf = np.log(adjusted / NO_COLOCALIZATION_LIKELIHOOD)
    np.fill_diagonal(log_bf, 0.0)
    return log_bf


def safe_hpa_column(location: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9]+", "_", location).strip("_").lower()
    return f"location__{safe}"


def matrix_qc(matrix: np.ndarray) -> dict[str, object]:
    off_diagonal = matrix[~np.eye(len(matrix), dtype=bool)]
    return {
        "shape": list(matrix.shape),
        "symmetric": bool(np.allclose(matrix, matrix.T, atol=1e-12)),
        "diagonal_all_zero": bool(np.allclose(np.diag(matrix), 0.0)),
        "off_diagonal_minimum": float(off_diagonal.min()),
        "off_diagonal_maximum": float(off_diagonal.max()),
    }


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    hpa_profile_path = (
        project
        / "data/edge_characterization/localization/hpa/v25.1/processed/"
        "node_hpa_localization_profiles.tsv"
    )
    compartments_raw_path = (
        project
        / "data/colocalization/compartments/raw/"
        "mouse_compartment_knowledge_filtered.tsv"
    )
    string_mapping_path = (
        project
        / "data/edge_characterization/string/v12.0/processed/"
        "node_to_string_mapping.tsv"
    )
    mpk_profile_path = (
        project
        / "data/edge_characterization/localization/processed/"
        "node_localization_profiles.tsv"
    )
    c_dir = project / "data/colocalization/compatibility_matrices"
    c_hpa_path = c_dir / "C_hpa_locations.tsv"
    c_compartments_path = c_dir / "C_compartments_go_locations.tsv"
    c_mpk_primary_path = c_dir / "C_mpkccd_fractions_primary.tsv"
    c_mpk_sensitivity_path = c_dir / "C_mpkccd_fractions_sensitivity.tsv"
    old_result_dir = project / "results/colocalization"
    old_final_path = old_result_dir / "colocalization_adjacency_matrix.tsv"
    old_hpa_likelihood_path = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa/"
        "hpa_localization_likelihood_matrix.tsv"
    )
    old_hpa_similarity_path = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa/"
        "hpa_localization_similarity_matrix.tsv"
    )
    old_comp_likelihood_path = (
        old_result_dir / "compartments_likelihood_matrix.tsv"
    )
    old_comp_similarity_path = (
        old_result_dir / "compartments_similarity_matrix.tsv"
    )
    old_mpk_likelihood_path = (
        project
        / "results/edge_characterization/localization/"
        "localization_evidence_matrix.tsv"
    )
    output_dir = project / "results/colocalization_adjacency_aware"
    output_dir.mkdir(parents=True, exist_ok=True)

    required = [
        universe_path,
        hpa_profile_path,
        compartments_raw_path,
        string_mapping_path,
        mpk_profile_path,
        c_hpa_path,
        c_compartments_path,
        c_mpk_primary_path,
        c_mpk_sensitivity_path,
        old_final_path,
        old_hpa_likelihood_path,
        old_hpa_similarity_path,
        old_comp_likelihood_path,
        old_comp_similarity_path,
        old_mpk_likelihood_path,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    symbols = universe["symbol"].astype(str).tolist()
    if len(symbols) != 868 or universe["symbol"].duplicated().any():
        raise ValueError("Expected 868 unique node symbols")
    protein_mask = universe["node_type"].astype(str).eq("protein").to_numpy()
    protein_indices = np.flatnonzero(protein_mask)
    protein_symbols = np.asarray(symbols)[protein_indices].tolist()
    if len(protein_symbols) != 848:
        raise ValueError("Expected 848 protein nodes")
    node_count = len(symbols)

    # HPA: align the 49 native labels to their exact processed profile columns.
    hpa_labels, c_hpa = read_c_matrix(c_hpa_path)
    hpa_profiles_frame = pd.read_csv(
        hpa_profile_path, sep="\t", dtype=str
    ).fillna("")
    if hpa_profiles_frame["symbol"].astype(str).tolist() != symbols:
        raise ValueError("HPA profile order differs from node universe")
    hpa_columns = [safe_hpa_column(label) for label in hpa_labels]
    missing_hpa = [
        column for column in hpa_columns if column not in hpa_profiles_frame
    ]
    if missing_hpa:
        raise ValueError(f"Missing HPA profile columns: {missing_hpa}")
    hpa_profiles = (
        hpa_profiles_frame[hpa_columns]
        .apply(pd.to_numeric, errors="raise")
        .to_numpy(dtype=float)
    )
    hpa_similarity = adjacency_aware_similarity(hpa_profiles, c_hpa)
    hpa_indices = np.asarray(hpa_similarity["observed_indices"], dtype=int)
    hpa_score = score_similarity(
        np.asarray(hpa_similarity["adjacency_observed"]),
        hpa_indices,
        node_count,
    )
    hpa_exact_full = full_similarity(
        np.asarray(hpa_similarity["exact_observed"]), hpa_indices, node_count
    )
    hpa_expected_full = full_similarity(
        np.asarray(hpa_similarity["expected_observed"]), hpa_indices, node_count
    )
    old_hpa_similarity = read_aligned_matrix(
        old_hpa_similarity_path, symbols
    )
    old_hpa_likelihood = read_aligned_matrix(
        old_hpa_likelihood_path, symbols
    )
    hpa_exact_reconstruction_error = float(
        np.max(np.abs(hpa_exact_full - old_hpa_similarity))
    )
    hpa_likelihood = np.maximum(
        np.asarray(hpa_score["likelihood"]), old_hpa_likelihood
    )
    np.fill_diagonal(hpa_likelihood, 0.0)

    # COMPARTMENTS: rebuild the same score>=4 node profiles, then apply C.
    compartments = load_compartments(compartments_raw_path)
    string_mapping = pd.read_csv(
        string_mapping_path, sep="\t", dtype=str
    ).fillna("")
    (
        compartments_profiles,
        compartments_mapping,
        _evidence_audit,
        compartments_dictionary,
    ) = build_compartments_profiles(
        universe,
        compartments,
        string_mapping,
        minimum_score=4,
        include_audit=False,
    )
    comp_labels, c_comp = read_c_matrix(c_compartments_path)
    comp_go_ids = [label.split(" | ", 1)[0] for label in comp_labels]
    expected_go_ids = compartments_dictionary["go_id"].astype(str).tolist()
    if comp_go_ids != expected_go_ids:
        raise ValueError("COMPARTMENTS C order differs from rebuilt profiles")
    comp_similarity = adjacency_aware_similarity(
        compartments_profiles, c_comp
    )
    comp_indices = np.asarray(comp_similarity["observed_indices"], dtype=int)
    comp_score = score_similarity(
        np.asarray(comp_similarity["adjacency_observed"]),
        comp_indices,
        node_count,
    )
    comp_exact_full = full_similarity(
        np.asarray(comp_similarity["exact_observed"]), comp_indices, node_count
    )
    comp_expected_full = full_similarity(
        np.asarray(comp_similarity["expected_observed"]),
        comp_indices,
        node_count,
    )
    old_comp_similarity = read_aligned_matrix(
        old_comp_similarity_path, symbols
    )
    old_comp_likelihood = read_aligned_matrix(
        old_comp_likelihood_path, symbols
    )
    comp_exact_reconstruction_error = float(
        np.max(np.abs(comp_exact_full - old_comp_similarity))
    )
    comp_likelihood = np.maximum(
        np.asarray(comp_score["likelihood"]), old_comp_likelihood
    )
    np.fill_diagonal(comp_likelihood, 0.0)

    # mpkCCD: primary C is identity; technical proximity is sensitivity only.
    mpk_primary_labels, c_mpk_primary = read_c_matrix(c_mpk_primary_path)
    mpk_sensitivity_labels, c_mpk_sensitivity = read_c_matrix(
        c_mpk_sensitivity_path
    )
    if mpk_primary_labels != FRACTIONS or mpk_sensitivity_labels != FRACTIONS:
        raise ValueError("mpkCCD C labels do not match fraction columns")
    mpk_frame = pd.read_csv(mpk_profile_path, sep="\t", dtype=str).fillna("")
    if mpk_frame["symbol"].astype(str).tolist() != symbols:
        raise ValueError("mpkCCD profile order differs from node universe")
    mpk_observed = boolean_series(mpk_frame, "localization_observed")
    mpk_indices = np.flatnonzero(mpk_observed)
    mpk_profiles = (
        mpk_frame[FRACTIONS]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    mpk_observed_profiles = mpk_profiles[mpk_observed]
    mpk_primary_similarity_observed = (
        mpk_observed_profiles @ c_mpk_primary @ mpk_observed_profiles.T
    )
    mpk_sensitivity_similarity_observed = (
        mpk_observed_profiles @ c_mpk_sensitivity @ mpk_observed_profiles.T
    )
    mpk_primary_score = score_similarity(
        mpk_primary_similarity_observed, mpk_indices, node_count
    )
    mpk_sensitivity_score = score_similarity(
        mpk_sensitivity_similarity_observed, mpk_indices, node_count
    )
    old_mpk_likelihood = read_aligned_matrix(
        old_mpk_likelihood_path, symbols
    )
    mpk_primary_likelihood = np.maximum(
        np.asarray(mpk_primary_score["likelihood"]), old_mpk_likelihood
    )
    mpk_sensitivity_likelihood = np.maximum(
        np.asarray(mpk_sensitivity_score["likelihood"]),
        mpk_primary_likelihood,
    )
    np.fill_diagonal(mpk_primary_likelihood, 0.0)
    np.fill_diagonal(mpk_sensitivity_likelihood, 0.0)
    mpk_primary_reconstruction_error = float(
        np.max(np.abs(mpk_primary_likelihood - old_mpk_likelihood))
    )

    # Primary integration.
    hpa_log_bf = log_bf_from_likelihood(hpa_likelihood)
    comp_log_bf = log_bf_from_likelihood(comp_likelihood)
    mpk_log_bf = log_bf_from_likelihood(mpk_primary_likelihood)
    combined_log_bf = hpa_log_bf + comp_log_bf + mpk_log_bf
    combined_posterior = posterior_from_log_bf(combined_log_bf)
    np.fill_diagonal(combined_posterior, 0.0)

    # Sensitivity: only mpkCCD's optional technical-proximity C changes.
    mpk_sensitivity_log_bf = log_bf_from_likelihood(
        mpk_sensitivity_likelihood
    )
    sensitivity_log_bf = (
        hpa_log_bf + comp_log_bf + mpk_sensitivity_log_bf
    )
    sensitivity_posterior = posterior_from_log_bf(sensitivity_log_bf)
    np.fill_diagonal(sensitivity_posterior, 0.0)

    # Reconstruct and compare the pre-C graph.
    old_log_bf = (
        log_bf_from_likelihood(old_hpa_likelihood)
        + log_bf_from_likelihood(old_comp_likelihood)
        + log_bf_from_likelihood(old_mpk_likelihood)
    )
    old_reconstructed = posterior_from_log_bf(old_log_bf)
    np.fill_diagonal(old_reconstructed, 0.0)
    old_protein = read_aligned_matrix(old_final_path, protein_symbols)
    old_reconstruction_error = float(
        np.max(
            np.abs(
                old_reconstructed[np.ix_(protein_indices, protein_indices)]
                - old_protein
            )
        )
    )

    protein_posterior = combined_posterior[
        np.ix_(protein_indices, protein_indices)
    ]
    protein_log_bf = combined_log_bf[np.ix_(protein_indices, protein_indices)]
    protein_sensitivity = sensitivity_posterior[
        np.ix_(protein_indices, protein_indices)
    ]
    upper_i, upper_j = np.triu_indices(len(protein_symbols), k=1)
    full_upper_i, full_upper_j = np.triu_indices(node_count, k=1)
    posterior_upper = protein_posterior[upper_i, upper_j]
    log_bf_upper = protein_log_bf[upper_i, upper_j]
    old_upper = old_protein[upper_i, upper_j]
    sensitivity_upper = protein_sensitivity[upper_i, upper_j]
    neutral = np.abs(log_bf_upper) <= NEUTRAL_TOLERANCE
    old_neutral = np.abs(old_upper - 0.5) <= NEUTRAL_TOLERANCE
    sensitivity_neutral = (
        np.abs(sensitivity_upper - 0.5) <= NEUTRAL_TOLERANCE
    )
    full_neutral = (
        np.abs(combined_log_bf[full_upper_i, full_upper_j])
        <= NEUTRAL_TOLERANCE
    )

    # Raw source compatibility and conservative post-C disjointness.
    protein_global_i = protein_indices[upper_i]
    protein_global_j = protein_indices[upper_j]
    mpk_raw = np.asarray(mpk_primary_score["similarity"])[
        protein_global_i, protein_global_j
    ]
    hpa_raw = np.asarray(hpa_score["similarity"])[
        protein_global_i, protein_global_j
    ]
    comp_raw = np.asarray(comp_score["similarity"])[
        protein_global_i, protein_global_j
    ]
    raw_all_zero = (
        np.maximum.reduce([mpk_raw, hpa_raw, comp_raw])
        <= NEUTRAL_TOLERANCE
    )
    hpa_pair_observed = np.outer(
        np.asarray(hpa_similarity["observed_mask"], dtype=bool),
        np.asarray(hpa_similarity["observed_mask"], dtype=bool),
    )[protein_global_i, protein_global_j]
    comp_pair_observed = np.outer(
        np.asarray(comp_similarity["observed_mask"], dtype=bool),
        np.asarray(comp_similarity["observed_mask"], dtype=bool),
    )[protein_global_i, protein_global_j]
    hpa_disjoint = hpa_pair_observed & (
        hpa_raw <= NEUTRAL_TOLERANCE
    )
    comp_disjoint = comp_pair_observed & (
        comp_raw <= NEUTRAL_TOLERANCE
    )
    any_adjacency_compatible = (
        (hpa_pair_observed & (hpa_raw > NEUTRAL_TOLERANCE))
        | (comp_pair_observed & (comp_raw > NEUTRAL_TOLERANCE))
    )
    recommended_exclude = (
        neutral
        & ~any_adjacency_compatible
        & (hpa_disjoint | comp_disjoint)
    )

    # Output matrices.
    matrix_outputs: list[Path] = []

    def output_matrix(
        filename: str,
        matrix: np.ndarray,
        labels: list[str] = symbols,
    ) -> Path:
        path = output_dir / filename
        write_matrix(path, labels, matrix)
        matrix_outputs.append(path)
        return path

    final_path = output_matrix(
        "colocalization_adjacency_matrix.tsv",
        protein_posterior,
        protein_symbols,
    )
    output_matrix(
        "colocalization_log_bayes_factor_matrix.tsv",
        protein_log_bf,
        protein_symbols,
    )
    output_matrix(
        "colocalization_adjacency_matrix_all_nodes.tsv",
        combined_posterior,
    )
    sensitivity_path = output_matrix(
        "colocalization_adjacency_matrix_mpkccd_sensitivity.tsv",
        protein_sensitivity,
        protein_symbols,
    )
    output_matrix(
        "hpa_exact_similarity_matrix.tsv", hpa_exact_full
    )
    output_matrix(
        "hpa_expected_C_compatibility_matrix.tsv", hpa_expected_full
    )
    output_matrix(
        "hpa_adjacency_aware_similarity_matrix.tsv",
        np.asarray(hpa_score["similarity"]),
    )
    output_matrix(
        "hpa_adjacency_aware_likelihood_matrix.tsv", hpa_likelihood
    )
    output_matrix(
        "compartments_exact_similarity_matrix.tsv", comp_exact_full
    )
    output_matrix(
        "compartments_expected_C_compatibility_matrix.tsv",
        comp_expected_full,
    )
    output_matrix(
        "compartments_adjacency_aware_similarity_matrix.tsv",
        np.asarray(comp_score["similarity"]),
    )
    output_matrix(
        "compartments_adjacency_aware_likelihood_matrix.tsv",
        comp_likelihood,
    )
    output_matrix(
        "mpkccd_primary_similarity_matrix.tsv",
        np.asarray(mpk_primary_score["similarity"]),
    )
    output_matrix(
        "mpkccd_primary_likelihood_matrix.tsv", mpk_primary_likelihood
    )
    output_matrix(
        "mpkccd_sensitivity_likelihood_matrix.tsv",
        mpk_sensitivity_likelihood,
    )

    threshold_audit = pd.DataFrame(
        {
            "symbol": symbols,
            "node_type": universe["node_type"].astype(str),
            "hpa_profile_observed": np.asarray(
                hpa_similarity["observed_mask"], dtype=bool
            ),
            "hpa_raw_Tq_q75": hpa_score["raw_threshold"],
            "hpa_Tq_q75": hpa_score["threshold"],
            "hpa_Tq_method": hpa_score["threshold_method"],
            "hpa_positive_background_size": hpa_score[
                "positive_background_size"
            ],
            "compartments_profile_observed": np.asarray(
                comp_similarity["observed_mask"], dtype=bool
            ),
            "compartments_raw_Tq_q75": comp_score["raw_threshold"],
            "compartments_Tq_q75": comp_score["threshold"],
            "compartments_Tq_method": comp_score["threshold_method"],
            "compartments_positive_background_size": comp_score[
                "positive_background_size"
            ],
            "mpkccd_profile_observed": mpk_observed,
            "mpkccd_primary_Tq_q75": mpk_primary_score["threshold"],
            "mpkccd_primary_Tq_method": mpk_primary_score[
                "threshold_method"
            ],
            "mpkccd_sensitivity_Tq_q75": mpk_sensitivity_score[
                "threshold"
            ],
            "mpkccd_sensitivity_Tq_method": mpk_sensitivity_score[
                "threshold_method"
            ],
        }
    )
    threshold_path = output_dir / "adjacency_aware_Tq_audit.tsv"
    threshold_audit.to_csv(threshold_path, sep="\t", index=False)

    edge_frame = pd.DataFrame(
        {
            "edge_id": [
                f"{a}|{b}"
                for a, b in zip(
                    np.asarray(protein_symbols)[upper_i],
                    np.asarray(protein_symbols)[upper_j],
                    strict=True,
                )
            ],
            "node_a": np.asarray(protein_symbols)[upper_i],
            "node_b": np.asarray(protein_symbols)[upper_j],
            "mpkccd_primary_similarity": mpk_raw,
            "mpkccd_likelihood": mpk_primary_likelihood[
                protein_global_i, protein_global_j
            ],
            "hpa_exact_similarity": hpa_exact_full[
                protein_global_i, protein_global_j
            ],
            "hpa_expected_C_compatibility": hpa_expected_full[
                protein_global_i, protein_global_j
            ],
            "hpa_adjacency_aware_similarity": hpa_raw,
            "hpa_likelihood": hpa_likelihood[
                protein_global_i, protein_global_j
            ],
            "compartments_exact_similarity": comp_exact_full[
                protein_global_i, protein_global_j
            ],
            "compartments_expected_C_compatibility": comp_expected_full[
                protein_global_i, protein_global_j
            ],
            "compartments_adjacency_aware_similarity": comp_raw,
            "compartments_likelihood": comp_likelihood[
                protein_global_i, protein_global_j
            ],
            "combined_log_bayes_factor": log_bf_upper,
            "colocalization_posterior": posterior_upper,
            "at_neutral_baseline": neutral,
            "all_raw_source_compatibilities_zero": raw_all_zero,
            "adjacency_rescued_from_old_neutral": old_neutral & ~neutral,
            "hpa_adjacency_disjoint": hpa_disjoint,
            "compartments_adjacency_disjoint": comp_disjoint,
            "recommended_exclude_from_alphafold": recommended_exclude,
            "recommended_retain_for_alphafold": ~recommended_exclude,
            "mpkccd_sensitivity_posterior": sensitivity_upper,
        }
    )
    all_pairs_path = output_dir / "colocalization_protein_pairs.tsv.gz"
    with gzip.open(
        all_pairs_path, "wt", encoding="utf-8", newline=""
    ) as handle:
        edge_frame.to_csv(
            handle, sep="\t", index=False, float_format="%.9g"
        )
    supported_path = output_dir / "supported_colocalization_edges.tsv.gz"
    with gzip.open(
        supported_path, "wt", encoding="utf-8", newline=""
    ) as handle:
        edge_frame.loc[~neutral].to_csv(
            handle, sep="\t", index=False, float_format="%.9g"
        )
    candidate_path = output_dir / "alphafold_candidate_pairs.tsv.gz"
    with gzip.open(
        candidate_path, "wt", encoding="utf-8", newline=""
    ) as handle:
        edge_frame.loc[~recommended_exclude].to_csv(
            handle, sep="\t", index=False, float_format="%.9g"
        )

    validation = {
        "hpa_C_aligned_to_49_profile_columns": len(hpa_columns) == 49,
        "compartments_C_aligned_to_profile_terms": (
            comp_go_ids == expected_go_ids
        ),
        "mpkccd_C_aligned_to_fractions": (
            mpk_primary_labels == FRACTIONS
            and mpk_sensitivity_labels == FRACTIONS
        ),
        "hpa_exact_similarity_reconstructs_previous": (
            hpa_exact_reconstruction_error <= 1e-8
        ),
        "compartments_exact_similarity_reconstructs_previous": (
            comp_exact_reconstruction_error <= 1e-8
        ),
        "mpkccd_primary_identity_reconstructs_previous": (
            mpk_primary_reconstruction_error <= 1e-8
        ),
        "old_combined_graph_reconstructed": (
            old_reconstruction_error <= 1e-8
        ),
        "new_graph_never_below_old_graph": bool(
            np.all(posterior_upper + 1e-8 >= old_upper)
        ),
        "new_minus_stored_old_minimum": float(
            np.min(posterior_upper - old_upper)
        ),
        "primary_graph": matrix_qc(protein_posterior),
        "full_graph": matrix_qc(combined_posterior),
        "sensitivity_graph": matrix_qc(protein_sensitivity),
        "edge_ids_unique": bool(edge_frame["edge_id"].is_unique),
        "pair_partition": bool(
            int(neutral.sum()) + int((~neutral).sum()) == len(edge_frame)
        ),
    }
    if not all(
        value
        for key, value in validation.items()
        if isinstance(value, bool)
    ):
        raise AssertionError(f"Validation failed: {validation}")

    summary = {
        "purpose": "Adjacency-aware three-stream Bayesian colocalization graph.",
        "protein_node_count": len(protein_symbols),
        "full_node_count": len(symbols),
        "possible_undirected_protein_pairs": len(edge_frame),
        "possible_undirected_all_node_pairs": len(full_upper_i),
        "neutral_baseline_definition": (
            f"abs(combined_log_BF) <= {NEUTRAL_TOLERANCE}; "
            "equivalent to posterior 0.5 in this positive-or-neutral model."
        ),
        "near_zero_baseline_protein_pairs": int(neutral.sum()),
        "near_zero_baseline_all_node_pairs": int(full_neutral.sum()),
        "non_neutral_supported_protein_pairs": int((~neutral).sum()),
        "raw_all_source_compatibilities_zero_protein_pairs": int(
            raw_all_zero.sum()
        ),
        "previous_neutral_baseline_protein_pairs": int(old_neutral.sum()),
        "protein_pairs_rescued_from_neutral_by_adjacency": int(
            (old_neutral & ~neutral).sum()
        ),
        "neutral_pair_reduction": int(old_neutral.sum() - neutral.sum()),
        "neutral_pair_reduction_fraction": float(
            (old_neutral.sum() - neutral.sum()) / old_neutral.sum()
        ),
        "source_edges_strengthened_by_C": {
            "hpa": int(
                (
                    hpa_likelihood[protein_global_i, protein_global_j]
                    > old_hpa_likelihood[
                        protein_global_i, protein_global_j
                    ]
                    + 1e-8
                ).sum()
            ),
            "compartments": int(
                (
                    comp_likelihood[protein_global_i, protein_global_j]
                    > old_comp_likelihood[
                        protein_global_i, protein_global_j
                    ]
                    + 1e-8
                ).sum()
            ),
            "mpkccd_primary": int(
                (
                    mpk_primary_likelihood[
                        protein_global_i, protein_global_j
                    ]
                    > old_mpk_likelihood[
                        protein_global_i, protein_global_j
                    ]
                    + 1e-8
                ).sum()
            ),
        },
        "recommended_alphafold_pairs_retained": int(
            (~recommended_exclude).sum()
        ),
        "recommended_alphafold_pairs_excluded": int(
            recommended_exclude.sum()
        ),
        "sensitivity_mpkccd_technical_proximity": {
            "near_zero_baseline_protein_pairs": int(
                sensitivity_neutral.sum()
            ),
            "non_neutral_supported_protein_pairs": int(
                (~sensitivity_neutral).sum()
            ),
            "pairs_rescued_relative_to_primary": int(
                (neutral & ~sensitivity_neutral).sum()
            ),
        },
        "adjacency_similarity_rule": (
            "max(original exact-overlap cosine, "
            "L1_profile_i.T @ C @ L1_profile_j)"
        ),
        "source_likelihood_rule": (
            "max(original source likelihood, adjacency-aware likelihood "
            "after recomputing node-specific q75 Tq backgrounds)"
        ),
        "combination_rule": (
            "logit(P_colocalized)=log(BF_mpkCCD)+log(BF_HPA)+"
            "log(BF_COMPARTMENTS), with prior 0.5"
        ),
        "mpkccd_primary_rule": (
            "Identity C; no biological adjacency inferred between "
            "centrifugation fractions."
        ),
        "interpretation_caveat": (
            "C weights are provisional compatibility priors, not "
            "probabilities or direct PPI evidence."
        ),
        "reconstruction_errors": {
            "hpa_exact_similarity": hpa_exact_reconstruction_error,
            "compartments_exact_similarity": (
                comp_exact_reconstruction_error
            ),
            "mpkccd_primary_likelihood": (
                mpk_primary_reconstruction_error
            ),
            "old_combined_graph": old_reconstruction_error,
        },
        "validation": validation,
        "input_sha256": {
            path.name: sha256_file(path)
            for path in required
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    readme = f"""ADJACENCY-AWARE BAYESIAN COLOCALIZATION

This analysis integrates the scientist-reviewable C matrices into the existing
mpkCCD, HPA, and COMPARTMENTS localization evidence streams.

PRIMARY RESULT
- Protein pairs: {len(edge_frame):,}
- Neutral/near-zero log-BF pairs: {int(neutral.sum()):,}
- Supported pairs: {int((~neutral).sum()):,}
- Previous neutral pairs: {int(old_neutral.sum()):,}
- Pairs rescued from neutral by adjacency: {int((old_neutral & ~neutral).sum()):,}

NEUTRAL DEFINITION
abs(combined log-BF) <= {NEUTRAL_TOLERANCE}, which is posterior 0.5.

SOURCE SIMILARITY
For HPA and COMPARTMENTS:
    exact = cosine similarity of the original location profiles
    expected_C = L1(profile_i)' C L1(profile_j)
    adjacency-aware = max(exact, expected_C)

Tq is recomputed for every observed node. The final likelihood for a source is
the maximum of the previous exact-location likelihood and the recomputed
adjacency-aware likelihood, so adding adjacency cannot erase previous support.

mpkCCD uses the identity C matrix in the primary analysis. The optional
fraction-proximity C is provided only as a sensitivity analysis.

KEY FILES
- colocalization_adjacency_matrix.tsv
- colocalization_adjacency_matrix_all_nodes.tsv
- colocalization_adjacency_matrix_mpkccd_sensitivity.tsv
- colocalization_protein_pairs.tsv.gz
- supported_colocalization_edges.tsv.gz
- alphafold_candidate_pairs.tsv.gz
- adjacency_aware_Tq_audit.tsv
- analysis_summary.json

IMPORTANT
C weights encode provisional physical compatibility. They are not direct
probabilities of interaction and should be used only after scientific review.
"""
    readme_path = output_dir / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")

    summary["output_sha256"] = {
        path.name: sha256_file(path)
        for path in [
            final_path,
            sensitivity_path,
            threshold_path,
            all_pairs_path,
            supported_path,
            candidate_path,
            readme_path,
            *matrix_outputs,
        ]
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
