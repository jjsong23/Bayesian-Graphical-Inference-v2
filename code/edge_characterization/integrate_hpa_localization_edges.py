#!/usr/bin/env python3
"""Integrate Human Protein Atlas localization into the undirected edge model.

Mouse protein nodes are mapped to human genes through MGI's stringent
one-to-one protein-coding orthology report. HPA Enhanced, Supported, and
Approved locations form a binary 49-location profile; Uncertain annotations
are audited but excluded from positive evidence. Profiles are unit-normalized,
so their dot product is cosine similarity and is not inflated merely by the
number of locations assigned to a protein.

The similarity is passed through the same per-node q75,
complement-of-minimum, reciprocal-symmetrization procedure used by the earlier
mpkCCD fractionation-localization analysis. Missing HPA evidence is neutral.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


PROJECT_DEFAULT = Path(__file__).resolve().parents[2]
CODE_DIR = PROJECT_DEFAULT / "code"
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from bayes_factors import binary_edge_bayes_update  # noqa: E402


HPA_VERSION = "25.1"
HPA_ENSEMBL_VERSION = "109"
Q = 0.75
MINIMUM_LIKELIHOOD = 0.5
NO_EDGE_LIKELIHOOD = 0.5
PRIMARY_TIERS = ("Enhanced", "Supported", "Approved")
HIGH_CONFIDENCE_TIERS = ("Enhanced", "Supported")
HPA_DATA_URL = (
    "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip"
)
HPA_PAGE_URL = "https://www.proteinatlas.org/humanproteome/subcellular"
HPA_METHOD_URL = (
    "https://www.proteinatlas.org/humanproteome/subcellular/method"
)
MGI_ORTHOLOGY_URL = (
    "https://www.informatics.jax.org/downloads/reports/HOM_ProteinCoding.rpt"
)
MGI_REPORT_PAGE = "https://www.informatics.jax.org/downloads/reports/index.html"

MGI_COLUMNS = [
    "mgi_id",
    "mouse_symbol",
    "mouse_entrez_gene_id",
    "hgnc_id",
    "human_symbol",
    "human_entrez_gene_id",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_DEFAULT,
        help=f"Project directory (default: {PROJECT_DEFAULT})",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_matrix(
    path: Path,
    symbols: list[str],
    matrix: np.ndarray,
    *,
    float_format: str = "%.9g",
) -> None:
    frame = pd.DataFrame(matrix, index=symbols, columns=symbols)
    frame.index.name = "symbol"
    frame.to_csv(path, sep="\t", float_format=float_format)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "-"} else text


def split_locations(value: object) -> list[str]:
    return [
        location.strip()
        for location in clean(value).split(";")
        if location.strip()
    ]


def load_hpa(path: Path) -> pd.DataFrame:
    with ZipFile(path) as archive:
        names = archive.namelist()
        if names != ["subcellular_location.tsv"]:
            raise ValueError(f"Unexpected HPA archive members: {names}")
        with archive.open(names[0]) as handle:
            frame = pd.read_csv(handle, sep="\t", dtype=str).fillna("")
    expected = {
        "Gene",
        "Gene name",
        "Reliability",
        "Main location",
        "Additional location",
        "Extracellular location",
        "Enhanced",
        "Supported",
        "Approved",
        "Uncertain",
        "Single-cell variation intensity",
        "Single-cell variation spatial",
        "Cell cycle dependency",
        "GO id",
    }
    if set(frame.columns) != expected:
        raise ValueError(f"Unexpected HPA columns: {frame.columns.tolist()}")
    if frame["Gene"].duplicated().any():
        raise ValueError("HPA table contains duplicate Ensembl gene identifiers")
    return frame


def load_one_to_one_orthology(
    path: Path,
) -> tuple[pd.DataFrame, dict[str, dict[str, str]], set[str]]:
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=MGI_COLUMNS,
        dtype=str,
    ).fillna("")
    if raw.shape[1] != len(MGI_COLUMNS):
        raise ValueError(f"Unexpected MGI report shape: {raw.shape}")

    # The report can contain technical duplicate rows with different Entrez
    # identifiers but the same asserted mouse/human pair. Collapse those.
    pair_level = raw.sort_values(MGI_COLUMNS, kind="stable").drop_duplicates(
        ["mouse_symbol", "human_symbol"],
        keep="first",
    )
    human_count = pair_level.groupby("mouse_symbol")["human_symbol"].nunique()
    ambiguous_mouse = set(human_count[human_count.ne(1)].index)
    unique = pair_level.loc[
        ~pair_level["mouse_symbol"].isin(ambiguous_mouse)
    ].drop_duplicates("mouse_symbol", keep="first")
    lookup = {
        str(row["mouse_symbol"]): {
            column: str(row[column])
            for column in MGI_COLUMNS
        }
        for row in unique.to_dict(orient="records")
    }
    return raw, lookup, ambiguous_mouse


def extract_location_dictionary(hpa: pd.DataFrame) -> pd.DataFrame:
    location_to_go: dict[str, set[str]] = defaultdict(set)
    location_to_gene_count: dict[str, set[str]] = defaultdict(set)
    pattern = re.compile(r"^(.*?)\s+\((GO:\d+)\)$")
    for row in hpa.to_dict(orient="records"):
        gene = str(row["Gene"])
        for entry in split_locations(row["GO id"]):
            match = pattern.match(entry)
            if match:
                location_to_go[match.group(1)].add(match.group(2))
            else:
                location_to_go[entry]
        for tier in (*PRIMARY_TIERS, "Uncertain"):
            for location in split_locations(row[tier]):
                location_to_gene_count[location].add(gene)

    primary_locations = {
        location
        for tier in PRIMARY_TIERS
        for value in hpa[tier]
        for location in split_locations(value)
    }
    uncertain_locations = {
        location
        for value in hpa["Uncertain"]
        for location in split_locations(value)
    }
    all_locations = sorted(primary_locations.union(uncertain_locations))
    return pd.DataFrame(
        [
            {
                "location": location,
                "go_ids": ";".join(sorted(location_to_go.get(location, set()))),
                "hpa_gene_count_any_tier": len(
                    location_to_gene_count.get(location, set())
                ),
                "included_in_primary_profile": location in primary_locations,
                "included_in_high_confidence_profile": any(
                    location in {
                        candidate
                        for value in hpa[tier]
                        for candidate in split_locations(value)
                    }
                    for tier in HIGH_CONFIDENCE_TIERS
                ),
            }
            for location in all_locations
        ]
    )


def build_profiles(
    universe: pd.DataFrame,
    hpa: pd.DataFrame,
    orthology_lookup: dict[str, dict[str, str]],
    ambiguous_mouse: set[str],
    locations: list[str],
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    location_to_index = {
        location: index for index, location in enumerate(locations)
    }
    hpa_name_counts = hpa["Gene name"].value_counts()
    hpa_unique = hpa.loc[
        hpa["Gene name"].map(hpa_name_counts).eq(1)
    ].set_index("Gene name")
    ambiguous_hpa_names = set(hpa_name_counts[hpa_name_counts.gt(1)].index)

    primary = np.zeros((len(universe), len(locations)), dtype=np.uint8)
    high_confidence = np.zeros_like(primary)
    audit_rows: list[dict[str, object]] = []

    for index, row in enumerate(universe.to_dict(orient="records")):
        symbol = str(row["symbol"])
        node_type = str(row["node_type"])
        orthology = orthology_lookup.get(symbol, {})
        human_symbol = clean(orthology.get("human_symbol", ""))
        hpa_row: pd.Series | None = None
        status = ""

        if node_type != "protein":
            status = "excluded_nonprotein_node"
        elif symbol in ambiguous_mouse:
            status = "ambiguous_mgi_human_orthology"
        elif not human_symbol:
            status = "no_mgi_one_to_one_human_ortholog"
        elif human_symbol in ambiguous_hpa_names:
            status = "ambiguous_hpa_gene_symbol"
        elif human_symbol not in hpa_unique.index:
            status = "human_ortholog_not_in_hpa"
        else:
            hpa_row = hpa_unique.loc[human_symbol]
            for tier in PRIMARY_TIERS:
                for location in split_locations(hpa_row[tier]):
                    primary[index, location_to_index[location]] = 1
            for tier in HIGH_CONFIDENCE_TIERS:
                for location in split_locations(hpa_row[tier]):
                    high_confidence[index, location_to_index[location]] = 1
            status = (
                "mapped_hpa_positive_profile"
                if primary[index].any()
                else "mapped_hpa_no_nonuncertain_location"
            )

        primary_locations = [
            locations[column]
            for column in np.flatnonzero(primary[index])
        ]
        high_locations = [
            locations[column]
            for column in np.flatnonzero(high_confidence[index])
        ]
        audit_rows.append(
            {
                "symbol": symbol,
                "display_symbol": clean(row.get("display_symbol", symbol)),
                "node_type": node_type,
                "mapping_status": status,
                "mgi_id": clean(orthology.get("mgi_id", "")),
                "mouse_entrez_gene_id": clean(
                    orthology.get("mouse_entrez_gene_id", "")
                ),
                "hgnc_id": clean(orthology.get("hgnc_id", "")),
                "human_symbol": human_symbol,
                "human_entrez_gene_id": clean(
                    orthology.get("human_entrez_gene_id", "")
                ),
                "hpa_ensembl_gene_id": (
                    clean(hpa_row["Gene"]) if hpa_row is not None else ""
                ),
                "hpa_gene_reliability": (
                    clean(hpa_row["Reliability"])
                    if hpa_row is not None
                    else ""
                ),
                "hpa_main_location": (
                    clean(hpa_row["Main location"])
                    if hpa_row is not None
                    else ""
                ),
                "hpa_additional_location": (
                    clean(hpa_row["Additional location"])
                    if hpa_row is not None
                    else ""
                ),
                "hpa_extracellular_location_audit_only": (
                    clean(hpa_row["Extracellular location"])
                    if hpa_row is not None
                    else ""
                ),
                "hpa_enhanced_locations": (
                    clean(hpa_row["Enhanced"]) if hpa_row is not None else ""
                ),
                "hpa_supported_locations": (
                    clean(hpa_row["Supported"]) if hpa_row is not None else ""
                ),
                "hpa_approved_locations": (
                    clean(hpa_row["Approved"]) if hpa_row is not None else ""
                ),
                "hpa_uncertain_locations_audit_only": (
                    clean(hpa_row["Uncertain"]) if hpa_row is not None else ""
                ),
                "hpa_profile_observed": bool(primary_locations),
                "hpa_location_count": len(primary_locations),
                "hpa_locations": ";".join(primary_locations),
                "hpa_high_confidence_profile_observed": bool(high_locations),
                "hpa_high_confidence_location_count": len(high_locations),
                "hpa_high_confidence_locations": ";".join(high_locations),
            }
        )

    audit = pd.DataFrame(audit_rows)
    for column_index, location in enumerate(locations):
        safe_name = re.sub(r"[^A-Za-z0-9]+", "_", location).strip("_").lower()
        audit[f"location__{safe_name}"] = primary[:, column_index]
    return audit, primary, high_confidence


def score_profiles(
    binary_profiles: np.ndarray,
) -> dict[str, object]:
    observed_mask = binary_profiles.sum(axis=1) > 0
    observed_indices = np.flatnonzero(observed_mask)
    observed_binary = binary_profiles[observed_mask].astype(float)
    norms = np.linalg.norm(observed_binary, axis=1, keepdims=True)
    unit_profiles = observed_binary / norms
    similarity_observed = unit_profiles @ unit_profiles.T
    threshold_background = similarity_observed.copy()
    np.fill_diagonal(threshold_background, np.nan)
    raw_thresholds = np.nanquantile(threshold_background, Q, axis=1)
    thresholds = raw_thresholds.copy()
    methods: list[str] = []
    positive_background_sizes: list[int] = []

    for index, raw_threshold in enumerate(raw_thresholds):
        positive = threshold_background[index][
            np.isfinite(threshold_background[index])
            & (threshold_background[index] > 0)
        ]
        positive_background_sizes.append(int(len(positive)))
        if np.isfinite(raw_threshold) and raw_threshold > 0:
            methods.append("all_observed_nodes_q75")
        elif len(positive):
            thresholds[index] = float(np.quantile(positive, Q))
            methods.append("positive_overlap_q75_fallback")
        else:
            thresholds[index] = np.nan
            methods.append("no_positive_overlap_neutral")

    directed = np.full_like(similarity_observed, MINIMUM_LIKELIHOOD)
    for index, threshold in enumerate(thresholds):
        if np.isfinite(threshold) and threshold > 0:
            z = similarity_observed[index] / threshold
            directed[index] = np.maximum(
                MINIMUM_LIKELIHOOD,
                1.0 - np.exp(-0.5 * np.square(z)),
            )
    np.fill_diagonal(directed, MINIMUM_LIKELIHOOD)
    symmetric = (directed + directed.T) / 2.0

    n = len(binary_profiles)
    similarity = np.zeros((n, n), dtype=float)
    likelihood = np.full((n, n), MINIMUM_LIKELIHOOD, dtype=float)
    factor_a = np.full((n, n), MINIMUM_LIKELIHOOD, dtype=float)
    similarity[np.ix_(observed_indices, observed_indices)] = (
        similarity_observed
    )
    likelihood[np.ix_(observed_indices, observed_indices)] = symmetric
    factor_a[np.ix_(observed_indices, observed_indices)] = directed
    np.fill_diagonal(similarity, 0.0)
    np.fill_diagonal(likelihood, 0.0)
    np.fill_diagonal(factor_a, 0.0)

    threshold_full = np.full(n, np.nan, dtype=float)
    raw_threshold_full = np.full(n, np.nan, dtype=float)
    method_full = np.full(n, "profile_not_observed", dtype=object)
    positive_size_full = np.zeros(n, dtype=int)
    threshold_full[observed_indices] = thresholds
    raw_threshold_full[observed_indices] = raw_thresholds
    method_full[observed_indices] = methods
    positive_size_full[observed_indices] = positive_background_sizes

    return {
        "observed_mask": observed_mask,
        "observed_indices": observed_indices,
        "similarity": similarity,
        "likelihood": likelihood,
        "directed_factor": factor_a,
        "threshold": threshold_full,
        "raw_threshold": raw_threshold_full,
        "threshold_method": method_full,
        "positive_background_size": positive_size_full,
    }


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    raw_dir = (
        project
        / "data/edge_characterization/localization/hpa/v25.1/raw"
    )
    hpa_path = raw_dir / "subcellular_location.tsv.zip"
    orthology_path = raw_dir / "HOM_ProteinCoding.rpt"
    current_matrix_path = (
        project
        / "results/edge_characterization/localization_kinase_predictor_string"
        / "combined_adjacency_matrix.tsv"
    )
    earlier_localization_path = (
        project
        / "results/edge_characterization/localization"
        / "localization_adjacency_matrix.tsv"
    )
    kinase_log_bf_path = (
        project
        / "results/edge_characterization/localization_kinase_predictor"
        / "kinase_predictor_log_bayes_factor_matrix.tsv"
    )
    string_observed_path = (
        project
        / "results/edge_characterization/localization_kinase_predictor_string"
        / "string_evidence_observed_matrix.tsv"
    )
    processed_dir = (
        project
        / "data/edge_characterization/localization/hpa/v25.1/processed"
    )
    output_dir = (
        project
        / "results/edge_characterization"
        / "localization_kinase_predictor_string_hpa"
    )
    processed_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in (
        universe_path,
        hpa_path,
        orthology_path,
        current_matrix_path,
        earlier_localization_path,
        kinase_log_bf_path,
        string_observed_path,
    ):
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    if universe["symbol"].duplicated().any():
        raise ValueError("Node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    symbol_array = np.asarray(symbols)
    if len(symbols) < 2:
        raise ValueError(f"Expected at least two nodes, found {len(symbols)}")

    hpa = load_hpa(hpa_path)
    orthology_raw, orthology_lookup, ambiguous_mouse = (
        load_one_to_one_orthology(orthology_path)
    )
    location_dictionary = extract_location_dictionary(hpa)
    primary_location_rows = location_dictionary.loc[
        location_dictionary["included_in_primary_profile"]
    ]
    locations = sorted(primary_location_rows["location"].tolist())
    if len(locations) != 49:
        raise ValueError(f"Expected 49 HPA locations, found {len(locations)}")

    profiles, primary_binary, high_binary = build_profiles(
        universe,
        hpa,
        orthology_lookup,
        ambiguous_mouse,
        locations,
    )
    primary_score = score_profiles(primary_binary)
    high_score = score_profiles(high_binary)

    for prefix, score in (
        ("hpa", primary_score),
        ("hpa_high_confidence", high_score),
    ):
        profiles[f"{prefix}_raw_Tq_q75"] = score["raw_threshold"]
        profiles[f"{prefix}_Tq_q75"] = score["threshold"]
        profiles[f"{prefix}_Tq_method"] = score["threshold_method"]
        profiles[f"{prefix}_positive_overlap_background_size"] = score[
            "positive_background_size"
        ]

    profiles_path = processed_dir / "node_hpa_localization_profiles.tsv"
    mapping_path = processed_dir / "mouse_human_hpa_mapping.tsv"
    location_dictionary_path = processed_dir / "hpa_location_dictionary.tsv"
    profiles.to_csv(profiles_path, sep="\t", index=False)
    mapping_columns = [
        column
        for column in profiles.columns
        if not column.startswith("location__")
    ]
    profiles[mapping_columns].to_csv(mapping_path, sep="\t", index=False)
    location_dictionary.to_csv(
        location_dictionary_path, sep="\t", index=False
    )

    current = pd.read_csv(
        current_matrix_path, sep="\t", index_col=0
    ).reindex(index=symbols, columns=symbols)
    if current.isna().any().any():
        raise ValueError("Current adjacency matrix does not align to universe")
    current_matrix = current.to_numpy(dtype=float)
    if not np.allclose(current_matrix, current_matrix.T, atol=1e-12):
        raise ValueError("Current adjacency matrix is not symmetric")
    if not np.allclose(np.diag(current_matrix), 0.0, atol=1e-12):
        raise ValueError("Current adjacency diagonal is not zero")

    def aligned_component(path: Path) -> np.ndarray:
        frame = pd.read_csv(path, sep="\t", index_col=0).reindex(
            index=symbols, columns=symbols
        )
        if frame.isna().any().any():
            raise ValueError(f"Evidence matrix does not align: {path}")
        return frame.to_numpy(dtype=float)

    earlier_localization = aligned_component(earlier_localization_path)
    kinase_log_bf = aligned_component(kinase_log_bf_path)
    string_observed = aligned_component(string_observed_path)

    n = len(symbols)
    upper_i, upper_j = np.triu_indices(n, k=1)
    edge_ids = pd.Index(
        [
            f"{symbols[i]}--{symbols[j]}"
            for i, j in zip(upper_i, upper_j, strict=True)
        ]
    )
    likelihood_matrix = np.asarray(primary_score["likelihood"], dtype=float)
    likelihood_upper = likelihood_matrix[upper_i, upper_j]
    current_upper = current_matrix[upper_i, upper_j]
    final_upper = binary_edge_bayes_update(
        pd.Series(current_upper, index=edge_ids),
        pd.Series(likelihood_upper, index=edge_ids),
        no_link_likelihood=NO_EDGE_LIKELIHOOD,
    ).to_numpy(dtype=float).copy()
    neutral_upper = likelihood_upper == MINIMUM_LIKELIHOOD
    final_upper[neutral_upper] = current_upper[neutral_upper]

    final_matrix = np.zeros((n, n), dtype=float)
    final_matrix[upper_i, upper_j] = final_upper
    final_matrix[upper_j, upper_i] = final_upper
    log_bf_matrix = np.zeros((n, n), dtype=float)
    log_bf_upper = np.log(likelihood_upper / NO_EDGE_LIKELIHOOD)
    log_bf_matrix[upper_i, upper_j] = log_bf_upper
    log_bf_matrix[upper_j, upper_i] = log_bf_upper
    observed_pair_matrix = np.outer(
        np.asarray(primary_score["observed_mask"], dtype=np.uint8),
        np.asarray(primary_score["observed_mask"], dtype=np.uint8),
    )
    np.fill_diagonal(observed_pair_matrix, 0)

    write_matrix(
        output_dir / "hpa_localization_similarity_matrix.tsv",
        symbols,
        np.asarray(primary_score["similarity"], dtype=float),
    )
    write_matrix(
        output_dir / "hpa_localization_likelihood_matrix.tsv",
        symbols,
        likelihood_matrix,
    )
    write_matrix(
        output_dir / "hpa_localization_log_bayes_factor_matrix.tsv",
        symbols,
        log_bf_matrix,
    )
    write_matrix(
        output_dir / "hpa_localization_observed_pair_matrix.tsv",
        symbols,
        observed_pair_matrix,
    )
    write_matrix(
        output_dir / "hpa_high_confidence_likelihood_matrix.tsv",
        symbols,
        np.asarray(high_score["likelihood"], dtype=float),
    )
    final_matrix_path = output_dir / "combined_adjacency_matrix.tsv"
    write_matrix(final_matrix_path, symbols, final_matrix)

    human_symbol = profiles["human_symbol"].to_numpy(dtype=object)
    reliability = profiles["hpa_gene_reliability"].to_numpy(dtype=object)
    primary_observed = np.asarray(
        primary_score["observed_mask"], dtype=bool
    )
    pair_observed = primary_observed[upper_i] & primary_observed[upper_j]
    similarity_upper = np.asarray(primary_score["similarity"])[
        upper_i, upper_j
    ]
    directed = np.asarray(primary_score["directed_factor"])
    factor_using_a = directed[upper_i, upper_j]
    factor_using_b = directed[upper_j, upper_i]
    threshold = np.asarray(primary_score["threshold"])
    threshold_method = np.asarray(
        primary_score["threshold_method"], dtype=object
    )
    strengthened = likelihood_upper > MINIMUM_LIKELIHOOD + 1e-12
    earlier_localization_non_neutral = (
        earlier_localization[upper_i, upper_j] > 0.5 + 1e-12
    )
    kinase_non_neutral = kinase_log_bf[upper_i, upper_j] > 1e-12
    string_reported = string_observed[upper_i, upper_j] > 0
    no_prior_stream_non_neutral = ~(
        earlier_localization_non_neutral
        | kinase_non_neutral
        | string_reported
    )

    all_edges = pd.DataFrame(
        {
            "node_a": symbol_array[upper_i],
            "node_b": symbol_array[upper_j],
            "human_ortholog_a": human_symbol[upper_i],
            "human_ortholog_b": human_symbol[upper_j],
            "hpa_gene_reliability_a": reliability[upper_i],
            "hpa_gene_reliability_b": reliability[upper_j],
            "hpa_profiles_observed_for_both_nodes": pair_observed,
            "hpa_cosine_similarity": np.where(
                pair_observed, similarity_upper, np.nan
            ),
            "hpa_Tq_node_a": threshold[upper_i],
            "hpa_Tq_method_node_a": threshold_method[upper_i],
            "hpa_Tq_node_b": threshold[upper_j],
            "hpa_Tq_method_node_b": threshold_method[upper_j],
            "hpa_factor_using_Tq_node_a": factor_using_a,
            "hpa_factor_using_Tq_node_b": factor_using_b,
            "hpa_undirected_likelihood": likelihood_upper,
            "hpa_log_bayes_factor": log_bf_upper,
            "hpa_bayes_factor": likelihood_upper / NO_EDGE_LIKELIHOOD,
            "pre_hpa_posterior": current_upper,
            "final_posterior": final_upper,
            "posterior_increase": final_upper - current_upper,
            "strengthened_by_hpa": strengthened,
        }
    )
    all_edges_path = output_dir / "combined_undirected_edges.tsv.gz"
    with gzip.open(
        all_edges_path, "wt", encoding="utf-8", newline=""
    ) as handle:
        all_edges.to_csv(
            handle, sep="\t", index=False, float_format="%.9g"
        )
    strengthened_edges = all_edges.loc[strengthened].sort_values(
        [
            "hpa_undirected_likelihood",
            "final_posterior",
            "node_a",
            "node_b",
        ],
        ascending=[False, False, True, True],
        kind="stable",
    )
    strengthened_path = (
        output_dir / "hpa_strengthened_undirected_edges.tsv"
    )
    strengthened_edges.to_csv(
        strengthened_path, sep="\t", index=False, float_format="%.9g"
    )
    all_edges.sort_values(
        [
            "final_posterior",
            "hpa_undirected_likelihood",
            "node_a",
            "node_b",
        ],
        ascending=[False, False, True, True],
        kind="stable",
    ).head(1000).to_csv(
        output_dir / "top_combined_undirected_edges.tsv",
        sep="\t",
        index=False,
        float_format="%.9g",
    )

    high_likelihood_upper = np.asarray(high_score["likelihood"])[
        upper_i, upper_j
    ]
    high_observed = np.asarray(high_score["observed_mask"], dtype=bool)
    high_pair_observed = high_observed[upper_i] & high_observed[upper_j]
    high_strengthened = (
        high_likelihood_upper > MINIMUM_LIKELIHOOD + 1e-12
    )
    compare_mask = high_pair_observed
    correlation = (
        float(
            np.corrcoef(
                likelihood_upper[compare_mask],
                high_likelihood_upper[compare_mask],
            )[0, 1]
        )
        if compare_mask.sum() > 1
        else float("nan")
    )
    union_strengthened = strengthened | high_strengthened
    jaccard = (
        float((strengthened & high_strengthened).sum())
        / float(union_strengthened.sum())
        if union_strengthened.any()
        else 1.0
    )

    molecule_mask = universe["node_type"].eq("molecule").to_numpy()
    molecule_pair = molecule_mask[upper_i] | molecule_mask[upper_j]
    validation = {
        "matrix_shape": list(final_matrix.shape),
        "matrix_symmetric": bool(
            np.array_equal(final_matrix, final_matrix.T)
        ),
        "diagonal_all_zero": bool(
            np.array_equal(np.diag(final_matrix), np.zeros(n))
        ),
        "probabilities_bounded_zero_one": bool(
            np.all((final_matrix >= 0.0) & (final_matrix <= 1.0))
        ),
        "neutral_pairs_exactly_unchanged": bool(
            np.array_equal(
                final_upper[neutral_upper],
                current_upper[neutral_upper],
            )
        ),
        "molecule_incident_pairs_exactly_unchanged": bool(
            np.array_equal(
                final_upper[molecule_pair],
                current_upper[molecule_pair],
            )
        ),
        "only_protein_pairs_strengthened": bool(
            not np.any(strengthened & molecule_pair)
        ),
        "strengthened_edge_table_has_no_duplicate_pairs": bool(
            not strengthened_edges.duplicated(["node_a", "node_b"]).any()
        ),
        "node_order_preserved": bool(
            current.index.astype(str).tolist() == symbols
            and current.columns.astype(str).tolist() == symbols
        ),
    }
    if not all(validation.values()):
        raise AssertionError(f"Validation failed: {validation}")

    profile_status_counts = {
        str(key): int(value)
        for key, value in profiles["mapping_status"].value_counts().items()
    }
    tq_method_counts = {
        str(key): int(value)
        for key, value in pd.Series(
            primary_score["threshold_method"]
        ).value_counts().items()
    }
    high_tq_method_counts = {
        str(key): int(value)
        for key, value in pd.Series(
            high_score["threshold_method"]
        ).value_counts().items()
    }
    final_minimum = float(final_upper.min())
    summary = {
        "hpa_version": HPA_VERSION,
        "hpa_ensembl_version": HPA_ENSEMBL_VERSION,
        "node_count": n,
        "protein_node_count": int(
            universe["node_type"].eq("protein").sum()
        ),
        "molecule_node_count": int(molecule_mask.sum()),
        "hpa_download_gene_count": int(len(hpa)),
        "hpa_location_count": len(locations),
        "mgi_report_row_count": int(len(orthology_raw)),
        "protein_nodes_with_mgi_one_to_one_human_ortholog": int(
            profiles["human_symbol"].ne("").sum()
        ),
        "protein_nodes_with_unique_hpa_row": int(
            profiles["hpa_ensembl_gene_id"].ne("").sum()
        ),
        "protein_nodes_with_primary_hpa_profile": int(
            primary_observed.sum()
        ),
        "protein_nodes_with_high_confidence_hpa_profile": int(
            high_observed.sum()
        ),
        "mapping_status_counts": profile_status_counts,
        "primary_tiers": list(PRIMARY_TIERS),
        "excluded_from_positive_evidence": [
            "Uncertain",
            "Extracellular location",
        ],
        "profile_encoding": (
            "Binary membership in each of 49 HPA locations, then L2 "
            "unit-normalization; dot product therefore equals cosine "
            "similarity."
        ),
        "possible_undirected_pairs": int(len(all_edges)),
        "pairs_with_two_primary_hpa_profiles": int(pair_observed.sum()),
        "pairs_with_positive_hpa_location_overlap": int(
            (pair_observed & (similarity_upper > 0)).sum()
        ),
        "pairs_strengthened_by_hpa": int(strengthened.sum()),
        "hpa_strengthened_pairs_previously_at_minimum": int(
            (strengthened & (current_upper == 0.5)).sum()
        ),
        "hpa_strengthened_pairs_previously_above_minimum": int(
            (strengthened & (current_upper > 0.5)).sum()
        ),
        "nonminimum_edge_count_before_hpa": int(
            (current_upper > 0.5).sum()
        ),
        "nonminimum_edge_count_after_hpa": int(
            (final_upper > 0.5).sum()
        ),
        "pairs_neutral_under_hpa": int(neutral_upper.sum()),
        "pairs_at_final_minimum_probability": int(
            (final_upper == final_minimum).sum()
        ),
        "final_minimum_off_diagonal_probability": final_minimum,
        "final_maximum_off_diagonal_probability": float(final_upper.max()),
        "hpa_likelihood_minimum": float(likelihood_upper.min()),
        "hpa_likelihood_maximum": float(likelihood_upper.max()),
        "hpa_bayes_factor_minimum": float(
            (likelihood_upper / NO_EDGE_LIKELIHOOD).min()
        ),
        "hpa_bayes_factor_maximum": float(
            (likelihood_upper / NO_EDGE_LIKELIHOOD).max()
        ),
        "q": Q,
        "minimum_likelihood": MINIMUM_LIKELIHOOD,
        "no_edge_likelihood": NO_EDGE_LIKELIHOOD,
        "Tq_method_counts": tq_method_counts,
        "background_definition": (
            "For each HPA-profiled node, cosine similarities to every other "
            f"HPA-profiled node in the {n}-node universe, self excluded. If "
            "the all-node q75 is zero, q75 of positive overlaps is used and "
            "flagged; a node with no positive overlap is neutral."
        ),
        "directed_likelihood_formula": (
            "max(0.5, 1-exp(-0.5*(cosine_similarity/Tq_node)^2))"
        ),
        "symmetrization": (
            "Arithmetic mean of the factors calculated with Tq_node_a and "
            "Tq_node_b."
        ),
        "hpa_bayes_factor_formula": (
            "BF_HPA=hpa_undirected_likelihood/0.5"
        ),
        "integration_formula": (
            "Independent Bernoulli update of the current localization + "
            "KinasePredictor + STRING posterior with the HPA likelihood; "
            "equivalently logit(P_final)=logit(P_current)+log(BF_HPA)."
        ),
        "missing_evidence_rule": (
            "No ortholog, no unique HPA row, no non-uncertain location, no "
            "shared location, or a score at the 0.5 floor is neutral and "
            "leaves the current edge posterior exactly unchanged."
        ),
        "cross_species_rule": (
            "Only MGI stringent one-to-one protein-coding mouse-human "
            "orthologs; no symbol-only or paralog fallback."
        ),
        "high_confidence_sensitivity": {
            "included_tiers": list(HIGH_CONFIDENCE_TIERS),
            "profiled_node_count": int(high_observed.sum()),
            "pairs_with_two_profiles": int(high_pair_observed.sum()),
            "pairs_strengthened": int(high_strengthened.sum()),
            "Tq_method_counts": high_tq_method_counts,
            "likelihood_correlation_on_high_confidence_observed_pairs": (
                correlation
            ),
            "strengthened_edge_set_jaccard": jaccard,
        },
        "overlap_of_hpa_strengthened_edges_with_prior_streams": {
            "earlier_localization": int(
                (strengthened & earlier_localization_non_neutral).sum()
            ),
            "kinase_predictor": int(
                (strengthened & kinase_non_neutral).sum()
            ),
            "string_reported": int(
                (strengthened & string_reported).sum()
            ),
            "hpa_only_among_these_streams": int(
                (strengthened & no_prior_stream_non_neutral).sum()
            ),
        },
        "dependence_caveat": (
            "HPA reliability uses external localization evidence including "
            "UniProt comparisons, and human cell-line localization may not "
            "match mouse renal principal cells. HPA is therefore a separate "
            "supporting stream, not a replacement for mpkCCD fractionation."
        ),
        "validation": validation,
        "source_urls": {
            "hpa_subcellular_page": HPA_PAGE_URL,
            "hpa_methods": HPA_METHOD_URL,
            "hpa_data": HPA_DATA_URL,
            "mgi_report_page": MGI_REPORT_PAGE,
            "mgi_one_to_one_report": MGI_ORTHOLOGY_URL,
        },
        "input_sha256": {
            path.name: sha256_file(path)
            for path in (
                universe_path,
                hpa_path,
                orthology_path,
                current_matrix_path,
                earlier_localization_path,
                kinase_log_bf_path,
                string_observed_path,
            )
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    readme = f"""HPA v{HPA_VERSION} LOCALIZATION EDGE INTEGRATION

PURPOSE
This directory extends the UNDIRECTED {n}-node adjacency model with Human
Protein Atlas subcellular localization. HPA is added after localization,
KinasePredictor, and STRING; it does not replace the mouse mpkCCD
fractionation stream.

INPUT AND MAPPING
- HPA v{HPA_VERSION}, Ensembl {HPA_ENSEMBL_VERSION}: {len(hpa):,} human genes
  and {len(locations)} subcellular locations.
- Mouse proteins are mapped only through MGI's stringent one-to-one
  protein-coding orthology report.
- No case-converted-symbol or paralog fallback is used.
- Small-molecule second messengers receive no HPA evidence.

RESULT
- {int(primary_observed.sum()):,} protein nodes have a primary HPA profile.
- {int(strengthened.sum()):,} edges receive non-neutral HPA evidence.
- {int((strengthened & (current_upper == 0.5)).sum()):,} of those edges were
  previously at 0.5.
- The number of non-minimum edges changes from
  {int((current_upper > 0.5).sum()):,} to {int((final_upper > 0.5).sum()):,}.

PROFILE DEFINITION
Enhanced, Supported, and Approved HPA locations are included as binary
annotations. Uncertain annotations are retained in the mapping audit but
excluded from positive evidence. The separate Extracellular location field is
also audit-only because it is not one of the tiered ICC-IF assignments used by
this model.

Each 49-location binary vector is divided by its L2 norm. A pair's dot product
is therefore cosine similarity, preventing broadly annotated proteins from
receiving larger scores solely because they have more location labels.

BACKGROUND AND LIKELIHOOD
For each profiled node, T_q is the 75th percentile of its cosine similarities
to every other HPA-profiled node in the current universe, excluding itself. If
that percentile is zero, T_q is the 75th percentile of positive overlaps and
the fallback is flagged. A node with no positive overlap remains neutral.

For each endpoint-specific threshold:

    L_i = max(0.5, 1 - exp(-0.5 * (similarity / T_q_i)^2))

The undirected HPA likelihood is (L_a + L_b)/2. Its Bayes factor is:

    BF_HPA = L_undirected / 0.5

The current edge posterior is updated independently:

    logit(P_final) = logit(P_current) + log(BF_HPA)

Missing coverage, zero overlap, and likelihoods at the 0.5 floor are neutral;
the prior edge value is preserved exactly.

SENSITIVITY ANALYSIS
hpa_high_confidence_likelihood_matrix.tsv repeats the analysis using only
Enhanced and Supported annotations. It is provided for robustness checking and
is not additionally integrated, which would double count HPA.

FILES
- combined_adjacency_matrix.tsv: final symmetric {n} x {n} probabilities.
- combined_undirected_edges.tsv.gz: all {len(all_edges):,} unordered pairs and audit data.
- hpa_strengthened_undirected_edges.tsv: pairs with BF_HPA > 1.
- hpa_localization_similarity_matrix.tsv: primary cosine similarities.
- hpa_localization_likelihood_matrix.tsv: primary undirected likelihoods.
- hpa_localization_log_bayes_factor_matrix.tsv: additive HPA log evidence.
- hpa_localization_observed_pair_matrix.tsv: coverage indicator.
- hpa_high_confidence_likelihood_matrix.tsv: Enhanced+Supported sensitivity.
- top_combined_undirected_edges.tsv: top 1,000 final posterior edges.
- analysis_summary.json: counts, formulas, hashes, sensitivity, validation.

PROCESSED AUDITS
{mapping_path}
{profiles_path}
{location_dictionary_path}

LIMITATIONS
HPA profiles human proteins in up to three cultured cell lines and uses
antibody-based ICC-IF/confocal imaging. Localization can be species-, tissue-,
cell-state-, and isoform-dependent. HPA reliability also considers external
evidence including UniProt, so strict conditional independence from every
future localization source should not be assumed.

OFFICIAL SOURCES
{HPA_PAGE_URL}
{HPA_METHOD_URL}
{HPA_DATA_URL}
{MGI_REPORT_PAGE}
{MGI_ORTHOLOGY_URL}
"""
    readme_path = output_dir / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")
    summary["output_sha256"] = {
        path.name: sha256_file(path)
        for path in (
            mapping_path,
            profiles_path,
            location_dictionary_path,
            final_matrix_path,
            strengthened_path,
            all_edges_path,
            readme_path,
        )
    }
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
