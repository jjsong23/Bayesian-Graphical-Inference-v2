#!/usr/bin/env python3
"""Integrate OmniPath core signaling interactions as undirected edge evidence.

The official mouse OmniPath core post-translational network is directed and
retains stimulation/inhibition annotations. This workflow preserves those
fields for audit, but collapses all records onto one unordered node pair for
the modeled graph. Small-molecule nodes and absent OmniPath pairs are neutral.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

CODE_DIR = Path(__file__).resolve().parents[1]
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from bayes_factors import binary_edge_bayes_update


DOWNLOAD_DATE = "2026-07-30"
DATASET = "omnipath"
ORGANISM = 10090
INTERACTION_TYPE = "post_translational"
ENTITY_TYPE = "protein"
Q = 0.75
NO_EDGE_LIKELIHOOD = 0.5
API_QUERY = (
    "https://omnipathdb.org/interactions?"
    "genesymbols=yes&datasets=omnipath&organisms=10090"
    "&types=post_translational&entity_types=protein&directed=yes"
    "&fields=sources,references,curation_effort,entity_type,type"
    "&format=tsv&license=academic"
)
API_QUERY_DOCS = "https://omnipathdb.org/queries/interactions"
OMNIPATH_HOME = "https://omnipathdb.org/"
OMNIPATH_R_DOCS = (
    "https://r.omnipathdb.org/reference/omnipath-interactions.html"
)


def parse_args() -> argparse.Namespace:
    default_project = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=default_project,
        help=f"Project directory (default: {default_project})",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_matrix(path: Path, symbols: list[str], matrix: np.ndarray) -> None:
    frame = pd.DataFrame(matrix, index=symbols, columns=symbols)
    frame.index.name = "symbol"
    frame.to_csv(path, sep="\t", float_format="%.9g")


def parse_bool(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes"}


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "-"} else text


def split_semicolon(value: object) -> set[str]:
    return {
        item.strip()
        for item in clean(value).split(";")
        if item.strip()
    }


def canonical_pair(
    node_a: str,
    node_b: str,
    symbol_to_index: dict[str, int],
) -> tuple[str, str]:
    if symbol_to_index[node_a] < symbol_to_index[node_b]:
        return node_a, node_b
    return node_b, node_a


def endpoint_mapping(
    api_symbol: object,
    api_uniprot: object,
    protein_symbols: set[str],
    accession_to_symbol: dict[str, str],
) -> tuple[str, str]:
    symbol = clean(api_symbol)
    accession = clean(api_uniprot)
    if symbol in protein_symbols:
        return symbol, "exact_mouse_gene_symbol"
    mapped = accession_to_symbol.get(accession, "")
    if mapped:
        return mapped, "exact_mouse_uniprot_fallback"
    return "", "unmapped"


def full_background(raw: pd.DataFrame) -> pd.DataFrame:
    """Collapse the full mouse query to the same unordered scoring unit."""

    background = raw.loc[
        raw["source_genesymbol"].ne("")
        & raw["target_genesymbol"].ne("")
        & raw["source_genesymbol"].ne(raw["target_genesymbol"])
    ].copy()
    background["node_a"] = background[
        ["source_genesymbol", "target_genesymbol"]
    ].min(axis=1)
    background["node_b"] = background[
        ["source_genesymbol", "target_genesymbol"]
    ].max(axis=1)
    collapsed = (
        background.groupby(["node_a", "node_b"], as_index=False, sort=False)
        .agg(
            omnipath_curation_effort=(
                "curation_effort",
                "max",
            ),
            directed_record_count=("curation_effort", "size"),
        )
    )
    return collapsed


def collapse_universe_pairs(
    mapped: pd.DataFrame,
    symbol_to_index: dict[str, int],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (node_a, node_b), group in mapped.groupby(
        ["node_a", "node_b"],
        sort=False,
    ):
        directions = sorted(
            {
                f"{source}->{target}"
                for source, target in zip(
                    group["mapped_source"],
                    group["mapped_target"],
                    strict=True,
                )
            }
        )
        resources: set[str] = set()
        references: set[str] = set()
        source_accessions: set[str] = set()
        target_accessions: set[str] = set()
        for record in group.to_dict(orient="records"):
            resources.update(split_semicolon(record["sources"]))
            references.update(split_semicolon(record["references"]))
            source_accessions.add(clean(record["source"]))
            target_accessions.add(clean(record["target"]))
        resources.discard("")
        references.discard("")
        source_accessions.discard("")
        target_accessions.discard("")

        any_stimulation = bool(group["is_stimulation_bool"].any())
        any_inhibition = bool(group["is_inhibition_bool"].any())
        any_unsigned = bool(
            (
                ~group["is_stimulation_bool"]
                & ~group["is_inhibition_bool"]
            ).any()
        )
        rows.append(
            {
                "node_a": node_a,
                "node_b": node_b,
                "directed_record_count": int(len(group)),
                "directional_interactions": ";".join(directions),
                "omnipath_bidirectional": len(directions) > 1,
                "any_stimulation_annotation": any_stimulation,
                "any_inhibition_annotation": any_inhibition,
                "any_unsigned_annotation": any_unsigned,
                "stimulation_and_inhibition_both_present": (
                    any_stimulation and any_inhibition
                ),
                "any_consensus_direction": bool(
                    group["consensus_direction_bool"].any()
                ),
                "any_consensus_stimulation": bool(
                    group["consensus_stimulation_bool"].any()
                ),
                "any_consensus_inhibition": bool(
                    group["consensus_inhibition_bool"].any()
                ),
                "omnipath_curation_effort": float(
                    group["curation_effort"].max()
                ),
                "curation_effort_sum_across_directed_records": float(
                    group["curation_effort"].sum()
                ),
                "curation_effort_values": ";".join(
                    map(
                        str,
                        sorted(
                            {
                                int(value)
                                for value in group["curation_effort"]
                            }
                        ),
                    )
                ),
                "unique_resource_count": len(resources),
                "resources": ";".join(sorted(resources)),
                "unique_reference_count": len(references),
                "references": ";".join(sorted(references)),
                "source_uniprot_accessions": ";".join(
                    sorted(source_accessions)
                ),
                "target_uniprot_accessions": ";".join(
                    sorted(target_accessions)
                ),
                "endpoint_mapping_methods": ";".join(
                    sorted(
                        set(group["source_mapping_method"])
                        | set(group["target_mapping_method"])
                    )
                ),
            }
        )
    collapsed = pd.DataFrame(rows)
    collapsed["_order_a"] = collapsed["node_a"].map(symbol_to_index)
    collapsed["_order_b"] = collapsed["node_b"].map(symbol_to_index)
    collapsed = (
        collapsed.sort_values(
            ["_order_a", "_order_b"],
            kind="stable",
        )
        .drop(columns=["_order_a", "_order_b"])
        .reset_index(drop=True)
    )
    return collapsed


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()

    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    protein_index_path = (
        project
        / "data/edge_characterization/kinase_predictor"
        / "phosphosite_database/protein_index.tsv"
    )
    raw_path = (
        project
        / "data/edge_characterization/omnipath"
        / DOWNLOAD_DATE
        / "raw/omnipath_mouse_core_post_translational.tsv"
    )
    processed_dir = (
        project
        / "data/edge_characterization/omnipath"
        / DOWNLOAD_DATE
        / "processed"
    )
    current_dir = (
        project
        / "results/edge_characterization"
        / "localization_kinase_predictor_string_hpa"
    )
    current_matrix_path = current_dir / "combined_adjacency_matrix.tsv"
    current_edges_path = current_dir / "combined_undirected_edges.tsv.gz"
    output_dir = (
        project
        / "results/edge_characterization"
        / "localization_kinase_predictor_string_hpa_omnipath"
    )
    backend_dir = (
        project
        / "results/backend_bayes_factor_catalogs/edge_factors_891"
    )
    localization_matrix_path = (
        project
        / "results/edge_characterization/localization"
        / "localization_adjacency_matrix.tsv"
    )
    kinase_matrix_path = (
        project
        / "results/edge_characterization/localization_kinase_predictor"
        / "kinase_predictor_log_bayes_factor_matrix.tsv"
    )
    string_matrix_path = (
        project
        / "results/edge_characterization/localization_kinase_predictor_string"
        / "string_evidence_observed_matrix.tsv"
    )
    hpa_matrix_path = (
        current_dir / "hpa_localization_log_bayes_factor_matrix.tsv"
    )
    processed_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    backend_dir.mkdir(parents=True, exist_ok=True)

    required = [
        universe_path,
        protein_index_path,
        raw_path,
        current_matrix_path,
        current_edges_path,
        localization_matrix_path,
        kinase_matrix_path,
        string_matrix_path,
        hpa_matrix_path,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    universe = pd.read_csv(
        universe_path,
        sep="\t",
        dtype=str,
    ).fillna("")
    if universe["symbol"].duplicated().any():
        raise ValueError("Node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    symbol_to_index = {
        symbol: index for index, symbol in enumerate(symbols)
    }
    protein_symbols = set(
        universe.loc[universe["node_type"].eq("protein"), "symbol"]
    )
    molecule_symbols = set(
        universe.loc[universe["node_type"].eq("molecule"), "symbol"]
    )
    n = len(symbols)
    if n != 891:
        raise ValueError(f"Expected 891 nodes, found {n}")

    protein_index = pd.read_csv(
        protein_index_path,
        sep="\t",
        dtype=str,
    ).fillna("")
    if set(protein_index["symbol"]) != protein_symbols:
        raise ValueError("Protein index does not match protein node universe")
    accession_to_symbol: dict[str, str] = {}
    for symbol, accession in zip(
        protein_index["symbol"],
        protein_index["mapped_uniprot"],
        strict=True,
    ):
        accession = clean(accession)
        if not accession:
            continue
        previous = accession_to_symbol.get(accession)
        if previous is not None and previous != symbol:
            raise ValueError(
                f"UniProt accession maps to two nodes: "
                f"{accession} -> {previous}, {symbol}"
            )
        accession_to_symbol[accession] = str(symbol)

    raw = pd.read_csv(raw_path, sep="\t", dtype=str).fillna("")
    required_columns = {
        "source",
        "target",
        "source_genesymbol",
        "target_genesymbol",
        "is_directed",
        "is_stimulation",
        "is_inhibition",
        "consensus_direction",
        "consensus_stimulation",
        "consensus_inhibition",
        "sources",
        "references",
        "type",
        "curation_effort",
        "entity_type_source",
        "entity_type_target",
    }
    missing_columns = required_columns.difference(raw.columns)
    if missing_columns:
        raise ValueError(
            f"Raw OmniPath table is missing columns: "
            f"{sorted(missing_columns)}"
        )
    raw["curation_effort"] = pd.to_numeric(
        raw["curation_effort"],
        errors="raise",
    ).astype(float)
    if (
        (~np.isfinite(raw["curation_effort"]))
        | (raw["curation_effort"] <= 0)
    ).any():
        raise ValueError("OmniPath curation_effort must be positive finite")
    for column in (
        "is_directed",
        "is_stimulation",
        "is_inhibition",
        "consensus_direction",
        "consensus_stimulation",
        "consensus_inhibition",
    ):
        raw[f"{column}_bool"] = raw[column].map(parse_bool)
    if not raw["is_directed_bool"].all():
        raise ValueError("Downloaded OmniPath query contains undirected rows")
    if set(raw["type"]) != {INTERACTION_TYPE}:
        raise ValueError(
            f"Unexpected OmniPath interaction types: {set(raw['type'])}"
        )
    if set(raw["entity_type_source"]) != {ENTITY_TYPE}:
        raise ValueError("Unexpected OmniPath source entity type")
    if set(raw["entity_type_target"]) != {ENTITY_TYPE}:
        raise ValueError("Unexpected OmniPath target entity type")

    background = full_background(raw)
    t_q = float(background["omnipath_curation_effort"].quantile(Q))
    if not np.isfinite(t_q) or t_q <= 0:
        raise ValueError(f"Invalid OmniPath T_q: {t_q}")
    background_path = (
        processed_dir / "omnipath_background_undirected_pairs.tsv.gz"
    )
    with gzip.open(
        background_path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        background.to_csv(
            handle,
            sep="\t",
            index=False,
            float_format="%.9g",
        )

    source_mappings = [
        endpoint_mapping(
            symbol,
            accession,
            protein_symbols,
            accession_to_symbol,
        )
        for symbol, accession in zip(
            raw["source_genesymbol"],
            raw["source"],
            strict=True,
        )
    ]
    target_mappings = [
        endpoint_mapping(
            symbol,
            accession,
            protein_symbols,
            accession_to_symbol,
        )
        for symbol, accession in zip(
            raw["target_genesymbol"],
            raw["target"],
            strict=True,
        )
    ]
    raw["mapped_source"] = [item[0] for item in source_mappings]
    raw["source_mapping_method"] = [item[1] for item in source_mappings]
    raw["mapped_target"] = [item[0] for item in target_mappings]
    raw["target_mapping_method"] = [item[1] for item in target_mappings]
    in_universe = raw.loc[
        raw["mapped_source"].ne("")
        & raw["mapped_target"].ne("")
        & raw["mapped_source"].ne(raw["mapped_target"])
    ].copy()
    canonical = [
        canonical_pair(source, target, symbol_to_index)
        for source, target in zip(
            in_universe["mapped_source"],
            in_universe["mapped_target"],
            strict=True,
        )
    ]
    in_universe["node_a"] = [item[0] for item in canonical]
    in_universe["node_b"] = [item[1] for item in canonical]
    supported = collapse_universe_pairs(
        in_universe,
        symbol_to_index,
    )
    if supported.empty:
        raise ValueError("No in-universe OmniPath interactions were retained")

    # Sparse curated presence is explicitly positive. This zero-anchored
    # complement transform assigns missing pairs L=0.5 exactly while every
    # observed positive curation effort receives L>0.5 and BF>1.
    effort = supported["omnipath_curation_effort"].to_numpy(dtype=float)
    support_component = 1.0 - np.exp(
        -0.5 * np.square(effort / t_q)
    )
    likelihood = (
        NO_EDGE_LIKELIHOOD
        + (1.0 - NO_EDGE_LIKELIHOOD) * support_component
    )
    bayes_factor = likelihood / NO_EDGE_LIKELIHOOD
    log_bayes_factor = np.log(bayes_factor)
    supported["omnipath_Tq"] = t_q
    supported["omnipath_likelihood"] = likelihood
    supported["omnipath_bayes_factor"] = bayes_factor
    supported["omnipath_log_bayes_factor"] = log_bayes_factor

    current = pd.read_csv(
        current_matrix_path,
        sep="\t",
        index_col=0,
    ).reindex(index=symbols, columns=symbols)
    if current.isna().any().any():
        raise ValueError("Current HPA matrix does not align to universe")
    current_matrix = current.to_numpy(dtype=float)
    if not np.allclose(current_matrix, current_matrix.T, atol=1e-12):
        raise ValueError("Current HPA matrix is not symmetric")
    if not np.allclose(np.diag(current_matrix), 0.0, atol=1e-12):
        raise ValueError("Current HPA matrix diagonal is not zero")

    observed_matrix = np.zeros((n, n), dtype=np.uint8)
    effort_matrix = np.zeros((n, n), dtype=float)
    likelihood_matrix = np.full((n, n), NO_EDGE_LIKELIHOOD)
    log_bf_matrix = np.zeros((n, n), dtype=float)
    np.fill_diagonal(likelihood_matrix, 0.0)

    pair_ids: list[str] = []
    pair_current: list[float] = []
    for row in supported.to_dict(orient="records"):
        i = symbol_to_index[str(row["node_a"])]
        j = symbol_to_index[str(row["node_b"])]
        observed_matrix[i, j] = observed_matrix[j, i] = 1
        effort_matrix[i, j] = effort_matrix[j, i] = float(
            row["omnipath_curation_effort"]
        )
        likelihood_matrix[i, j] = likelihood_matrix[j, i] = float(
            row["omnipath_likelihood"]
        )
        log_bf_matrix[i, j] = log_bf_matrix[j, i] = float(
            row["omnipath_log_bayes_factor"]
        )
        pair_ids.append(f"{row['node_a']}--{row['node_b']}")
        pair_current.append(float(current_matrix[i, j]))

    updated = binary_edge_bayes_update(
        pd.Series(pair_current, index=pair_ids, dtype=float),
        pd.Series(likelihood, index=pair_ids, dtype=float),
        no_link_likelihood=NO_EDGE_LIKELIHOOD,
    )
    final_matrix = current_matrix.copy()
    final_values = updated.to_numpy(dtype=float)
    for row, final_probability in zip(
        supported.to_dict(orient="records"),
        final_values,
        strict=True,
    ):
        i = symbol_to_index[str(row["node_a"])]
        j = symbol_to_index[str(row["node_b"])]
        final_matrix[i, j] = final_matrix[j, i] = final_probability
    np.fill_diagonal(final_matrix, 0.0)
    supported["pre_omnipath_posterior"] = pair_current
    supported["final_posterior"] = final_values
    supported["posterior_increase"] = (
        supported["final_posterior"]
        - supported["pre_omnipath_posterior"]
    )
    supported["newly_nonbaseline_edge"] = np.isclose(
        supported["pre_omnipath_posterior"],
        0.5,
        atol=1e-12,
        rtol=0.0,
    )

    supported_path = (
        output_dir / "omnipath_supported_undirected_edges.tsv"
    )
    supported.to_csv(
        supported_path,
        sep="\t",
        index=False,
        float_format="%.9g",
    )
    write_matrix(
        output_dir / "omnipath_evidence_observed_matrix.tsv",
        symbols,
        observed_matrix,
    )
    write_matrix(
        output_dir / "omnipath_curation_effort_matrix.tsv",
        symbols,
        effort_matrix,
    )
    write_matrix(
        output_dir / "omnipath_likelihood_matrix.tsv",
        symbols,
        likelihood_matrix,
    )
    write_matrix(
        output_dir / "omnipath_log_bayes_factor_matrix.tsv",
        symbols,
        log_bf_matrix,
    )
    final_matrix_path = output_dir / "combined_adjacency_matrix.tsv"
    write_matrix(final_matrix_path, symbols, final_matrix)

    upper_i, upper_j = np.triu_indices(n, k=1)
    current_upper = current_matrix[upper_i, upper_j]
    final_upper = final_matrix[upper_i, upper_j]
    observed_upper = observed_matrix[upper_i, upper_j].astype(bool)
    molecule_mask = universe["node_type"].eq("molecule").to_numpy()
    molecule_pair = molecule_mask[upper_i] | molecule_mask[upper_j]
    current_nonbaseline = current_upper > 0.5 + 1e-12
    final_nonbaseline = final_upper > 0.5 + 1e-12
    newly_nonbaseline = final_nonbaseline & ~current_nonbaseline

    previous_edges = pd.read_csv(
        current_edges_path,
        sep="\t",
        compression="gzip",
        dtype={"node_a": str, "node_b": str},
    )
    expected_pairs = pd.DataFrame(
        {
            "node_a": np.asarray(symbols, dtype=object)[upper_i],
            "node_b": np.asarray(symbols, dtype=object)[upper_j],
        }
    )
    if not previous_edges[["node_a", "node_b"]].equals(expected_pairs):
        raise ValueError(
            "Previous HPA edge table is not in strict upper-triangle order"
        )
    previous_edges = previous_edges.rename(
        columns={
            "final_posterior": "posterior_after_hpa",
            "posterior_increase": "hpa_posterior_increase",
        }
    )
    support_by_pair = supported.set_index(["node_a", "node_b"])
    support_columns = [
        "directed_record_count",
        "directional_interactions",
        "omnipath_bidirectional",
        "any_stimulation_annotation",
        "any_inhibition_annotation",
        "any_unsigned_annotation",
        "stimulation_and_inhibition_both_present",
        "any_consensus_direction",
        "any_consensus_stimulation",
        "any_consensus_inhibition",
        "omnipath_curation_effort",
        "curation_effort_sum_across_directed_records",
        "curation_effort_values",
        "unique_resource_count",
        "resources",
        "unique_reference_count",
        "references",
        "source_uniprot_accessions",
        "target_uniprot_accessions",
        "endpoint_mapping_methods",
        "omnipath_Tq",
        "omnipath_likelihood",
        "omnipath_bayes_factor",
        "omnipath_log_bayes_factor",
    ]
    all_pair_index = pd.MultiIndex.from_frame(
        previous_edges[["node_a", "node_b"]]
    )
    aligned_support = support_by_pair.reindex(all_pair_index)
    previous_edges["omnipath_evidence_observed"] = observed_upper
    for column in support_columns:
        values = aligned_support[column].reset_index(drop=True)
        if column in {
            "directed_record_count",
            "unique_resource_count",
            "unique_reference_count",
        }:
            values = values.fillna(0).astype(int)
        elif column in {
            "omnipath_bidirectional",
            "any_stimulation_annotation",
            "any_inhibition_annotation",
            "any_unsigned_annotation",
            "stimulation_and_inhibition_both_present",
            "any_consensus_direction",
            "any_consensus_stimulation",
            "any_consensus_inhibition",
        }:
            values = values.fillna(False).astype(bool)
        elif column == "omnipath_likelihood":
            values = values.fillna(NO_EDGE_LIKELIHOOD).astype(float)
        elif column == "omnipath_bayes_factor":
            values = values.fillna(1.0).astype(float)
        elif column == "omnipath_log_bayes_factor":
            values = values.fillna(0.0).astype(float)
        elif column in {
            "omnipath_curation_effort",
            "curation_effort_sum_across_directed_records",
            "omnipath_Tq",
        }:
            values = values.astype(float)
        else:
            values = values.fillna("").astype(str)
        previous_edges[column] = values.to_numpy()
    previous_edges["pre_omnipath_posterior"] = current_upper
    previous_edges["final_posterior"] = final_upper
    previous_edges["omnipath_posterior_increase"] = (
        final_upper - current_upper
    )
    previous_edges["newly_nonbaseline_from_omnipath"] = (
        newly_nonbaseline
    )
    all_edges_path = output_dir / "combined_undirected_edges.tsv.gz"
    with gzip.open(
        all_edges_path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        previous_edges.to_csv(
            handle,
            sep="\t",
            index=False,
            float_format="%.9g",
        )
    previous_edges.sort_values(
        [
            "final_posterior",
            "omnipath_bayes_factor",
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

    backend_columns = [
        "node_a",
        "node_b",
        "omnipath_curation_effort",
        "omnipath_Tq",
        "omnipath_likelihood",
        "omnipath_bayes_factor",
        "omnipath_log_bayes_factor",
        "directed_record_count",
        "directional_interactions",
        "omnipath_bidirectional",
        "any_stimulation_annotation",
        "any_inhibition_annotation",
        "any_unsigned_annotation",
        "unique_resource_count",
        "resources",
        "unique_reference_count",
        "references",
    ]
    backend = supported[backend_columns].copy()
    backend.insert(
        0,
        "edge_id",
        backend["node_a"] + "|" + backend["node_b"],
    )
    backend_path = backend_dir / "omnipath_core_bf_gt1.tsv.gz"
    with gzip.open(
        backend_path,
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        backend.to_csv(
            handle,
            sep="\t",
            index=False,
            float_format="%.9g",
        )

    raw_nodes = set(raw["source_genesymbol"]) | set(
        raw["target_genesymbol"]
    )
    in_universe_degree = Counter(
        supported["node_a"].tolist() + supported["node_b"].tolist()
    )
    mapping_rows = []
    source_symbol_counts = Counter(raw["source_genesymbol"])
    target_symbol_counts = Counter(raw["target_genesymbol"])
    source_accession_counts = Counter(raw["source"])
    target_accession_counts = Counter(raw["target"])
    for row in protein_index.to_dict(orient="records"):
        symbol = str(row["symbol"])
        accession = clean(row["mapped_uniprot"])
        mapping_rows.append(
            {
                "symbol": symbol,
                "mapped_uniprot": accession,
                "exact_symbol_present_anywhere": symbol in raw_nodes,
                "source_rows_by_exact_symbol": source_symbol_counts[symbol],
                "target_rows_by_exact_symbol": target_symbol_counts[symbol],
                "source_rows_by_exact_uniprot": (
                    source_accession_counts[accession] if accession else 0
                ),
                "target_rows_by_exact_uniprot": (
                    target_accession_counts[accession] if accession else 0
                ),
                "in_universe_undirected_degree": (
                    in_universe_degree[symbol]
                ),
            }
        )
    mapping_path = processed_dir / "node_omnipath_coverage.tsv"
    pd.DataFrame(mapping_rows).to_csv(
        mapping_path,
        sep="\t",
        index=False,
    )

    def aligned_matrix(path: Path) -> np.ndarray:
        frame = pd.read_csv(
            path,
            sep="\t",
            index_col=0,
        ).reindex(index=symbols, columns=symbols)
        if frame.isna().any().any():
            raise ValueError(f"Matrix does not align to universe: {path}")
        return frame.to_numpy(dtype=float)

    localization = aligned_matrix(localization_matrix_path)
    kinase = aligned_matrix(kinase_matrix_path)
    string_observed = aligned_matrix(string_matrix_path)
    hpa = aligned_matrix(hpa_matrix_path)
    supported_mask = observed_upper
    overlap = {
        "mpkccd_localization": int(
            (
                supported_mask
                & (localization[upper_i, upper_j] > 0.5 + 1e-12)
            ).sum()
        ),
        "kinase_predictor": int(
            (
                supported_mask
                & (kinase[upper_i, upper_j] > 1e-12)
            ).sum()
        ),
        "string_v12": int(
            (
                supported_mask
                & (string_observed[upper_i, upper_j] > 0)
            ).sum()
        ),
        "hpa_primary": int(
            (
                supported_mask
                & (hpa[upper_i, upper_j] > 1e-12)
            ).sum()
        ),
        "previously_nonbaseline_under_any_stream": int(
            (supported_mask & current_nonbaseline).sum()
        ),
        "newly_nonbaseline_from_omnipath": int(
            newly_nonbaseline.sum()
        ),
    }

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
        "unsupported_pairs_exactly_unchanged": bool(
            np.array_equal(
                final_upper[~observed_upper],
                current_upper[~observed_upper],
            )
        ),
        "molecule_incident_pairs_exactly_unchanged": bool(
            np.array_equal(
                final_upper[molecule_pair],
                current_upper[molecule_pair],
            )
        ),
        "only_protein_pairs_receive_omnipath_evidence": bool(
            not np.any(observed_upper & molecule_pair)
        ),
        "all_supported_pairs_have_bf_gt1": bool(
            (supported["omnipath_bayes_factor"] > 1.0).all()
        ),
        "supported_table_has_unique_undirected_pairs": bool(
            not supported.duplicated(["node_a", "node_b"]).any()
        ),
        "raw_query_contains_only_directed_records": bool(
            raw["is_directed_bool"].all()
        ),
        "node_order_preserved": bool(
            current.index.astype(str).tolist() == symbols
            and current.columns.astype(str).tolist() == symbols
        ),
        "unique_pairs_partition_total": bool(
            int(final_nonbaseline.sum())
            + int(
                np.isclose(
                    final_upper,
                    0.5,
                    atol=1e-12,
                    rtol=0.0,
                ).sum()
            )
            == len(final_upper)
        ),
    }
    if not all(validation.values()):
        raise AssertionError(f"Validation failed: {validation}")

    query_metadata = {
        "download_date": DOWNLOAD_DATE,
        "api_query": API_QUERY,
        "dataset": DATASET,
        "organism_ncbi_taxonomy_id": ORGANISM,
        "interaction_type": INTERACTION_TYPE,
        "entity_type": ENTITY_TYPE,
        "directed_query": True,
        "model_directionality": "undirected",
        "sign_used_in_model": False,
        "license_parameter": "academic",
        "raw_sha256": sha256_file(raw_path),
        "raw_row_count": int(len(raw)),
    }
    metadata_path = processed_dir / "query_metadata.json"
    metadata_path.write_text(
        json.dumps(query_metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    sign_counts = {
        "pairs_with_any_stimulation_annotation": int(
            supported["any_stimulation_annotation"].sum()
        ),
        "pairs_with_any_inhibition_annotation": int(
            supported["any_inhibition_annotation"].sum()
        ),
        "pairs_with_any_unsigned_annotation": int(
            supported["any_unsigned_annotation"].sum()
        ),
        "pairs_with_both_stimulation_and_inhibition": int(
            supported[
                "stimulation_and_inhibition_both_present"
            ].sum()
        ),
        "pairs_with_bidirectional_records": int(
            supported["omnipath_bidirectional"].sum()
        ),
        "pairs_with_any_consensus_direction": int(
            supported["any_consensus_direction"].sum()
        ),
    }
    baseline_final = int(
        np.isclose(
            final_upper,
            0.5,
            atol=1e-12,
            rtol=0.0,
        ).sum()
    )
    summary = {
        "download_date": DOWNLOAD_DATE,
        "omnipath_dataset": DATASET,
        "organism_ncbi_taxonomy_id": ORGANISM,
        "interaction_type": INTERACTION_TYPE,
        "entity_type": ENTITY_TYPE,
        "node_count": n,
        "protein_node_count": len(protein_symbols),
        "molecule_node_count": len(molecule_symbols),
        "possible_unique_undirected_pairs": int(len(final_upper)),
        "raw_directed_interaction_rows": int(len(raw)),
        "full_mouse_background_unique_undirected_pairs": int(
            len(background)
        ),
        "full_mouse_background_curation_effort_quantiles": {
            str(key): float(value)
            for key, value in background[
                "omnipath_curation_effort"
            ]
            .quantile([0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1])
            .items()
        },
        "q": Q,
        "Tq_curation_effort": t_q,
        "in_universe_directed_rows": int(len(in_universe)),
        "in_universe_unique_undirected_pairs": int(len(supported)),
        "symbol_endpoint_mappings": int(
            (
                raw["source_mapping_method"].eq(
                    "exact_mouse_gene_symbol"
                ).sum()
                + raw["target_mapping_method"].eq(
                    "exact_mouse_gene_symbol"
                ).sum()
            )
        ),
        "uniprot_fallback_endpoint_mappings": int(
            (
                raw["source_mapping_method"].eq(
                    "exact_mouse_uniprot_fallback"
                ).sum()
                + raw["target_mapping_method"].eq(
                    "exact_mouse_uniprot_fallback"
                ).sum()
            )
        ),
        "protein_nodes_incident_to_omnipath_pair": int(
            len(set(supported["node_a"]) | set(supported["node_b"]))
        ),
        "sign_and_direction_audit": sign_counts,
        "direction_and_sign_model_rule": (
            "All source-target directions and stimulation/inhibition "
            "annotations are retained in the audit table but ignored in the "
            "modeled graph. Every supported pair is stored once in node-"
            "universe order and receives one undirected Bayes factor."
        ),
        "pair_collapse_rule": (
            "Rows mapping to the same unordered node pair are collapsed; "
            "the maximum API curation_effort is used as the conservative "
            "score. Resources, references, directions, and signs are unioned "
            "for audit."
        ),
        "background_definition": (
            "All unique unordered non-self gene-symbol pairs in the full "
            "mouse OmniPath core post-translational query, before node-"
            "universe restriction; score is maximum curation_effort across "
            "directional rows for each unordered pair."
        ),
        "likelihood_formula": (
            "L_OmniPath = 0.5 + 0.5 * "
            "(1 - exp(-0.5 * (curation_effort / Tq)^2))"
        ),
        "bayes_factor_formula": "BF_OmniPath = L_OmniPath / 0.5",
        "integration_formula": (
            "logit(P_final) = logit(P_after_HPA) + log(BF_OmniPath)"
        ),
        "missing_evidence_rule": (
            "A pair absent from the mapped OmniPath query receives "
            "likelihood 0.5, BF=1, and is left exactly unchanged."
        ),
        "supported_bayes_factor_minimum": float(
            supported["omnipath_bayes_factor"].min()
        ),
        "supported_bayes_factor_median": float(
            supported["omnipath_bayes_factor"].median()
        ),
        "supported_bayes_factor_maximum": float(
            supported["omnipath_bayes_factor"].max()
        ),
        "overlap_with_prior_streams": overlap,
        "nonbaseline_edge_count_before_omnipath": int(
            current_nonbaseline.sum()
        ),
        "nonbaseline_edge_count_after_omnipath": int(
            final_nonbaseline.sum()
        ),
        "newly_nonbaseline_edge_count": int(newly_nonbaseline.sum()),
        "pairs_at_final_baseline_probability": baseline_final,
        "final_minimum_off_diagonal_probability": float(
            final_upper.min()
        ),
        "final_maximum_off_diagonal_probability": float(
            final_upper.max()
        ),
        "dependence_caveat": (
            "OmniPath integrates many literature and database resources. "
            "Some resources can overlap evidence already represented by "
            "STRING or kinase/phosphosite streams, so the independent "
            "Bayesian update can overstate evidence when sources are not "
            "conditionally independent. The sparse OmniPath factor table is "
            "retained so this stream can be toggled off."
        ),
        "validation": validation,
        "source_urls": {
            "api_query": API_QUERY,
            "api_query_documentation": API_QUERY_DOCS,
            "omnipath_home": OMNIPATH_HOME,
            "interaction_documentation": OMNIPATH_R_DOCS,
        },
        "input_sha256": {
            path.name: sha256_file(path)
            for path in (
                universe_path,
                protein_index_path,
                raw_path,
                current_matrix_path,
                current_edges_path,
            )
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    readme = f"""OMNIPATH CORE UNDIRECTED EDGE INTEGRATION

PURPOSE
This directory extends the current {n}-node UNDIRECTED adjacency model with
the OmniPath core mouse post-translational signaling network downloaded on
{DOWNLOAD_DATE}. The API records are directed and can be activating,
inhibitory, or unsigned. Those annotations are preserved for audit but are
not used as graph direction or edge sign in this analysis.

SCOPE
- Dataset: {DATASET} (curated core only; no pathwayextra, kinaseextra, or
  ligrecextra expansion).
- Organism query: mouse, NCBI taxonomy {ORGANISM}. OmniPath documents mouse
  interactions as orthology-translated from its human knowledge base.
- Type: {INTERACTION_TYPE}; entity type: protein.
- Exact mouse gene symbols are mapped first. Exact mouse UniProt accessions
  from the project protein index are used only as a fallback.
- Small molecules receive no OmniPath evidence.

UNDIRECTED COLLAPSE
All {len(in_universe):,} mapped directed records are collapsed onto
{len(supported):,} unique unordered node pairs. If both A->B and B->A are
present, the pair is still counted once. The maximum curation_effort across
directional records is the conservative pair score; directions, signs,
resources, and references are unioned in the audit table.

BACKGROUND AND BAYES FACTOR
T_q is q75 of maximum curation_effort across all {len(background):,} unique
unordered non-self pairs in the full mouse query before restriction to the
project node universe. Here, T_q = {t_q:g}.

For an observed OmniPath pair with score s:

  L = 0.5 + 0.5 * (1 - exp(-0.5 * (s / T_q)^2))
  BF_OmniPath = L / 0.5
  logit(P_final) = logit(P_after_HPA) + log(BF_OmniPath)

The zero-anchored form makes every curated positive score weakly or strongly
supportive while retaining the same 0.5 neutral likelihood and a maximum BF
of 2. A pair absent from OmniPath has L=0.5 and BF=1 and remains unchanged.

RESULT
- OmniPath-supported unique undirected pairs: {len(supported):,}
- Previously non-baseline supported pairs: {int((supported['pre_omnipath_posterior'] > 0.5 + 1e-12).sum()):,}
- Newly non-baseline pairs: {int(supported['newly_nonbaseline_edge'].sum()):,}
- Cumulative non-baseline edges: {int(final_nonbaseline.sum()):,}
- Remaining baseline pairs: {baseline_final:,}

FILES
- combined_adjacency_matrix.tsv: final symmetric {n} x {n} matrix.
- combined_undirected_edges.tsv.gz: all {len(final_upper):,} unique pairs.
- omnipath_supported_undirected_edges.tsv: complete supported-pair audit.
- omnipath_evidence_observed_matrix.tsv: binary coverage.
- omnipath_curation_effort_matrix.tsv: conservative pair score.
- omnipath_likelihood_matrix.tsv: edge likelihood for this stream.
- omnipath_log_bayes_factor_matrix.tsv: additive log evidence.
- top_combined_undirected_edges.tsv: top 1,000 final edges.
- analysis_summary.json: definitions, counts, hashes, and validation.

PROCESSED AUDITS
- {background_path}
- {mapping_path}
- {metadata_path}

BACKEND FACTOR TABLE
- {backend_path}
Absent rows in this sparse table mean BF=1 and log BF=0.

DEPENDENCE CAVEAT
OmniPath integrates many literature and database resources. Some evidence may
overlap sources already represented by STRING or kinase/phosphosite streams.
Treating all streams as conditionally independent can overstate combined
support; retain OmniPath as a separately toggleable stream.

OFFICIAL SOURCES
{API_QUERY}
{API_QUERY_DOCS}
{OMNIPATH_HOME}
{OMNIPATH_R_DOCS}
"""
    readme_path = output_dir / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")
    summary["output_sha256"] = {
        path.name: sha256_file(path)
        for path in (
            background_path,
            mapping_path,
            metadata_path,
            supported_path,
            final_matrix_path,
            all_edges_path,
            backend_path,
            readme_path,
        )
    }
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
