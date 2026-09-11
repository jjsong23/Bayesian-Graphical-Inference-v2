#!/usr/bin/env python3
"""Integrate STRING v12.0 mouse associations into the undirected edge model.

The script maps protein nodes to STRING proteins, retains functional
associations whose two endpoints are in the current node universe, converts
STRING's combined confidence score to a Bayes factor using STRING's 0.041
prior, and updates the existing localization-plus-KinasePredictor adjacency
matrix. Missing STRING relationships are treated as missing evidence, not as
evidence against an edge.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


STRING_VERSION = "12.0"
STRING_SPECIES = "10090"
STRING_PRIOR = 0.041
MEDIUM_CONFIDENCE = 0.4
HIGH_CONFIDENCE = 0.7
STRING_PAPER = "https://pubmed.ncbi.nlm.nih.gov/39558183/"
STRING_SCORE_DOCS = "https://string-db.org/help/scores/"
STRING_MOUSE_DOWNLOADS = (
    "https://string-db.org/cgi/download?species_text=Mus+musculus"
)

CHANNEL_COLUMNS = [
    "neighborhood_score",
    "fusion_score",
    "cooccurrence_score",
    "coexpression_score",
    "experimental_score",
    "database_score",
    "textmining_score",
]


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


def stable_expit(value: np.ndarray | float) -> np.ndarray | float:
    array = np.asarray(value, dtype=float)
    result = 1.0 / (1.0 + np.exp(-np.clip(array, -700.0, 700.0)))
    if np.isscalar(value):
        return float(result)
    return result


def logit(value: np.ndarray | float) -> np.ndarray | float:
    array = np.asarray(value, dtype=float)
    clipped = np.clip(array, 1e-15, 1.0 - 1e-15)
    result = np.log(clipped / (1.0 - clipped))
    if np.isscalar(value):
        return float(result)
    return result


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "-"} else text


def read_string_info(path: Path) -> dict[str, dict[str, object]]:
    info: dict[str, dict[str, object]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected = [
            "#string_protein_id",
            "preferred_name",
            "protein_size",
            "annotation",
        ]
        if header != expected:
            raise ValueError(f"Unexpected STRING info header: {header}")
        for line in handle:
            string_id, preferred, size, annotation = line.rstrip("\n").split(
                "\t", 3
            )
            info[string_id] = {
                "preferred_name": preferred,
                "protein_size": int(size),
                "annotation": annotation,
            }
    return info


def collect_alias_hits(
    path: Path,
    aliases_of_interest: set[str],
) -> dict[str, dict[str, set[str]]]:
    wanted = {alias.casefold() for alias in aliases_of_interest if alias}
    hits: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["#string_protein_id", "alias", "source"]:
            raise ValueError(f"Unexpected STRING alias header: {header}")
        for line in handle:
            string_id, alias, source = line.rstrip("\n").split("\t", 2)
            key = alias.casefold()
            if key in wanted:
                hits[key][string_id].add(source)
    return hits


def choose_mapping(
    symbol: str,
    accessions: list[str],
    info: dict[str, dict[str, object]],
    alias_hits: dict[str, dict[str, set[str]]],
) -> dict[str, object]:
    """Choose an auditable STRING identifier without arbitrary synonym picks."""

    accession_candidates: dict[str, set[str]] = defaultdict(set)
    for accession in accessions:
        for string_id, sources in alias_hits.get(
            accession.casefold(), {}
        ).items():
            accession_candidates[string_id].update(sources)

    symbol_candidates = {
        string_id: set(sources)
        for string_id, sources in alias_hits.get(
            symbol.casefold(), {}
        ).items()
    }
    preferred_exact = {
        string_id
        for string_id, metadata in info.items()
        if str(metadata["preferred_name"]).casefold() == symbol.casefold()
    }

    method = ""
    selected = ""
    candidates: dict[str, set[str]] = accession_candidates

    accession_preferred = set(accession_candidates).intersection(preferred_exact)
    if len(accession_preferred) == 1:
        selected = next(iter(accession_preferred))
        method = "uniprot_accession_plus_preferred_name"
    elif len(accession_candidates) == 1:
        selected = next(iter(accession_candidates))
        method = "unique_uniprot_accession"
    elif len(accession_candidates) > 1:
        # Exact gene-name agreement is required when a UniProt accession is
        # shared by several STRING/Ensembl protein records.
        candidates = accession_candidates
        method = "ambiguous_uniprot_accession_unmapped"
    else:
        symbol_preferred = set(symbol_candidates).intersection(preferred_exact)
        candidates = symbol_candidates
        if len(symbol_preferred) == 1:
            selected = next(iter(symbol_preferred))
            method = "gene_symbol_preferred_name_fallback"
        elif len(preferred_exact) == 1:
            selected = next(iter(preferred_exact))
            method = "preferred_name_fallback"
            candidates = {selected: set()}
        elif len(symbol_candidates) == 1:
            selected = next(iter(symbol_candidates))
            method = "unique_gene_symbol_alias_fallback"
        elif len(symbol_candidates) > 1:
            method = "ambiguous_gene_symbol_unmapped"
        else:
            method = "no_string_mapping"

    selected_sources = sorted(candidates.get(selected, set()))
    candidate_ids = sorted(candidates)
    return {
        "symbol": symbol,
        "selected_uniprot": ";".join(accessions),
        "string_id": selected,
        "string_preferred_name": (
            str(info[selected]["preferred_name"]) if selected else ""
        ),
        "mapping_status": "mapped" if selected else "unmapped",
        "mapping_method": method,
        "candidate_count": len(candidate_ids),
        "candidate_string_ids": ";".join(candidate_ids),
        "selected_alias_sources": ";".join(selected_sources),
    }


def build_mapping(
    protein_index: pd.DataFrame,
    info: dict[str, dict[str, object]],
    alias_hits: dict[str, dict[str, set[str]]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in protein_index.to_dict(orient="records"):
        symbol = clean(row["symbol"])
        accessions = list(
            dict.fromkeys(
                accession
                for accession in (
                    clean(row.get("mapped_uniprot")),
                    clean(row.get("selected_uniprot")),
                )
                if accession
            )
        )
        rows.append(
            choose_mapping(symbol, accessions, info, alias_hits)
        )
    mapping = pd.DataFrame(rows)
    mapped = mapping.loc[mapping["mapping_status"].eq("mapped")]
    collisions = mapped.loc[
        mapped["string_id"].duplicated(keep=False)
    ].sort_values(["string_id", "symbol"])
    if not collisions.empty:
        detail = collisions[["symbol", "string_id"]].to_dict(orient="records")
        raise ValueError(f"Multiple nodes map to one STRING protein: {detail}")
    return mapping


def read_in_universe_links(
    path: Path,
    string_id_to_symbol: dict[str, str],
    symbol_to_index: dict[str, int],
) -> tuple[dict[tuple[str, str], dict[str, object]], dict[str, int]]:
    """Stream the large detailed network and retain only in-universe pairs."""

    records: dict[tuple[str, str], dict[str, object]] = {}
    stats = {
        "raw_network_rows": 0,
        "rows_with_both_endpoints_mapped": 0,
        "self_rows_excluded": 0,
        "duplicate_unordered_rows": 0,
        "duplicate_rows_with_score_disagreement": 0,
        "replacement_rows_with_higher_combined_score": 0,
    }
    expected_header = [
        "protein1",
        "protein2",
        "neighborhood",
        "fusion",
        "cooccurence",
        "coexpression",
        "experimental",
        "database",
        "textmining",
        "combined_score",
    ]
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().strip().split()
        if header != expected_header:
            raise ValueError(f"Unexpected STRING links header: {header}")
        for line in handle:
            stats["raw_network_rows"] += 1
            fields = line.split()
            if len(fields) != 10:
                raise ValueError(f"Malformed STRING link row: {line[:200]}")
            string_a, string_b = fields[:2]
            node_a = string_id_to_symbol.get(string_a)
            node_b = string_id_to_symbol.get(string_b)
            if node_a is None or node_b is None:
                continue
            stats["rows_with_both_endpoints_mapped"] += 1
            if node_a == node_b:
                stats["self_rows_excluded"] += 1
                continue
            if symbol_to_index[node_a] > symbol_to_index[node_b]:
                node_a, node_b = node_b, node_a
                string_a, string_b = string_b, string_a

            scores = [int(value) / 1000.0 for value in fields[2:]]
            record = {
                "node_a": node_a,
                "node_b": node_b,
                "string_id_a": string_a,
                "string_id_b": string_b,
                "neighborhood_score": scores[0],
                "fusion_score": scores[1],
                "cooccurrence_score": scores[2],
                "coexpression_score": scores[3],
                "experimental_score": scores[4],
                "database_score": scores[5],
                "textmining_score": scores[6],
                "string_combined_score": scores[7],
            }
            key = (node_a, node_b)
            previous = records.get(key)
            if previous is None:
                records[key] = record
                continue

            stats["duplicate_unordered_rows"] += 1
            score_columns = CHANNEL_COLUMNS + ["string_combined_score"]
            disagrees = any(
                not math.isclose(
                    float(previous[column]),
                    float(record[column]),
                    abs_tol=1e-12,
                )
                for column in score_columns
            )
            if disagrees:
                stats["duplicate_rows_with_score_disagreement"] += 1
            if float(record["string_combined_score"]) > float(
                previous["string_combined_score"]
            ):
                records[key] = record
                stats["replacement_rows_with_higher_combined_score"] += 1
    return records, stats


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()

    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    protein_index_path = (
        project
        / "data/edge_characterization/kinase_predictor/phosphosite_database"
        / "protein_index.tsv"
    )
    current_matrix_path = (
        project
        / "results/edge_characterization/localization_kinase_predictor"
        / "combined_adjacency_matrix.tsv"
    )
    raw_dir = (
        project / "data/edge_characterization/string/v12.0/raw"
    )
    info_path = raw_dir / "10090.protein.info.v12.0.txt.gz"
    aliases_path = raw_dir / "10090.protein.aliases.v12.0.txt.gz"
    links_path = raw_dir / "10090.protein.links.detailed.v12.0.txt.gz"
    processed_dir = (
        project / "data/edge_characterization/string/v12.0/processed"
    )
    output_dir = (
        project
        / "results/edge_characterization"
        / "localization_kinase_predictor_string"
    )
    processed_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    required = [
        universe_path,
        protein_index_path,
        current_matrix_path,
        info_path,
        aliases_path,
        links_path,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    if universe["symbol"].duplicated().any():
        raise ValueError("Node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}
    protein_symbols = set(
        universe.loc[universe["node_type"].eq("protein"), "symbol"]
    )
    molecule_symbols = set(
        universe.loc[universe["node_type"].eq("molecule"), "symbol"]
    )
    n = len(symbols)
    if n < 2:
        raise ValueError(f"Expected at least two nodes, found {n}")

    protein_index = pd.read_csv(
        protein_index_path, sep="\t", dtype=str
    ).fillna("")
    if set(protein_index["symbol"]) != protein_symbols:
        missing = sorted(protein_symbols.difference(protein_index["symbol"]))
        extra = sorted(set(protein_index["symbol"]).difference(protein_symbols))
        raise ValueError(
            f"Protein index mismatch; missing={missing}, extra={extra}"
        )

    info = read_string_info(info_path)
    aliases_of_interest = set(protein_index["symbol"])
    aliases_of_interest.update(protein_index["mapped_uniprot"])
    aliases_of_interest.update(protein_index["selected_uniprot"])
    alias_hits = collect_alias_hits(aliases_path, aliases_of_interest)
    mapping = build_mapping(protein_index, info, alias_hits)
    mapping = mapping.sort_values(
        "symbol",
        key=lambda series: series.map(symbol_to_index),
        kind="stable",
    )
    mapping_path = processed_dir / "node_to_string_mapping.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)
    mapping.loc[mapping["mapping_status"].eq("unmapped")].to_csv(
        processed_dir / "unmapped_protein_nodes.tsv",
        sep="\t",
        index=False,
    )

    mapped = mapping.loc[mapping["mapping_status"].eq("mapped")]
    string_id_to_symbol = dict(zip(mapped["string_id"], mapped["symbol"]))
    link_records, link_stats = read_in_universe_links(
        links_path,
        string_id_to_symbol,
        symbol_to_index,
    )
    supported = pd.DataFrame(link_records.values())
    if supported.empty:
        raise ValueError("No in-universe STRING associations were retained")
    supported = supported.sort_values(
        ["string_combined_score", "node_a", "node_b"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)

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

    n = len(symbols)
    observed_matrix = np.zeros((n, n), dtype=np.uint8)
    score_matrix = np.zeros((n, n), dtype=float)
    log_bf_matrix = np.zeros((n, n), dtype=float)
    prior_log_odds = float(logit(STRING_PRIOR))

    log_bfs: list[float] = []
    current_probabilities: list[float] = []
    final_probabilities: list[float] = []
    for record in supported.to_dict(orient="records"):
        i = symbol_to_index[str(record["node_a"])]
        j = symbol_to_index[str(record["node_b"])]
        score = float(record["string_combined_score"])
        if not 0.0 < score < 1.0:
            raise ValueError(f"Invalid STRING probability {score}: {record}")
        string_log_bf = float(logit(score) - prior_log_odds)
        current_probability = float(current_matrix[i, j])
        final_probability = float(
            stable_expit(logit(current_probability) + string_log_bf)
        )
        observed_matrix[i, j] = observed_matrix[j, i] = 1
        score_matrix[i, j] = score_matrix[j, i] = score
        log_bf_matrix[i, j] = log_bf_matrix[j, i] = string_log_bf
        log_bfs.append(string_log_bf)
        current_probabilities.append(current_probability)
        final_probabilities.append(final_probability)

    final_matrix = stable_expit(logit(np.where(
        np.eye(n, dtype=bool),
        0.5,
        current_matrix,
    )) + log_bf_matrix)
    final_matrix[observed_matrix == 0] = current_matrix[
        observed_matrix == 0
    ]
    np.fill_diagonal(final_matrix, 0.0)

    supported["string_log_bayes_factor"] = log_bfs
    supported["string_bayes_factor"] = np.exp(
        np.clip(np.asarray(log_bfs), -700.0, 700.0)
    )
    supported["pre_string_posterior"] = current_probabilities
    supported["final_posterior"] = final_probabilities
    supported["posterior_increase"] = (
        supported["final_posterior"] - supported["pre_string_posterior"]
    )
    supported["medium_confidence_string"] = (
        supported["string_combined_score"] >= MEDIUM_CONFIDENCE
    )
    supported["high_confidence_string"] = (
        supported["string_combined_score"] >= HIGH_CONFIDENCE
    )

    supported_path = output_dir / "string_supported_undirected_edges.tsv"
    supported.to_csv(
        supported_path, sep="\t", index=False, float_format="%.9g"
    )
    write_matrix(
        output_dir / "string_evidence_observed_matrix.tsv",
        symbols,
        observed_matrix,
    )
    write_matrix(
        output_dir / "string_combined_score_matrix.tsv",
        symbols,
        score_matrix,
    )
    write_matrix(
        output_dir / "string_log_bayes_factor_matrix.tsv",
        symbols,
        log_bf_matrix,
    )
    write_matrix(
        output_dir / "combined_adjacency_matrix.tsv",
        symbols,
        final_matrix,
    )

    support_by_pair = {
        (str(row["node_a"]), str(row["node_b"])): row
        for row in supported.to_dict(orient="records")
    }
    upper_i, upper_j = np.triu_indices(n, k=1)
    all_rows: list[dict[str, object]] = []
    for i, j in zip(upper_i, upper_j):
        node_a, node_b = symbols[i], symbols[j]
        string_row = support_by_pair.get((node_a, node_b))
        row: dict[str, object] = {
            "node_a": node_a,
            "node_b": node_b,
            "pre_string_posterior": current_matrix[i, j],
            "string_evidence_observed": string_row is not None,
            "string_id_a": "",
            "string_id_b": "",
            "string_combined_score": "",
            **{column: "" for column in CHANNEL_COLUMNS},
            "string_log_bayes_factor": log_bf_matrix[i, j],
            "string_bayes_factor": (
                math.exp(log_bf_matrix[i, j])
                if string_row is not None
                else 1.0
            ),
            "final_posterior": final_matrix[i, j],
            "posterior_increase": (
                final_matrix[i, j] - current_matrix[i, j]
            ),
            "medium_confidence_string": False,
            "high_confidence_string": False,
        }
        if string_row is not None:
            for column in (
                "string_id_a",
                "string_id_b",
                "string_combined_score",
                *CHANNEL_COLUMNS,
                "medium_confidence_string",
                "high_confidence_string",
            ):
                row[column] = string_row[column]
        all_rows.append(row)
    all_edges = pd.DataFrame(all_rows)
    all_edges_path = output_dir / "combined_undirected_edges.tsv.gz"
    with gzip.open(
        all_edges_path, "wt", encoding="utf-8", newline=""
    ) as handle:
        all_edges.to_csv(
            handle, sep="\t", index=False, float_format="%.9g"
        )
    all_edges.sort_values(
        [
            "final_posterior",
            "string_evidence_observed",
            "string_combined_score",
            "node_a",
            "node_b",
        ],
        ascending=[False, False, False, True, True],
        kind="stable",
    ).head(1000).to_csv(
        output_dir / "top_combined_undirected_edges.tsv",
        sep="\t",
        index=False,
        float_format="%.9g",
    )

    supported_mask = observed_matrix.astype(bool)
    unsupported_mask = ~supported_mask
    np.fill_diagonal(unsupported_mask, False)
    molecule_indices = np.asarray(
        [symbol_to_index[symbol] for symbol in molecule_symbols], dtype=int
    )
    molecule_pair_mask = np.zeros((n, n), dtype=bool)
    molecule_pair_mask[molecule_indices, :] = True
    molecule_pair_mask[:, molecule_indices] = True
    np.fill_diagonal(molecule_pair_mask, False)

    validation = {
        "matrix_shape": list(final_matrix.shape),
        "matrix_symmetric": bool(
            np.allclose(final_matrix, final_matrix.T, atol=1e-12)
        ),
        "diagonal_all_zero": bool(
            np.allclose(np.diag(final_matrix), 0.0, atol=1e-12)
        ),
        "all_probabilities_in_zero_one": bool(
            np.all((final_matrix >= 0.0) & (final_matrix <= 1.0))
        ),
        "unsupported_pairs_exactly_unchanged": bool(
            np.array_equal(
                final_matrix[unsupported_mask],
                current_matrix[unsupported_mask],
            )
        ),
        "all_molecule_incident_pairs_exactly_unchanged": bool(
            np.array_equal(
                final_matrix[molecule_pair_mask],
                current_matrix[molecule_pair_mask],
            )
        ),
        "string_supported_pairs_are_protein_protein": bool(
            set(supported["node_a"]).union(supported["node_b"]).issubset(
                protein_symbols
            )
        ),
        "duplicate_unordered_supported_pairs": int(
            supported.duplicated(["node_a", "node_b"]).sum()
        ),
    }
    required_validation_checks = [
        validation["matrix_shape"] == [n, n],
        validation["matrix_symmetric"],
        validation["diagonal_all_zero"],
        validation["all_probabilities_in_zero_one"],
        validation["unsupported_pairs_exactly_unchanged"],
        validation["all_molecule_incident_pairs_exactly_unchanged"],
        validation["string_supported_pairs_are_protein_protein"],
        validation["duplicate_unordered_supported_pairs"] == 0,
    ]
    if not all(required_validation_checks):
        raise AssertionError(f"Validation failed: {validation}")

    current_upper = current_matrix[upper_i, upper_j]
    final_upper = final_matrix[upper_i, upper_j]
    summary = {
        "string_version": STRING_VERSION,
        "species_taxon": STRING_SPECIES,
        "network_type": "functional association",
        "edge_directionality": "undirected",
        "node_count": n,
        "protein_node_count": len(protein_symbols),
        "molecule_node_count": len(molecule_symbols),
        "protein_nodes_mapped_to_string": int(len(mapped)),
        "protein_nodes_unmapped_from_string": int(
            mapping["mapping_status"].eq("unmapped").sum()
        ),
        "unmapped_symbols": mapping.loc[
            mapping["mapping_status"].eq("unmapped"), "symbol"
        ].tolist(),
        "possible_undirected_pairs": len(all_edges),
        "unique_string_supported_pairs": len(supported),
        "string_supported_medium_confidence_pairs": int(
            supported["medium_confidence_string"].sum()
        ),
        "string_supported_high_confidence_pairs": int(
            supported["high_confidence_string"].sum()
        ),
        "pairs_strengthened_by_string": int(
            (final_upper > current_upper).sum()
        ),
        "pairs_unchanged_by_string": int(
            (final_upper == current_upper).sum()
        ),
        "pre_string_minimum_off_diagonal_probability": float(
            current_upper.min()
        ),
        "pre_string_maximum_off_diagonal_probability": float(
            current_upper.max()
        ),
        "final_minimum_off_diagonal_probability": float(final_upper.min()),
        "final_maximum_off_diagonal_probability": float(final_upper.max()),
        "string_prior_probability": STRING_PRIOR,
        "string_bayes_factor_formula": (
            "BF_STRING=(s/(1-s))/(0.041/(1-0.041))"
        ),
        "integration_formula": (
            "logit(P_final)=logit(P_localization+KinasePredictor)"
            "+log(BF_STRING)"
        ),
        "missing_evidence_rule": (
            "No reported STRING relationship or unmapped endpoint is neutral; "
            "the pre-STRING posterior is retained exactly."
        ),
        "score_interpretation": (
            "STRING combined score is an approximate probability that a "
            "functional association exists, not interaction strength."
        ),
        "channel_rule": (
            "Individual STRING channels are retained for audit only. They are "
            "not added separately because they already contribute to the "
            "combined score."
        ),
        "confidence_reporting_thresholds": {
            "medium": MEDIUM_CONFIDENCE,
            "high": HIGH_CONFIDENCE,
        },
        "network_streaming_statistics": link_stats,
        "validation": validation,
        "sources": {
            "paper": STRING_PAPER,
            "score_documentation": STRING_SCORE_DOCS,
            "mouse_downloads": STRING_MOUSE_DOWNLOADS,
        },
        "input_sha256": {
            path.name: sha256_file(path)
            for path in (
                universe_path,
                protein_index_path,
                current_matrix_path,
                info_path,
                aliases_path,
                links_path,
            )
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    readme = f"""STRING v{STRING_VERSION} UNDIRECTED EDGE INTEGRATION

SCOPE
- Species: Mus musculus (NCBI taxonomy {STRING_SPECIES}).
- Network: STRING functional-association network.
- Current universe: {n} nodes ({len(protein_symbols)} proteins and
  {len(molecule_symbols)} small-molecule second messengers).
- STRING evidence is restricted to pairs whose two endpoints map to protein
  nodes already present in the node universe.
- The graph remains UNDIRECTED. STRING endpoint order does not encode
  information flow.

IDENTIFIER MAPPING
Each protein is mapped using its project-selected UniProt accession first.
When an accession is absent or shared, the node gene symbol must agree with the
STRING preferred name. A gene-symbol/preferred-name fallback is recorded when
needed. No arbitrary choice is made among unresolved synonyms.

BAYESIAN UPDATE
Let s be STRING's combined confidence score and pi={STRING_PRIOR} be STRING's
documented prior probability. The STRING evidence Bayes factor is:

    BF_STRING = (s / (1-s)) / (pi / (1-pi))

The existing localization-plus-KinasePredictor posterior is updated as:

    logit(P_final) = logit(P_current) + log(BF_STRING)

The continuous combined score is used without a hard cutoff. Medium (s >=
{MEDIUM_CONFIDENCE}) and high (s >= {HIGH_CONFIDENCE}) labels are descriptive
only.

MISSING EVIDENCE
An absent STRING row is not interpreted as evidence against an edge. STRING's
downloaded network is score-filtered, and lack of a row can also reflect
coverage. Unsupported and molecule-incident pairs therefore retain their
pre-STRING posterior exactly.

CHANNELS AND DEPENDENCE
Neighborhood, fusion, co-occurrence, co-expression, experimental, database,
and text-mining scores are retained in the edge audit. They are not added as
separate Bayes factors because STRING's combined score already integrates
them. Future evidence streams should be checked for overlap with these
channels before assuming conditional independence.

FILES
- combined_adjacency_matrix.tsv: final symmetric {n} x {n} probabilities.
- combined_undirected_edges.tsv.gz: all {len(all_edges):,} unordered pairs.
- string_supported_undirected_edges.tsv: only pairs with a STRING row.
- string_combined_score_matrix.tsv: STRING scores; zero means unreported, not
  negative evidence.
- string_evidence_observed_matrix.tsv: distinguishes reported from unreported
  pairs.
- string_log_bayes_factor_matrix.tsv: neutral zero for missing evidence.
- top_combined_undirected_edges.tsv: top 1,000 final posterior edges.
- analysis_summary.json: provenance, counts, hashes, and validation.
- ../../../../data/edge_characterization/string/v12.0/processed/
  node_to_string_mapping.tsv: identifier-mapping audit.

OFFICIAL SOURCES
- STRING 2025 paper: {STRING_PAPER}
- Score documentation: {STRING_SCORE_DOCS}
- Mouse downloads: {STRING_MOUSE_DOWNLOADS}
"""
    readme_path = output_dir / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")

    summary["output_sha256"] = {
        path.name: sha256_file(path)
        for path in (
            mapping_path,
            supported_path,
            output_dir / "combined_adjacency_matrix.tsv",
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
