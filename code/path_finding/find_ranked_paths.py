#!/usr/bin/env python3
"""Find and rank likely simple paths in the current partially directed graph.

An external mouse protein target is automatically added with
``build_target_adjacency_vector``. Edges at the neutral 0.5 baseline are
excluded by default. Ontology rules and mapped OmniPath source-target records
conservatively remove uniquely disallowed reverse traversals. Retained paths
are ranked by the product of their edge probabilities, equivalently by the sum
of ``-log(probability)`` edge costs.
"""

from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from build_target_adjacency_vector import (
    PROJECT_DEFAULT,
    TargetAdjacencyResult,
    build_target_adjacency_vector,
    resolve_target,
    safe_name,
    sha256_file,
)
from ontology_directionality import (
    RULE_CATALOG_PATH,
    apply_ontology_directionality,
    build_omnipath_direction_evidence,
    build_complete_class_pair_catalog,
    load_direction_rule_catalog,
)


CURRENT_MATRIX_RELATIVE = Path(
    "results/edge_characterization/"
    "localization_kinase_predictor_string_hpa_omnipath/"
    "combined_adjacency_matrix.tsv"
)
UNIVERSE_RELATIVE = Path("data/node_selection/node_universe_combined_nonzero.tsv")
UNIPROT_RELATIVE = Path(
    "data/edge_characterization/kinase_predictor/"
    "phosphosite_database/raw/uniprot_mouse_reference_proteome.tsv.gz"
)
OMNIPATH_DIRECTION_RAW_RELATIVE = Path(
    "data/edge_characterization/omnipath/2026-07-30/raw/"
    "omnipath_mouse_core_post_translational.tsv"
)

# These aliases are intentionally narrow. They make the motivating PKA example
# convenient without attempting an unreliable general natural-language gene
# resolver.
CURATED_START_ALIASES = {
    "pka subunit a": "Prkaca",
    "pka subunit alpha": "Prkaca",
    "pka catalytic subunit a": "Prkaca",
    "pka catalytic subunit alpha": "Prkaca",
    "pka alpha": "Prkaca",
    "pkaca": "Prkaca",
    "pka subunit b": "Prkacb",
    "pka subunit beta": "Prkacb",
    "pka catalytic subunit b": "Prkacb",
    "pka catalytic subunit beta": "Prkacb",
    "pka beta": "Prkacb",
    "pkacb": "Prkacb",
}

DEFAULT_SIGNAL_RELAY_CLASSES = (
    "receptor",
    "receptor_regulator",
    "ligand",
    "kinase",
    "phosphatase",
    "cyclase",
    "phosphodiesterase",
    "phospholipase",
    "nos",
    "gtpase",
    "gtpase_regulator",
    "kinase_phosphatase_binding",
    "second_messenger",
)


@dataclass(frozen=True)
class SearchPath:
    """One ranked path represented by matrix indices."""

    cost: float
    nodes: tuple[int, ...]


@dataclass
class RankedPathResult:
    """In-memory result returned by :func:`find_ranked_paths`."""

    start_symbol: str
    target_symbol: str
    paths: pd.DataFrame
    path_edges: pd.DataFrame
    summary: dict[str, Any]
    output_directory: Path | None
    matrix_path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("start", help="Starting node symbol, display symbol, or supported alias.")
    parser.add_argument("target", help="Existing node or external mouse gene symbol/UniProt accession.")
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_DEFAULT,
        help=f"Project directory (default: {PROJECT_DEFAULT})",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=25,
        help="Maximum number of ranked paths to return (default: 25).",
    )
    parser.add_argument(
        "--max-hops",
        type=int,
        default=6,
        help="Maximum number of edges in a path (default: 6).",
    )
    parser.add_argument(
        "--minimum-edge-probability",
        type=float,
        default=0.5,
        help="Retain edges strictly above this value (default: 0.5).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory.",
    )
    parser.add_argument(
        "--rebuild-target-extension",
        action="store_true",
        help="Recompute an external target vector even when a validated current cache exists.",
    )
    parser.add_argument(
        "--no-auto-extend",
        action="store_true",
        help="Reject a target outside the current matrix instead of extending it.",
    )
    parser.add_argument(
        "--allow-all-intermediates",
        action="store_true",
        help="Disable the signal-relay class filter and reproduce unconstrained path search.",
    )
    parser.add_argument(
        "--exclude-any-scaffold",
        action="store_true",
        help=(
            "Exclude every adaptor_scaffold-tagged intermediate, even when it also has "
            "a mechanistic relay class such as kinase."
        ),
    )
    parser.add_argument(
        "--disable-ontology-directionality",
        action="store_true",
        help="Retain both traversals for every edge and reproduce undirected path search.",
    )
    parser.add_argument(
        "--disable-omnipath-directionality",
        action="store_true",
        help="Use ontology rules but ignore mapped OmniPath source-target directions.",
    )
    return parser.parse_args()


def normalized_identifier(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def resolve_existing_node(
    identifier: str,
    metadata: pd.DataFrame,
    matrix_symbols: list[str],
    *,
    allow_curated_alias: bool,
) -> tuple[str, str] | None:
    """Resolve an identifier only against nodes already in the matrix."""

    key = normalized_identifier(identifier)
    if not key:
        raise ValueError("node identifier is empty")
    symbols_ci = {normalized_identifier(symbol): symbol for symbol in matrix_symbols}
    if key in symbols_ci:
        return symbols_ci[key], "exact_matrix_symbol"

    if allow_curated_alias and key in CURATED_START_ALIASES:
        symbol = CURATED_START_ALIASES[key]
        if symbol in matrix_symbols:
            return symbol, "curated_common_name_alias"

    candidates: dict[str, set[str]] = {}
    for column in ("display_symbol", "name", "stable_id"):
        if column not in metadata.columns:
            continue
        for symbol, value in zip(metadata["symbol"], metadata[column]):
            value_key = normalized_identifier(value)
            if value_key:
                candidates.setdefault(value_key, set()).add(str(symbol))
    matches = candidates.get(key, set()).intersection(matrix_symbols)
    if len(matches) == 1:
        return next(iter(matches)), "exact_node_metadata"
    if len(matches) > 1:
        raise ValueError(
            f"Identifier {identifier!r} is ambiguous among existing nodes: "
            + ", ".join(sorted(matches))
        )
    return None


def split_classes(value: object) -> set[str]:
    return {item.strip() for item in str(value).split(";") if item.strip()}


def build_intermediate_eligibility(
    metadata: pd.DataFrame,
    matrix_symbols: list[str],
    start_symbol: str,
    target_symbol: str,
    *,
    signaling_intermediates_only: bool,
    relay_classes: Iterable[str],
    exclude_any_scaffold: bool,
) -> pd.DataFrame:
    """Audit which nodes may occur between the selected endpoints."""

    relay_set = {str(value).strip() for value in relay_classes if str(value).strip()}
    metadata_by_symbol = metadata.set_index("symbol", drop=False).to_dict("index")
    rows: list[dict[str, Any]] = []
    endpoints = {start_symbol, target_symbol}
    for matrix_index, symbol in enumerate(matrix_symbols):
        record = metadata_by_symbol.get(symbol, {})
        classes = split_classes(record.get("classes", ""))
        matching = sorted(classes.intersection(relay_set))
        scaffold_tagged = "adaptor_scaffold" in classes
        if not signaling_intermediates_only:
            eligible = True
            reason = "all_intermediates_policy"
        elif matching and not (exclude_any_scaffold and scaffold_tagged):
            eligible = True
            reason = "mechanistic_relay_class:" + ";".join(matching)
        elif exclude_any_scaffold and scaffold_tagged:
            eligible = False
            reason = "excluded_adaptor_scaffold_tag"
        else:
            eligible = False
            reason = "no_mechanistic_relay_class"
        endpoint = symbol in endpoints
        permitted = eligible or endpoint
        if endpoint and not eligible:
            search_reason = "endpoint_exemption;" + reason
        elif endpoint:
            search_reason = "endpoint_and_" + reason
        else:
            search_reason = reason
        rows.append(
            {
                "matrix_index": matrix_index,
                "symbol": symbol,
                "display_symbol": record.get("display_symbol", symbol),
                "name": record.get("name", ""),
                "classes": record.get("classes", ""),
                "node_type": record.get("node_type", ""),
                "matching_relay_classes": ";".join(matching),
                "adaptor_scaffold_tagged": scaffold_tagged,
                "allowed_as_intermediate": eligible,
                "is_selected_endpoint": endpoint,
                "permitted_in_search": permitted,
                "eligibility_reason": search_reason,
            }
        )
    return pd.DataFrame(rows)


def validate_matrix(matrix: pd.DataFrame) -> None:
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("adjacency matrix is not square")
    if matrix.index.astype(str).tolist() != matrix.columns.astype(str).tolist():
        raise ValueError("adjacency matrix row and column labels differ")
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("adjacency matrix contains non-finite values")
    if ((values < 0.0) | (values > 1.0)).any():
        raise ValueError("adjacency matrix contains values outside [0, 1]")
    if not np.array_equal(values, values.T):
        raise ValueError("adjacency matrix is not exactly symmetric")
    if not np.array_equal(np.diag(values), np.zeros(len(values))):
        raise ValueError("adjacency matrix diagonal is not zero")


def cached_target_extension(
    project: Path,
    target: str,
    target_symbol: str,
) -> tuple[pd.DataFrame, Path, dict[str, Any]] | None:
    """Load a cached extension only when its inputs and output hash are current."""

    output = project / "results/path_finding/target_extensions" / safe_name(target_symbol)
    matrix_path = output / "extended_adjacency_matrix.tsv"
    summary_path = output / "analysis_summary.json"
    if not matrix_path.exists() or not summary_path.exists():
        return None
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        expected_inputs = summary["input_sha256"]
        current_inputs = {
            "node_universe_combined_nonzero.tsv": sha256_file(project / UNIVERSE_RELATIVE),
            "combined_adjacency_matrix.tsv": sha256_file(project / CURRENT_MATRIX_RELATIVE),
            "uniprot_mouse_reference_proteome.tsv.gz": sha256_file(project / UNIPROT_RELATIVE),
        }
        if summary.get("target_symbol", "").casefold() != target_symbol.casefold():
            return None
        if expected_inputs != current_inputs:
            return None
        expected_matrix_hash = summary["output_sha256"]["extended_adjacency_matrix.tsv"]
        if sha256_file(matrix_path) != expected_matrix_hash:
            return None
        if not all(summary.get("validation", {}).values()):
            return None
        matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
        validate_matrix(matrix)
        return matrix, matrix_path, summary
    except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError):
        return None


def load_or_extend_graph(
    project: Path,
    target: str,
    *,
    auto_extend_target: bool,
    reuse_target_extension: bool,
) -> tuple[pd.DataFrame, Path, str, str, dict[str, Any] | None]:
    """Load the active graph and add an external target when necessary."""

    universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    current_path = project / CURRENT_MATRIX_RELATIVE
    matrix = pd.read_csv(current_path, sep="\t", index_col=0)
    validate_matrix(matrix)
    symbols = matrix.index.astype(str).tolist()
    existing = resolve_existing_node(
        target,
        universe,
        symbols,
        allow_curated_alias=False,
    )
    if existing is not None:
        symbol, method = existing
        return matrix, current_path, symbol, method, None
    if not auto_extend_target:
        raise ValueError(
            f"Target {target!r} is not in the active graph and automatic extension is disabled"
        )

    target_symbol, _, target_mapping_method = resolve_target(
        target,
        project / UNIPROT_RELATIVE,
    )
    if reuse_target_extension:
        cached = cached_target_extension(project, target, target_symbol)
        if cached is not None:
            extended, matrix_path, summary = cached
            return (
                extended,
                matrix_path,
                target_symbol,
                f"external_target_cached:{target_mapping_method}",
                summary,
            )

    result: TargetAdjacencyResult = build_target_adjacency_vector(
        target,
        project_root=project,
        write_outputs=True,
        write_extended_matrix=True,
    )
    if result.extended_matrix is None or result.output_directory is None:
        raise RuntimeError("target extension did not produce an extended matrix")
    matrix_path = result.output_directory / "extended_adjacency_matrix.tsv"
    return (
        result.extended_matrix,
        matrix_path,
        result.target_symbol,
        f"external_target_recomputed:{target_mapping_method}",
        result.summary,
    )


def build_adjacency(
    matrix: pd.DataFrame,
    minimum_edge_probability: float,
    *,
    directed: bool = False,
) -> tuple[list[list[tuple[int, float, float]]], int]:
    """Create an adjacency list of traversals strictly above the chosen cutoff.

    Symmetric matrices use one unordered edge count by default. For a partially
    directed propagation matrix, ``directed=True`` scans every ordered pair and
    returns the number of retained directed traversals.
    """

    values = matrix.to_numpy(dtype=float)
    size = len(values)
    adjacency: list[list[tuple[int, float, float]]] = [[] for _ in range(size)]
    edge_count = 0
    if directed:
        for i in range(size):
            for j in range(size):
                if i == j:
                    continue
                probability = float(values[i, j])
                if probability <= minimum_edge_probability:
                    continue
                adjacency[i].append((j, probability, -math.log(probability)))
                edge_count += 1
    else:
        for i in range(size):
            for j in range(i + 1, size):
                probability = float(values[i, j])
                if probability <= minimum_edge_probability:
                    continue
                cost = -math.log(probability)
                adjacency[i].append((j, probability, cost))
                adjacency[j].append((i, probability, cost))
                edge_count += 1
    for neighbors in adjacency:
        neighbors.sort(key=lambda item: (item[2], item[0]))
    return adjacency, edge_count


def shortest_hop_limited_path(
    adjacency: list[list[tuple[int, float, float]]],
    start: int,
    target: int,
    max_hops: int,
    *,
    banned_nodes: frozenset[int] = frozenset(),
    banned_transitions: frozenset[tuple[int, int]] = frozenset(),
    permitted_nodes: frozenset[int] | None = None,
) -> SearchPath | None:
    """Dijkstra search on a hop-layered graph with optional exclusions."""

    if (
        max_hops < 0
        or start in banned_nodes
        or target in banned_nodes
        or (permitted_nodes is not None and start not in permitted_nodes)
        or (permitted_nodes is not None and target not in permitted_nodes)
    ):
        return None
    if start == target:
        return SearchPath(0.0, (start,))

    counter = itertools.count()
    heap: list[tuple[float, int, int, int]] = [(0.0, 0, next(counter), start)]
    distances: dict[tuple[int, int], float] = {(start, 0): 0.0}
    predecessor: dict[tuple[int, int], tuple[int, int]] = {}

    while heap:
        cost, hops, _, node = heapq.heappop(heap)
        state = (node, hops)
        if cost > distances.get(state, math.inf) + 1e-15:
            continue
        if node == target:
            states = [state]
            while states[-1] in predecessor:
                states.append(predecessor[states[-1]])
            nodes = tuple(item[0] for item in reversed(states))
            # Positive edge costs make an optimal path loopless. Retain a
            # defensive check for inputs containing exact probability-one cycles.
            if len(nodes) == len(set(nodes)):
                return SearchPath(cost, nodes)
        if hops >= max_hops:
            continue
        next_hops = hops + 1
        for neighbor, _, edge_cost in adjacency[node]:
            if (
                neighbor in banned_nodes
                or (node, neighbor) in banned_transitions
                or (permitted_nodes is not None and neighbor not in permitted_nodes)
            ):
                continue
            next_state = (neighbor, next_hops)
            next_cost = cost + edge_cost
            if next_cost + 1e-15 < distances.get(next_state, math.inf):
                distances[next_state] = next_cost
                predecessor[next_state] = state
                heapq.heappush(
                    heap,
                    (next_cost, next_hops, next(counter), neighbor),
                )
    return None


def path_cost(nodes: tuple[int, ...], probabilities: np.ndarray) -> float:
    return float(
        sum(-math.log(float(probabilities[a, b])) for a, b in zip(nodes, nodes[1:]))
    )


def k_shortest_simple_paths(
    adjacency: list[list[tuple[int, float, float]]],
    probabilities: np.ndarray,
    start: int,
    target: int,
    *,
    top_k: int,
    max_hops: int,
    permitted_nodes: frozenset[int] | None = None,
) -> list[SearchPath]:
    """Return hop-limited loopless paths using Yen's ranking algorithm."""

    first = shortest_hop_limited_path(
        adjacency,
        start,
        target,
        max_hops,
        permitted_nodes=permitted_nodes,
    )
    if first is None:
        return []
    accepted = [first]
    accepted_nodes = {first.nodes}
    candidate_heap: list[tuple[float, tuple[int, ...]]] = []
    candidate_nodes: set[tuple[int, ...]] = set()

    while len(accepted) < top_k:
        previous = accepted[-1].nodes
        for spur_index in range(len(previous) - 1):
            root = previous[: spur_index + 1]
            banned_transitions: set[tuple[int, int]] = set()
            for accepted_path in accepted:
                nodes = accepted_path.nodes
                if len(nodes) > spur_index and nodes[: spur_index + 1] == root:
                    banned_transitions.add((nodes[spur_index], nodes[spur_index + 1]))

            remaining_hops = max_hops - spur_index
            spur = shortest_hop_limited_path(
                adjacency,
                root[-1],
                target,
                remaining_hops,
                banned_nodes=frozenset(root[:-1]),
                banned_transitions=frozenset(banned_transitions),
                permitted_nodes=permitted_nodes,
            )
            if spur is None:
                continue
            total_nodes = root[:-1] + spur.nodes
            if (
                len(total_nodes) != len(set(total_nodes))
                or len(total_nodes) - 1 > max_hops
                or total_nodes in accepted_nodes
                or total_nodes in candidate_nodes
            ):
                continue
            total_cost = path_cost(total_nodes, probabilities)
            heapq.heappush(candidate_heap, (total_cost, total_nodes))
            candidate_nodes.add(total_nodes)

        while candidate_heap:
            cost, nodes = heapq.heappop(candidate_heap)
            candidate_nodes.discard(nodes)
            if nodes not in accepted_nodes:
                accepted.append(SearchPath(cost, nodes))
                accepted_nodes.add(nodes)
                break
        else:
            break
    return accepted


def connected_component_size(
    adjacency: list[list[tuple[int, float, float]]],
    start: int,
    permitted_nodes: frozenset[int] | None = None,
) -> int:
    visited = {start}
    stack = [start]
    while stack:
        node = stack.pop()
        for neighbor, _, _ in adjacency[node]:
            if (
                neighbor not in visited
                and (permitted_nodes is None or neighbor in permitted_nodes)
            ):
                visited.add(neighbor)
                stack.append(neighbor)
    return len(visited)


def serialize_paths(
    ranked: Iterable[SearchPath],
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    symbols = matrix.index.astype(str).tolist()
    values = matrix.to_numpy(dtype=float)
    metadata_by_symbol = metadata.set_index("symbol", drop=False).to_dict("index")
    path_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    for rank, found in enumerate(ranked, start=1):
        node_symbols = [symbols[index] for index in found.nodes]
        edge_probabilities = [
            float(values[a, b]) for a, b in zip(found.nodes, found.nodes[1:])
        ]
        probability_product = math.exp(-found.cost)
        hops = len(found.nodes) - 1
        path_rows.append(
            {
                "rank": rank,
                "start_symbol": node_symbols[0],
                "target_symbol": node_symbols[-1],
                "hop_count": hops,
                "node_count": len(node_symbols),
                "path_probability_product": probability_product,
                "negative_log_path_probability": found.cost,
                "geometric_mean_edge_probability": (
                    probability_product ** (1.0 / hops) if hops else 1.0
                ),
                "minimum_edge_probability": min(edge_probabilities) if edge_probabilities else 1.0,
                "maximum_edge_probability": max(edge_probabilities) if edge_probabilities else 1.0,
                "path_symbols": " -> ".join(node_symbols),
                "intermediate_symbols": ";".join(node_symbols[1:-1]),
                "intermediate_classes": " | ".join(
                    str(metadata_by_symbol.get(symbol, {}).get("classes", ""))
                    for symbol in node_symbols[1:-1]
                ),
                "edge_probabilities": ";".join(f"{value:.9g}" for value in edge_probabilities),
            }
        )
        for step, (source_index, target_index, probability) in enumerate(
            zip(found.nodes, found.nodes[1:], edge_probabilities),
            start=1,
        ):
            source = symbols[source_index]
            target = symbols[target_index]
            edge_rows.append(
                {
                    "path_rank": rank,
                    "step": step,
                    "source_symbol": source,
                    "target_symbol": target,
                    "edge_probability": probability,
                    "negative_log_edge_probability": -math.log(probability),
                    "source_name": metadata_by_symbol.get(source, {}).get("name", ""),
                    "target_name": metadata_by_symbol.get(target, {}).get("name", ""),
                    "source_classes": metadata_by_symbol.get(source, {}).get("classes", ""),
                    "target_classes": metadata_by_symbol.get(target, {}).get("classes", ""),
                }
            )
    return pd.DataFrame(path_rows), pd.DataFrame(edge_rows)


def find_ranked_paths(
    start: str,
    target: str,
    *,
    project_root: Path | str = PROJECT_DEFAULT,
    top_k: int = 25,
    max_hops: int = 6,
    minimum_edge_probability: float = 0.5,
    signaling_intermediates_only: bool = True,
    relay_classes: Iterable[str] = DEFAULT_SIGNAL_RELAY_CLASSES,
    exclude_any_scaffold: bool = False,
    ontology_directionality_enabled: bool = True,
    omnipath_directionality_enabled: bool = True,
    auto_extend_target: bool = True,
    reuse_target_extension: bool = True,
    output_dir: Path | str | None = None,
    write_outputs: bool = True,
) -> RankedPathResult:
    """Find the top simple paths between a current node and chosen target.

    Paths are ranked by descending product of edge probabilities. The equivalent
    additive search cost is the sum of ``-log(p_edge)``. By default, exact
    neutral-baseline edges are excluded, no path may exceed six hops, and every
    intermediate must have a mechanistic signal-relay class.
    """

    if top_k < 1:
        raise ValueError("top_k must be at least 1")
    if max_hops < 1:
        raise ValueError("max_hops must be at least 1")
    if not 0.0 <= minimum_edge_probability < 1.0:
        raise ValueError("minimum_edge_probability must be in [0, 1)")
    relay_classes = tuple(
        str(value).strip() for value in relay_classes if str(value).strip()
    )
    if signaling_intermediates_only and not relay_classes:
        raise ValueError("at least one relay class is required for constrained search")

    project = Path(project_root).resolve()
    universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    matrix, matrix_path, target_symbol, target_method, extension_summary = load_or_extend_graph(
        project,
        target,
        auto_extend_target=auto_extend_target,
        reuse_target_extension=reuse_target_extension,
    )
    symbols = matrix.index.astype(str).tolist()
    start_result = resolve_existing_node(
        start,
        universe,
        symbols,
        allow_curated_alias=True,
    )
    if start_result is None:
        raise ValueError(
            f"Starting node {start!r} could not be resolved in the active graph. "
            "Use its mouse gene symbol."
        )
    start_symbol, start_method = start_result
    if start_symbol == target_symbol:
        raise ValueError("start and target resolve to the same node")

    target_metadata = {
        "symbol": target_symbol,
        "display_symbol": target_symbol,
        "name": (
            extension_summary.get("target_protein_name", "")
            if extension_summary is not None
            else ""
        ),
        "classes": "external_target" if target_symbol not in set(universe["symbol"]) else "",
        "node_type": "protein" if extension_summary is not None else "",
    }
    metadata = universe.copy()
    if target_symbol not in set(metadata["symbol"]):
        metadata = pd.concat([metadata, pd.DataFrame([target_metadata])], ignore_index=True)

    propagation_matrix = matrix
    directionality_audit = pd.DataFrame()
    directionality_class_catalog = pd.DataFrame()
    directionality_summary: dict[str, Any] | None = None
    directionality_catalog: dict[str, Any] | None = None
    omnipath_direction_evidence = pd.DataFrame()
    if ontology_directionality_enabled:
        directionality_catalog = load_direction_rule_catalog(RULE_CATALOG_PATH)
        directionality_class_catalog = build_complete_class_pair_catalog(
            directionality_catalog,
            directionality_catalog["ontology_classes"],
        )
        if omnipath_directionality_enabled:
            omnipath_direction_path = project / OMNIPATH_DIRECTION_RAW_RELATIVE
            if not omnipath_direction_path.is_file():
                raise FileNotFoundError(
                    "OmniPath directionality is enabled but its raw directed "
                    f"interaction table is missing: {omnipath_direction_path}"
                )
            omnipath_direction_raw = pd.read_csv(
                omnipath_direction_path,
                sep="\t",
                dtype=str,
                keep_default_na=False,
            )
            omnipath_direction_evidence = build_omnipath_direction_evidence(
                omnipath_direction_raw,
                symbols,
            )
        (
            propagation_matrix,
            directionality_audit,
            directionality_summary,
        ) = apply_ontology_directionality(
            matrix,
            metadata,
            directionality_catalog,
            audit_probability_cutoff=minimum_edge_probability,
            edge_output_cutoff=minimum_edge_probability,
            path_probability_cutoff=minimum_edge_probability,
            omnipath_directions=(
                omnipath_direction_evidence
                if omnipath_directionality_enabled
                else None
            ),
        )
        directionality_summary["omnipath_direction_source"] = (
            str(OMNIPATH_DIRECTION_RAW_RELATIVE)
            if omnipath_directionality_enabled
            else None
        )

    eligibility = build_intermediate_eligibility(
        metadata,
        symbols,
        start_symbol,
        target_symbol,
        signaling_intermediates_only=signaling_intermediates_only,
        relay_classes=relay_classes,
        exclude_any_scaffold=exclude_any_scaffold,
    )
    permitted_nodes = frozenset(
        eligibility.loc[eligibility["permitted_in_search"], "matrix_index"].astype(int)
    )
    eligible_intermediates = frozenset(
        eligibility.loc[eligibility["allowed_as_intermediate"], "matrix_index"].astype(int)
    )

    adjacency, retained_edge_count = build_adjacency(
        propagation_matrix,
        minimum_edge_probability,
        directed=ontology_directionality_enabled,
    )
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}
    start_index = symbol_to_index[start_symbol]
    target_index = symbol_to_index[target_symbol]
    ranked = k_shortest_simple_paths(
        adjacency,
        propagation_matrix.to_numpy(dtype=float),
        start_index,
        target_index,
        top_k=top_k,
        max_hops=max_hops,
        permitted_nodes=permitted_nodes,
    )

    paths, path_edges = serialize_paths(ranked, propagation_matrix, metadata)

    probabilities = propagation_matrix.to_numpy(dtype=float)
    retained_search_edge_count = sum(
        1
        for source, neighbors in enumerate(adjacency)
        if source in permitted_nodes
        for target_node, _, _ in neighbors
        if target_node in permitted_nodes
        and (ontology_directionality_enabled or source < target_node)
    )
    validation = {
        "matrix_validated_before_search": True,
        "all_paths_start_at_requested_node": bool(
            paths.empty or (paths["start_symbol"] == start_symbol).all()
        ),
        "all_paths_end_at_requested_target": bool(
            paths.empty or (paths["target_symbol"] == target_symbol).all()
        ),
        "all_paths_are_simple": all(
            len(path.nodes) == len(set(path.nodes)) for path in ranked
        ),
        "all_paths_respect_hop_limit": all(
            len(path.nodes) - 1 <= max_hops for path in ranked
        ),
        "all_internal_nodes_satisfy_relay_policy": all(
            all(node in eligible_intermediates for node in path.nodes[1:-1])
            for path in ranked
        ),
        "all_edges_exceed_cutoff": bool(
            path_edges.empty
            or (path_edges["edge_probability"] > minimum_edge_probability).all()
        ),
        "ranking_is_nonincreasing_by_probability_product": bool(
            paths.empty
            or np.all(np.diff(paths["path_probability_product"].to_numpy(float)) <= 1e-15)
        ),
        "reported_products_reconstruct_from_edges": all(
            math.isclose(
                math.exp(-path.cost),
                math.prod(
                    float(probabilities[a, b])
                    for a, b in zip(path.nodes, path.nodes[1:])
                ),
                rel_tol=1e-12,
                abs_tol=1e-15,
            )
            for path in ranked
        ),
    }
    if not all(validation.values()):
        raise AssertionError(f"path validation failed: {validation}")

    upper_triangle = np.triu_indices(len(propagation_matrix), k=1)
    retained_unique_edge_count = int(
        (
            np.maximum(probabilities, probabilities.T)[upper_triangle]
            > minimum_edge_probability
        ).sum()
    )
    start_reachable_node_count = connected_component_size(
        adjacency, start_index, permitted_nodes
    )

    summary: dict[str, Any] = {
        "start_input": start,
        "start_symbol": start_symbol,
        "start_resolution_method": start_method,
        "target_input": target,
        "target_symbol": target_symbol,
        "target_resolution_method": target_method,
        "target_was_external": extension_summary is not None,
        "matrix_path": str(matrix_path),
        "matrix_node_count": len(matrix),
        "retained_unique_edge_count_before_intermediate_filter": (
            retained_unique_edge_count
        ),
        "retained_transition_count_before_intermediate_filter": retained_edge_count,
        "retained_transition_count_in_search_subgraph": retained_search_edge_count,
        "ontology_directionality_enabled": ontology_directionality_enabled,
        "omnipath_directionality_enabled": bool(
            ontology_directionality_enabled and omnipath_directionality_enabled
        ),
        "ontology_directionality": directionality_summary,
        "minimum_edge_probability_exclusive": minimum_edge_probability,
        "neutral_baseline_edges_excluded": minimum_edge_probability >= 0.5,
        "signaling_intermediates_only": signaling_intermediates_only,
        "signal_relay_classes": sorted(
            {str(value).strip() for value in relay_classes if str(value).strip()}
        ),
        "exclude_any_adaptor_scaffold_tag": exclude_any_scaffold,
        "eligible_intermediate_node_count": int(
            eligibility["allowed_as_intermediate"].sum()
        ),
        "excluded_intermediate_node_count": int(
            (~eligibility["allowed_as_intermediate"]).sum()
        ),
        "endpoint_exemptions": eligibility.loc[
            eligibility["is_selected_endpoint"]
            & ~eligibility["allowed_as_intermediate"],
            "symbol",
        ].tolist(),
        "intermediate_policy": (
            "Sensitivity-oriented policy: internal nodes require at least one "
            "signaling class other than adaptor_scaffold. Ligand, generic "
            "signaling-process, kinase/phosphatase-binding, second-messenger-binding, "
            "and catalytic signaling classes all qualify. The selected endpoints are "
            "always permitted. "
            + (
                "Every adaptor_scaffold-tagged internal node is excluded, including "
                "multi-role proteins with another relay class."
                if exclude_any_scaffold
                else "A multi-role scaffold remains eligible through another relay class."
            )
            if signaling_intermediates_only
            else "All matrix nodes may be internal nodes."
        ),
        "maximum_hops": max_hops,
        "top_k_requested": top_k,
        "paths_found": len(paths),
        "ranking_rule": (
            "Descending product of edge probabilities; search minimizes the "
            "equivalent sum of -log(edge_probability)."
        ),
        "path_constraint": (
            "Simple paths following allowed ontology/OmniPath-directed traversals; no "
            "node may repeat within a path. Unresolved edges retain both traversals."
            if ontology_directionality_enabled
            else "Simple undirected paths only; no node may repeat within a path."
        ),
        "start_reachable_node_count_after_filtering": start_reachable_node_count,
        "start_component_size_after_filtering": (
            None
            if ontology_directionality_enabled
            else start_reachable_node_count
        ),
        "target_degree_after_filtering": sum(
            neighbor in permitted_nodes for neighbor, _, _ in adjacency[target_index]
        ),
        "target_in_degree_after_filtering": sum(
            source in permitted_nodes
            for source, neighbors in enumerate(adjacency)
            if any(neighbor == target_index for neighbor, _, _ in neighbors)
        ),
        "direct_edge_retained": any(
            neighbor == target_index for neighbor, _, _ in adjacency[start_index]
        ),
        "validation": validation,
        "input_sha256": {
            matrix_path.name: sha256_file(matrix_path),
            (project / UNIVERSE_RELATIVE).name: sha256_file(project / UNIVERSE_RELATIVE),
            **(
                {
                    (project / OMNIPATH_DIRECTION_RAW_RELATIVE).name: sha256_file(
                        project / OMNIPATH_DIRECTION_RAW_RELATIVE
                    )
                }
                if ontology_directionality_enabled and omnipath_directionality_enabled
                else {}
            ),
        },
    }

    resolved_output: Path | None = None
    if write_outputs:
        resolved_output = (
            Path(output_dir).resolve()
            if output_dir is not None
            else project
            / "results/path_finding/ranked_paths"
            / (
                (
                    f"{safe_name(start_symbol)}_to_{safe_name(target_symbol)}_"
                    + ("ontology_directed_" if ontology_directionality_enabled else "")
                    + (
                        "signal_relay_no_scaffold"
                        if exclude_any_scaffold
                        else "signal_relay"
                    )
                )
                if signaling_intermediates_only
                else f"{safe_name(start_symbol)}_to_{safe_name(target_symbol)}_all_nodes"
            )
        )
        resolved_output.mkdir(parents=True, exist_ok=True)
        paths_path = resolved_output / "ranked_paths.tsv"
        edges_path = resolved_output / "ranked_path_edges.tsv"
        eligibility_path = resolved_output / "intermediate_node_eligibility.tsv"
        summary_path = resolved_output / "analysis_summary.json"
        readme_path = resolved_output / "README.md"
        propagation_path = resolved_output / "propagation_adjacency_matrix.tsv"
        directionality_audit_path = resolved_output / "ontology_directionality_audit.tsv.gz"
        class_catalog_path = resolved_output / "ontology_class_pair_catalog.tsv"
        rules_path = resolved_output / "ontology_direction_rules.json"
        omnipath_direction_path = resolved_output / "omnipath_direction_evidence.tsv.gz"
        paths.to_csv(paths_path, sep="\t", index=False, float_format="%.12g")
        path_edges.to_csv(edges_path, sep="\t", index=False, float_format="%.12g")
        eligibility.to_csv(eligibility_path, sep="\t", index=False)
        if ontology_directionality_enabled and directionality_catalog is not None:
            propagation_matrix.to_csv(propagation_path, sep="\t", float_format="%.9g")
            directionality_audit.to_csv(
                directionality_audit_path,
                sep="\t",
                index=False,
                compression="gzip",
            )
            directionality_class_catalog.to_csv(class_catalog_path, sep="\t", index=False)
            rules_path.write_text(
                json.dumps(
                    {
                        key: value
                        for key, value in directionality_catalog.items()
                        if key != "catalog_path"
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            if omnipath_directionality_enabled:
                omnipath_direction_evidence.to_csv(
                    omnipath_direction_path,
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
        policy_readme = (
            "Internal nodes are permitted through this sensitivity-oriented signaling "
            f"class set: `{';'.join(summary['signal_relay_classes'])}`. Ligand, binding, "
            "and generic signaling-process labels qualify. Endpoints remain allowed even "
            "when they do not qualify as intermediates. "
            + (
                "Every adaptor_scaffold-tagged internal node was excluded."
                if exclude_any_scaffold
                else "Multi-role scaffolds were permitted through another relay class."
            )
            if signaling_intermediates_only
            else "All matrix nodes were permitted as intermediates for this run."
        )
        readme_path.write_text(
            f"""# Ranked paths: {start_symbol} to {target_symbol}

Paths are simple, contain at most {max_hops} edges, and use only
edges with probability strictly greater than {minimum_edge_probability}. They
are ranked by the product of their edge probabilities. This is equivalent to
minimizing the sum of `-log(edge_probability)`.

{"Ontology rules are applied first; mapped OmniPath source-target records add unique directions only for pairs ontology left unresolved. OmniPath cannot reopen or reverse an ontology-disallowed traversal. Remaining unresolved pairs retain both traversals." if ontology_directionality_enabled and omnipath_directionality_enabled else "Ontology rules partially orient the propagation graph. Uniquely disallowed reverse traversals are removed; unresolved pairs retain both traversals." if ontology_directionality_enabled else "Directionality was disabled, so every retained edge can be traversed both ways."}

`ranked_paths.tsv` contains one row per complete path. `ranked_path_edges.tsv`
contains one row per edge per path, with node names and signaling classes.
`intermediate_node_eligibility.tsv` records the class-based decision for every
matrix node. {policy_readme}

The product score is a model-based ranking quantity. It should not be reported
as a calibrated probability that the entire biological route is true because
the edge evidence streams and adjacent edges are not necessarily independent.
""",
            encoding="utf-8",
        )
        summary["output_sha256"] = {
            paths_path.name: sha256_file(paths_path),
            edges_path.name: sha256_file(edges_path),
            eligibility_path.name: sha256_file(eligibility_path),
            readme_path.name: sha256_file(readme_path),
            **(
                {
                    propagation_path.name: sha256_file(propagation_path),
                    directionality_audit_path.name: sha256_file(directionality_audit_path),
                    class_catalog_path.name: sha256_file(class_catalog_path),
                    rules_path.name: sha256_file(rules_path),
                    **(
                        {
                            omnipath_direction_path.name: sha256_file(
                                omnipath_direction_path
                            )
                        }
                        if omnipath_directionality_enabled
                        else {}
                    ),
                }
                if ontology_directionality_enabled
                else {}
            ),
        }
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    return RankedPathResult(
        start_symbol=start_symbol,
        target_symbol=target_symbol,
        paths=paths,
        path_edges=path_edges,
        summary=summary,
        output_directory=resolved_output,
        matrix_path=matrix_path,
    )


def main() -> int:
    args = parse_args()
    result = find_ranked_paths(
        args.start,
        args.target,
        project_root=args.project_root,
        top_k=args.top_k,
        max_hops=args.max_hops,
        minimum_edge_probability=args.minimum_edge_probability,
        signaling_intermediates_only=not args.allow_all_intermediates,
        exclude_any_scaffold=args.exclude_any_scaffold,
        ontology_directionality_enabled=not args.disable_ontology_directionality,
        omnipath_directionality_enabled=not args.disable_omnipath_directionality,
        auto_extend_target=not args.no_auto_extend,
        reuse_target_extension=not args.rebuild_target_extension,
        output_dir=args.output_dir,
        write_outputs=True,
    )
    print(json.dumps(result.summary, indent=2))
    if not result.paths.empty:
        print("\nTop paths:")
        print(
            result.paths[
                ["rank", "hop_count", "path_probability_product", "path_symbols"]
            ]
            .head(10)
            .to_string(index=False)
        )
    if result.output_directory is not None:
        print(f"\nOutput directory: {result.output_directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
