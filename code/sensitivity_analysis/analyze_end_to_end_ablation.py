#!/usr/bin/env python3
"""Run Tq-aware, end-to-end leave-one-evidence-stream-out analysis.

By default, every registered node and edge evidence stream is tested, even
when it is disabled in the generic GUI profile.  A configured stream is
removed from the supplied model.  An optional stream is first enabled on top
of that model, then removed again in its matched comparator.  Node-stream
ablation is genuinely end to end: node selection is repeated, the graph is
rebuilt on the resulting nodes, ontology directionality is reapplied, and
paths are reranked.  Edge-stream ablation keeps node selection fixed but
rebuilds the edge graph and paths.

Mutually exclusive alternatives cannot coexist in one valid model.  An
optional member therefore replaces the configured member of its exclusive
group in its matched reference context.

The comparison is made both at the configured Tq/reference multiplier and
across a bounded stream-specific multiplier grid.  Only the focal stream's
multiplier changes; all other Tq values, weights, priors, cutoffs, and path
settings remain fixed.  The ablated model is therefore constant across that
stream's sweep and is evaluated once.
"""

from __future__ import annotations

import argparse
import copy
import gc
import json
import math
import sys
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
GUI_DIR = PROJECT_ROOT / "gui"
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

from workflow_engine import (  # noqa: E402
    RULE_CATALOG_PATH,
    UNIVERSE_RELATIVE,
    _calibration_tq_bounds,
    append_requested_graph_endpoint,
    apply_ontology_directionality,
    combine_edge_factors,
    default_configuration,
    ensure_incremental_pairs,
    load_direction_rule_catalog,
    load_registry,
    normalize_configuration,
    run_paths,
    scaffold_triadic_closure_factors,
    select_nodes,
    stable_expit,
    supported_edge_table,
)


PROBABILITY_EPSILON = 1e-15
REFERENCE_TOLERANCE = 1e-12


@dataclass
class Snapshot:
    """In-memory output required for one model comparison."""

    config: dict[str, Any]
    node_factors: pd.DataFrame
    selected_nodes: pd.DataFrame
    matrix: pd.DataFrame
    propagation_matrix: pd.DataFrame
    supported_edges: pd.DataFrame
    paths: pd.DataFrame
    path_status: str
    target_symbol: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Project root containing gui/evidence_registry.json.",
    )
    parser.add_argument(
        "--configuration",
        type=Path,
        default=None,
        help=(
            "Optional GUI configuration JSON. If omitted, the registered default "
            "configuration is analyzed. Calibration is not refit during ablation."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: dated results/sensitivity_analysis folder).",
    )
    parser.add_argument(
        "--stream-scope",
        choices=("all", "configured"),
        default="all",
        help=(
            "Streams to ablate. 'all' (default) tests every registered stream; "
            "streams that are currently off are evaluated in matched add-one "
            "contexts, with exclusive conflicts replaced. "
            "'configured' reproduces the earlier behavior and analyzes only "
            "streams enabled in the supplied configuration."
        ),
    )
    parser.add_argument(
        "--resume-compatible-reference",
        type=Path,
        default=None,
        help=(
            "Optional earlier completed ablation directory. Its analyzed "
            "configuration must exactly match the supplied baseline; matching "
            "per-stream results are imported as checkpoints before the run."
        ),
    )
    parser.add_argument(
        "--optional-tq-grid",
        choices=("reference", "compact"),
        default="reference",
        help=(
            "End-to-end Tq grid for streams disabled in the supplied profile. "
            "'reference' (default) evaluates their configured Tq only; use the "
            "separate all-stream Tq-response analysis for dense BF/information "
            "curves. 'compact' adds registered low/high bounds and a distinct "
            "preferred value, but may create very large graphs."
        ),
    )
    return parser.parse_args()


def _exclusive_groups(registry: dict[str, Any]) -> dict[str, list[str]]:
    groups: dict[str, list[str]] = {}
    for definition in registry["edge_streams"]:
        group = definition.get("exclusive_group")
        if group:
            groups.setdefault(str(group), []).append(str(definition["id"]))
    return groups


def focal_context_configuration(
    baseline_config: dict[str, Any],
    registry: dict[str, Any],
    *,
    stage: str,
    stream_id: str,
) -> tuple[dict[str, Any], str]:
    """Return the supplied or matched add-one context containing the focal stream."""
    config = copy.deepcopy(baseline_config)
    definitions = {
        item["id"]: item for item in registry[f"{stage}_streams"]
    }
    focal_state = config[f"{stage}_streams"][stream_id]
    already_active = bool(focal_state["enabled"] and focal_state["weight"] > 0)
    if float(focal_state["weight"]) <= 0:
        focal_state["weight"] = float(definitions[stream_id].get("default_weight", 1.0))
    focal_state["enabled"] = True
    if stage == "edge":
        group = definitions[stream_id].get("exclusive_group")
        if group:
            members = _exclusive_groups(registry)[str(group)]
            for member in members:
                config["edge_streams"][member]["enabled"] = member == stream_id
    if already_active:
        return config, "supplied_configuration"
    return config, f"supplied_plus[{stage}:{stream_id}]"


def bernoulli_kl_bits(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Return D_KL(Bernoulli(p) || Bernoulli(q)) in bits."""
    p_safe = np.clip(np.asarray(p, dtype=float), PROBABILITY_EPSILON, 1 - PROBABILITY_EPSILON)
    q_safe = np.clip(np.asarray(q, dtype=float), PROBABILITY_EPSILON, 1 - PROBABILITY_EPSILON)
    return p_safe * np.log2(p_safe / q_safe) + (1 - p_safe) * np.log2(
        (1 - p_safe) / (1 - q_safe)
    )


def safe_spearman(left: np.ndarray, right: np.ndarray) -> float:
    finite = np.isfinite(left) & np.isfinite(right)
    if int(finite.sum()) < 3:
        return math.nan
    left_values = left[finite]
    right_values = right[finite]
    if np.ptp(left_values) <= REFERENCE_TOLERANCE or np.ptp(right_values) <= REFERENCE_TOLERANCE:
        return math.nan
    left_rank = pd.Series(left_values).rank(method="average").to_numpy(float)
    right_rank = pd.Series(right_values).rank(method="average").to_numpy(float)
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((str(left), str(right)), key=lambda item: (item.casefold(), item)))


def edge_probability_map(matrix: pd.DataFrame, symbols: list[str]) -> dict[tuple[str, str], float]:
    if len(symbols) < 2:
        return {}
    aligned = matrix.reindex(index=symbols, columns=symbols)
    values = aligned.to_numpy(float)
    left, right = np.triu_indices(len(symbols), 1)
    return {
        canonical_pair(symbols[i], symbols[j]): float(values[i, j])
        for i, j in zip(left, right)
    }


def supported_pair_set(edges: pd.DataFrame) -> set[tuple[str, str]]:
    if edges.empty:
        return set()
    return {
        canonical_pair(left, right)
        for left, right in edges[["node_a", "node_b"]].itertuples(index=False, name=None)
    }


def path_map(paths: pd.DataFrame) -> dict[str, dict[str, float | int]]:
    if paths.empty:
        return {}
    return {
        str(row.path_symbols): {
            "rank": int(row.rank),
            "score": float(row.path_probability_product),
        }
        for row in paths[["path_symbols", "rank", "path_probability_product"]].itertuples(index=False)
    }


def tq_grid(
    definition: dict[str, Any],
    reference: float,
    *,
    compact: bool = False,
) -> list[float]:
    """Return a compact grid spanning the registered calibration bounds."""
    lower, upper = _calibration_tq_bounds(definition)
    if compact:
        preferred = definition.get("normalization", {}).get(
            "calibration_preferred_multiplier"
        )
        candidates = [lower, reference, upper]
        if preferred is not None and lower <= float(preferred) <= upper:
            candidates.append(float(preferred))
    elif upper <= 1.0 + REFERENCE_TOLERANCE:
        candidates = [lower, 0.1, 0.25, 0.5, 1.0, upper]
    else:
        candidates = [lower, 0.5, 1.0, 2.0, upper]
    candidates.append(float(reference))
    return sorted(
        {
            round(float(value), 12)
            for value in candidates
            if lower - REFERENCE_TOLERANCE <= float(value) <= upper + REFERENCE_TOLERANCE
        }
    )


class Evaluator:
    """Evaluate workflow configurations without creating a GUI run directory."""

    def __init__(self, project: Path, registry: dict[str, Any]):
        self.project = project
        self.registry = registry
        self.universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
        self.seed_symbols = self.universe["symbol"].astype(str).tolist()
        self.seed_set = set(self.seed_symbols)
        self.direction_catalog = load_direction_rule_catalog(
            RULE_CATALOG_PATH,
            known_classes=[item["id"] for item in registry["path_ontology_classes"]],
        )
        self.cached_incremental_scopes: set[tuple[str, ...]] = set()
        self.edge_contribution_cache: dict[tuple[str, str, str], np.ndarray] = {}
        self.downstream_cache: dict[
            tuple[str, str, str],
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str],
        ] = {}

    @staticmethod
    def _json_key(value: Any) -> str:
        return json.dumps(value, sort_keys=True, separators=(",", ":"))

    def _edge_matrix(
        self,
        config: dict[str, Any],
        graph_symbols: list[str],
        graph_metadata: pd.DataFrame,
    ) -> pd.DataFrame:
        """Integrate cached exact per-stream log-odds contributions."""
        definitions = {item["id"]: item for item in self.registry["edge_streams"]}
        active = [
            (stream_id, state)
            for stream_id, state in config["edge_streams"].items()
            if state["enabled"] and state["weight"] > 0
        ]
        primary = [
            (stream_id, state)
            for stream_id, state in active
            if not definitions[stream_id].get("derived")
        ]
        derived = [
            (stream_id, state)
            for stream_id, state in active
            if definitions[stream_id].get("derived")
        ]
        integration_key = self._json_key(
            {
                key: config["edge_integration"][key]
                for key in (
                    "penalize_unsupported",
                    "unsupported_bayes_factor",
                    "continuous_negative_evidence",
                    "continuous_bayes_factor_floor",
                )
            }
        )
        graph_key = self._json_key(graph_symbols)
        prior = float(config["edge_integration"]["prior_probability"])
        log_odds = np.full(
            (len(graph_symbols), len(graph_symbols)),
            math.log(prior / (1.0 - prior)),
            dtype=float,
        )
        for stream_id, state in primary:
            state_key = self._json_key(state)
            cache_key = (graph_key, integration_key, f"{stream_id}:{state_key}")
            contribution = self.edge_contribution_cache.get(cache_key)
            if contribution is None:
                one_stream = copy.deepcopy(config)
                one_stream["edge_integration"]["prior_probability"] = 0.5
                for candidate_id in one_stream["edge_streams"]:
                    one_stream["edge_streams"][candidate_id]["enabled"] = candidate_id == stream_id
                collector: dict[str, np.ndarray] = {}
                combine_edge_factors(
                    self.project,
                    self.registry,
                    one_stream,
                    graph_symbols,
                    graph_metadata=graph_metadata,
                    contribution_collector=collector,
                )
                contribution = collector[stream_id]
                self.edge_contribution_cache[cache_key] = contribution
            log_odds += contribution
        if derived:
            preclosure_probabilities = stable_expit(log_odds)
            np.fill_diagonal(preclosure_probabilities, 0.0)
            index = {symbol: position for position, symbol in enumerate(graph_symbols)}
            for stream_id, state in derived:
                definition = definitions[stream_id]
                if definition["normalization"]["handler"] != "scaffold_triadic_closure":
                    raise ValueError(
                        f"unsupported derived edge handler: "
                        f"{definition['normalization']['handler']}"
                    )
                table, _ = scaffold_triadic_closure_factors(
                    self.project,
                    graph_symbols,
                    preclosure_probabilities,
                    state,
                    graph_metadata=graph_metadata,
                )
                if table.empty:
                    continue
                left = table["node_a"].map(index).astype(int).to_numpy()
                right = table["node_b"].map(index).astype(int).to_numpy()
                weighted_logs = float(state["weight"]) * np.log(
                    table["bayes_factor"].to_numpy(float)
                )
                np.add.at(log_odds, (left, right), weighted_logs)
                np.add.at(log_odds, (right, left), weighted_logs)
        probabilities = stable_expit(log_odds)
        np.fill_diagonal(probabilities, 0.0)
        matrix = pd.DataFrame(probabilities, index=graph_symbols, columns=graph_symbols)
        matrix.index.name = "symbol"
        return matrix

    def evaluate(self, supplied_config: dict[str, Any]) -> Snapshot:
        config = normalize_configuration(supplied_config, self.registry)
        # Ablation measures fixed fitted settings. Refitting after removing a stream
        # would answer a different question and allow remaining streams to compensate.
        config["calibration"]["enabled"] = False
        node_factors, selected_nodes, _ = select_nodes(self.project, self.registry, config)
        graph_symbols = selected_nodes["symbol"].astype(str).tolist()
        graph_metadata = selected_nodes.copy()
        target_symbol = ""
        if config["path"]["enabled"]:
            target_symbol, graph_metadata, _ = append_requested_graph_endpoint(
                config["path"]["target"],
                graph_symbols,
                graph_metadata,
                self.universe,
                node_factors,
                self.project,
                reason="ablation_path_endpoint",
            )
        downstream_key = (
            self._json_key(graph_symbols),
            self._json_key(
                {
                    "edge_streams": config["edge_streams"],
                    "edge_integration": config["edge_integration"],
                }
            ),
            self._json_key(config["path"]),
        )
        cached_downstream = self.downstream_cache.get(downstream_key)
        if cached_downstream is not None:
            matrix, propagation, supported, paths, path_status = cached_downstream
            return Snapshot(
                config=config,
                node_factors=node_factors,
                selected_nodes=selected_nodes,
                matrix=matrix,
                propagation_matrix=propagation,
                supported_edges=supported,
                paths=paths,
                path_status=path_status,
                target_symbol=target_symbol,
            )
        incremental = tuple(
            sorted(
                (symbol for symbol in graph_symbols if symbol not in self.seed_set),
                key=lambda value: (value.casefold(), value),
            )
        )
        if incremental and incremental not in self.cached_incremental_scopes:
            ensure_incremental_pairs(
                self.project,
                graph_metadata,
                self.seed_symbols,
            )
            self.cached_incremental_scopes.add(incremental)
        matrix = self._edge_matrix(config, graph_symbols, graph_metadata)
        cutoff = float(config["edge_integration"]["output_probability_cutoff"])
        supported = supported_edge_table(matrix, cutoff)
        propagation = matrix
        directed = False
        if config["path"]["enabled"] and config["path"]["ontology_directionality_enabled"]:
            orientation_cutoff = min(
                cutoff,
                float(config["path"]["minimum_edge_probability"]),
            )
            propagation, _, _ = apply_ontology_directionality(
                matrix,
                graph_metadata,
                self.direction_catalog,
                audit_probability_cutoff=orientation_cutoff,
                edge_output_cutoff=cutoff,
                path_probability_cutoff=float(config["path"]["minimum_edge_probability"]),
            )
            directed = True
        paths = pd.DataFrame()
        path_status = "disabled"
        if config["path"]["enabled"]:
            try:
                paths, _, _, _ = run_paths(
                    propagation,
                    graph_metadata,
                    config,
                    config["path"]["start"],
                    target_symbol,
                    directed_graph=directed,
                )
                path_status = "ok" if len(paths) else "no_supported_path"
            except ValueError as exc:
                path_status = str(exc)
        self.downstream_cache[downstream_key] = (
            matrix,
            propagation,
            supported,
            paths,
            path_status,
        )
        return Snapshot(
            config=config,
            node_factors=node_factors,
            selected_nodes=selected_nodes,
            matrix=matrix,
            propagation_matrix=propagation,
            supported_edges=supported,
            paths=paths,
            path_status=path_status,
            target_symbol=target_symbol,
        )


def compare_snapshots(
    full: Snapshot,
    reduced: Snapshot,
    *,
    stage: str,
    stream_id: str,
    stream_label: str,
    control_label: str,
    multiplier: float,
    reference_multiplier: float,
    analysis_context: str,
    enabled_in_supplied_configuration: bool,
) -> dict[str, Any]:
    full_nodes = full.node_factors.set_index("gene_symbol")["gui_posterior"].astype(float)
    reduced_nodes = reduced.node_factors.set_index("gene_symbol")["gui_posterior"].astype(float)
    common_candidates = full_nodes.index.intersection(reduced_nodes.index)
    full_node_values = full_nodes.loc[common_candidates].to_numpy(float)
    reduced_node_values = reduced_nodes.loc[common_candidates].to_numpy(float)
    node_delta = full_node_values - reduced_node_values
    node_kl = bernoulli_kl_bits(full_node_values, reduced_node_values)
    full_graph_nodes = set(full.matrix.index.astype(str))
    reduced_graph_nodes = set(reduced.matrix.index.astype(str))
    graph_union = full_graph_nodes | reduced_graph_nodes
    graph_intersection = full_graph_nodes & reduced_graph_nodes

    full_supported = supported_pair_set(full.supported_edges)
    reduced_supported = supported_pair_set(reduced.supported_edges)
    full_common_supported = {
        pair for pair in full_supported if pair[0] in graph_intersection and pair[1] in graph_intersection
    }
    reduced_common_supported = {
        pair for pair in reduced_supported if pair[0] in graph_intersection and pair[1] in graph_intersection
    }
    common_symbols = sorted(graph_intersection, key=lambda value: (value.casefold(), value))
    full_edge_map = edge_probability_map(full.matrix, common_symbols)
    reduced_edge_map = edge_probability_map(reduced.matrix, common_symbols)
    common_pairs = sorted(full_edge_map)
    full_edge_values = np.asarray([full_edge_map[pair] for pair in common_pairs], dtype=float)
    reduced_edge_values = np.asarray([reduced_edge_map[pair] for pair in common_pairs], dtype=float)
    edge_delta = full_edge_values - reduced_edge_values
    edge_kl = bernoulli_kl_bits(full_edge_values, reduced_edge_values)

    full_paths = path_map(full.paths)
    reduced_paths = path_map(reduced.paths)
    full_path_set = set(full_paths)
    reduced_path_set = set(reduced_paths)
    common_paths = full_path_set & reduced_path_set
    full_top_route = str(full.paths.iloc[0]["path_symbols"]) if len(full.paths) else ""
    reduced_top_route = str(reduced.paths.iloc[0]["path_symbols"]) if len(reduced.paths) else ""
    common_path_score_shift = [
        abs(float(full_paths[route]["score"]) - float(reduced_paths[route]["score"]))
        for route in common_paths
    ]
    original_top_rank = (
        int(reduced_paths[full_top_route]["rank"])
        if full_top_route and full_top_route in reduced_paths
        else math.nan
    )
    path_rank_correlation = safe_spearman(
        np.asarray([float(full_paths[route]["rank"]) for route in sorted(common_paths)]),
        np.asarray([float(reduced_paths[route]["rank"]) for route in sorted(common_paths)]),
    )

    return {
        "stage": stage,
        "stream_id": stream_id,
        "stream_label": stream_label,
        "analysis_context": analysis_context,
        "enabled_in_supplied_configuration": enabled_in_supplied_configuration,
        "normalization_control": control_label,
        "full_model_focal_multiplier": float(multiplier),
        "configured_reference_multiplier": float(reference_multiplier),
        "is_configured_reference_multiplier": bool(
            abs(multiplier - reference_multiplier) <= REFERENCE_TOLERANCE
        ),
        "node_candidate_count": int(len(common_candidates)),
        "full_selected_graph_nodes": int(len(full_graph_nodes)),
        "ablated_selected_graph_nodes": int(len(reduced_graph_nodes)),
        "selected_graph_nodes_lost": int(len(full_graph_nodes - reduced_graph_nodes)),
        "selected_graph_nodes_gained": int(len(reduced_graph_nodes - full_graph_nodes)),
        "selected_graph_node_jaccard": (
            float(len(graph_intersection) / len(graph_union)) if graph_union else 1.0
        ),
        "node_sum_absolute_posterior_shift": float(np.abs(node_delta).sum()),
        "node_mean_absolute_posterior_shift": float(np.abs(node_delta).mean()),
        "node_maximum_absolute_posterior_shift": float(np.abs(node_delta).max()),
        "node_sum_conditional_kl_bits": float(node_kl.sum()),
        "node_mean_conditional_kl_bits": float(node_kl.mean()),
        "node_decision_flip_count": int(
            ((full_node_values > full.config["node_integration"]["output_probability_cutoff"]) !=
             (reduced_node_values > reduced.config["node_integration"]["output_probability_cutoff"])).sum()
        ),
        "node_posterior_spearman": safe_spearman(full_node_values, reduced_node_values),
        "full_supported_edges": int(len(full_supported)),
        "ablated_supported_edges": int(len(reduced_supported)),
        "supported_edges_lost_total": int(len(full_supported - reduced_supported)),
        "supported_edges_gained_total": int(len(reduced_supported - full_supported)),
        "full_supported_edges_incident_to_nodes_absent_after_ablation": int(
            len(full_supported - full_common_supported)
        ),
        "supported_edges_lost_among_shared_nodes": int(
            len(full_common_supported - reduced_common_supported)
        ),
        "supported_edges_gained_among_shared_nodes": int(
            len(reduced_common_supported - full_common_supported)
        ),
        "supported_edge_jaccard_among_shared_nodes": (
            float(len(full_common_supported & reduced_common_supported) /
                  len(full_common_supported | reduced_common_supported))
            if (full_common_supported | reduced_common_supported)
            else 1.0
        ),
        "shared_node_unique_edge_hypotheses": int(len(common_pairs)),
        "edge_sum_absolute_probability_shift_shared_nodes": float(np.abs(edge_delta).sum()),
        "edge_mean_absolute_probability_shift_shared_nodes": (
            float(np.abs(edge_delta).mean()) if len(edge_delta) else math.nan
        ),
        "edge_maximum_absolute_probability_shift_shared_nodes": (
            float(np.abs(edge_delta).max()) if len(edge_delta) else math.nan
        ),
        "edge_sum_conditional_kl_bits_shared_nodes": float(edge_kl.sum()),
        "edge_mean_conditional_kl_bits_shared_nodes": (
            float(edge_kl.mean()) if len(edge_kl) else math.nan
        ),
        "edge_probability_spearman_shared_nodes": safe_spearman(full_edge_values, reduced_edge_values),
        "full_path_status": full.path_status,
        "ablated_path_status": reduced.path_status,
        "full_paths_found": int(len(full.paths)),
        "ablated_paths_found": int(len(reduced.paths)),
        "common_top_k_paths": int(len(common_paths)),
        "top_k_path_jaccard": (
            float(len(common_paths) / len(full_path_set | reduced_path_set))
            if (full_path_set | reduced_path_set)
            else 1.0
        ),
        "full_top_path": full_top_route,
        "ablated_top_path": reduced_top_route,
        "top_path_changed": bool(full_top_route != reduced_top_route),
        "full_top_path_score": (
            float(full.paths.iloc[0]["path_probability_product"]) if len(full.paths) else math.nan
        ),
        "ablated_top_path_score": (
            float(reduced.paths.iloc[0]["path_probability_product"]) if len(reduced.paths) else math.nan
        ),
        "full_top_path_rank_after_ablation": original_top_rank,
        "common_path_mean_absolute_score_shift": (
            float(np.mean(common_path_score_shift)) if common_path_score_shift else math.nan
        ),
        "common_path_rank_spearman": path_rank_correlation,
    }


def reference_node_change_table(
    full: Snapshot,
    reduced: Snapshot,
    stage: str,
    stream_id: str,
    label: str,
    analysis_context: str,
) -> pd.DataFrame:
    full_values = full.node_factors.set_index("gene_symbol")["gui_posterior"].astype(float)
    reduced_values = reduced.node_factors.set_index("gene_symbol")["gui_posterior"].astype(float)
    frame = pd.concat(
        [full_values.rename("full_posterior"), reduced_values.rename("ablated_posterior")],
        axis=1,
    )
    frame["posterior_change_full_minus_ablated"] = frame["full_posterior"] - frame["ablated_posterior"]
    frame["absolute_posterior_change"] = frame["posterior_change_full_minus_ablated"].abs()
    frame["stage"] = stage
    frame["stream_id"] = stream_id
    frame["stream_label"] = label
    frame["analysis_context"] = analysis_context
    return frame.reset_index().sort_values("absolute_posterior_change", ascending=False)


def reference_edge_change_table(
    full: Snapshot,
    reduced: Snapshot,
    stage: str,
    stream_id: str,
    label: str,
    analysis_context: str,
) -> pd.DataFrame:
    common = sorted(
        set(full.matrix.index.astype(str)) & set(reduced.matrix.index.astype(str)),
        key=lambda value: (value.casefold(), value),
    )
    full_map = edge_probability_map(full.matrix, common)
    reduced_map = edge_probability_map(reduced.matrix, common)
    cutoff = float(full.config["edge_integration"]["output_probability_cutoff"])
    rows = []
    for pair, full_probability in full_map.items():
        reduced_probability = reduced_map[pair]
        full_supported = full_probability > cutoff
        reduced_supported = reduced_probability > cutoff
        if full_supported != reduced_supported:
            rows.append(
                {
                    "stage": stage,
                    "stream_id": stream_id,
                    "stream_label": label,
                    "analysis_context": analysis_context,
                    "node_a": pair[0],
                    "node_b": pair[1],
                    "full_probability": full_probability,
                    "ablated_probability": reduced_probability,
                    "probability_change_full_minus_ablated": full_probability - reduced_probability,
                    "change_type": "lost_after_ablation" if full_supported else "gained_after_ablation",
                }
            )
    return pd.DataFrame(rows)


def reference_path_table(
    full: Snapshot,
    reduced: Snapshot,
    stage: str,
    stream_id: str,
    label: str,
    analysis_context: str,
) -> pd.DataFrame:
    full_map = path_map(full.paths)
    reduced_map = path_map(reduced.paths)
    routes = sorted(set(full_map) | set(reduced_map), key=lambda route: (
        full_map.get(route, {}).get("rank", 10**9),
        reduced_map.get(route, {}).get("rank", 10**9),
        route,
    ))
    return pd.DataFrame(
        [
            {
                "stage": stage,
                "stream_id": stream_id,
                "stream_label": label,
                "analysis_context": analysis_context,
                "path_symbols": route,
                "full_rank": full_map.get(route, {}).get("rank", math.nan),
                "ablated_rank": reduced_map.get(route, {}).get("rank", math.nan),
                "full_score": full_map.get(route, {}).get("score", math.nan),
                "ablated_score": reduced_map.get(route, {}).get("score", math.nan),
                "membership": (
                    "both" if route in full_map and route in reduced_map
                    else "full_only" if route in full_map else "ablated_only"
                ),
            }
            for route in routes
        ]
    )


def load_font(size: int, *, bold: bool = False) -> ImageFont.ImageFont:
    windows_font = Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf")
    try:
        return ImageFont.truetype(str(windows_font), size=size)
    except OSError:
        return ImageFont.load_default()


def abbreviate(label: str, maximum: int = 32) -> str:
    return label if len(label) <= maximum else label[: maximum - 1] + "…"


def draw_reference_impact(summary: pd.DataFrame, output: Path) -> None:
    """Render three aligned ablation-impact panels as a PNG."""
    if summary.empty:
        return
    data = summary.sort_values(
        ["stage", "node_sum_conditional_kl_bits"], ascending=[True, False]
    ).reset_index(drop=True)
    width = 1800
    row_height = 58
    top = 150
    bottom = 80
    height = top + bottom + row_height * len(data)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(34, bold=True)
    subtitle_font = load_font(20)
    label_font = load_font(18)
    small_font = load_font(16)
    draw.text((42, 28), "End-to-end evidence-stream ablation at configured Tq", fill="#111827", font=title_font)
    draw.text(
        (42, 78),
        "Bars show selected nodes lost, supported edges lost, and loss of top-k path agreement after omitting one stream.",
        fill="#4b5563",
        font=subtitle_font,
    )
    label_x = 42
    panel_x = [530, 940, 1350]
    panel_width = 320
    columns = [
        ("selected_graph_nodes_lost", "Nodes lost", "#2563eb"),
        ("supported_edges_lost_total", "Edges lost", "#7c3aed"),
        ("top_k_path_jaccard_loss", "1 − path Jaccard", "#dc2626"),
    ]
    data = data.copy()
    data["top_k_path_jaccard_loss"] = 1.0 - data["top_k_path_jaccard"].astype(float)
    for x, (_, heading, _) in zip(panel_x, columns):
        draw.text((x, 115), heading, fill="#111827", font=label_font)
    maxima = {
        column: max(float(data[column].max()), 1.0 if column != "top_k_path_jaccard_loss" else 0.01)
        for column, _, _ in columns
    }
    for row_index, row in data.iterrows():
        y = top + row_index * row_height
        if row_index % 2 == 0:
            draw.rectangle((20, y - 5, width - 20, y + row_height - 7), fill="#f8fafc")
        stage_tag = "NODE" if row["stage"] == "node" else "EDGE"
        draw.text((label_x, y + 7), f"{stage_tag}  {abbreviate(str(row['stream_label']))}", fill="#111827", font=label_font)
        for x, (column, _, color) in zip(panel_x, columns):
            value = float(row[column])
            bar_width = int(panel_width * value / maxima[column]) if maxima[column] else 0
            draw.rectangle((x, y + 9, x + bar_width, y + 35), fill=color)
            shown = f"{value:.3f}" if column == "top_k_path_jaccard_loss" else f"{int(value):,}"
            draw.text((x + min(bar_width + 8, panel_width - 5), y + 8), shown, fill="#111827", font=small_font)
    image.save(output, format="PNG")


def draw_tq_robustness(sweep: pd.DataFrame, output: Path, stage: str) -> None:
    subset = sweep.loc[sweep["stage"].eq(stage)].copy()
    if subset.empty:
        return
    streams = list(dict.fromkeys(subset["stream_id"].astype(str)))
    width = 1800
    panel_height = 230
    top = 130
    height = top + panel_height * len(streams) + 60
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(32, bold=True)
    heading_font = load_font(20, bold=True)
    small_font = load_font(15)
    draw.text((42, 25), f"{stage.title()}-stream ablation across focal Tq/reference values", fill="#111827", font=title_font)
    draw.text(
        (42, 72),
        "Each curve compares the full model at the shown focal multiplier with the same model lacking that stream; all other settings are fixed.",
        fill="#4b5563",
        font=load_font(18),
    )
    plot_left, plot_right = 470, 1730
    metrics = [
        ("selected_graph_nodes_lost", "nodes lost", "#2563eb"),
        ("supported_edges_lost_total", "edges lost", "#7c3aed"),
        ("top_k_path_jaccard", "path Jaccard", "#dc2626"),
    ]
    for stream_index, stream_id in enumerate(streams):
        group = subset.loc[subset["stream_id"].eq(stream_id)].sort_values("full_model_focal_multiplier")
        y0 = top + stream_index * panel_height
        draw.text((42, y0 + 15), abbreviate(str(group.iloc[0]["stream_label"]), 40), fill="#111827", font=heading_font)
        normalization_control = group.iloc[0]["normalization_control"]
        control_label = "derived" if pd.isna(normalization_control) else str(normalization_control)
        draw.text((42, y0 + 52), control_label, fill="#6b7280", font=small_font)
        x_values = group["full_model_focal_multiplier"].to_numpy(float)
        x_log = np.log10(x_values)
        x_min, x_max = float(x_log.min()), float(x_log.max())
        x_span = max(x_max - x_min, 1e-9)
        band_width = (plot_right - plot_left - 40) / len(metrics)
        for metric_index, (column, label, color) in enumerate(metrics):
            left = plot_left + metric_index * band_width
            right = left + band_width - 45
            top_metric = y0 + 24
            bottom_metric = y0 + 172
            values = group[column].to_numpy(float)
            if column == "top_k_path_jaccard":
                y_min, y_max = 0.0, 1.0
            else:
                y_min, y_max = 0.0, max(float(np.nanmax(values)), 1.0)
            draw.line((left, bottom_metric, right, bottom_metric), fill="#9ca3af", width=1)
            draw.line((left, top_metric, left, bottom_metric), fill="#9ca3af", width=1)
            draw.text((left, y0 + 184), label, fill=color, font=small_font)
            points = []
            for x_value, x_position, value in zip(x_values, x_log, values):
                px = left + (right - left) * (x_position - x_min) / x_span
                py = bottom_metric - (bottom_metric - top_metric) * (float(value) - y_min) / max(y_max - y_min, 1e-9)
                points.append((px, py))
                draw.ellipse((px - 4, py - 4, px + 4, py + 4), fill=color)
                if abs(float(x_value) - float(group.iloc[0]["configured_reference_multiplier"])) <= REFERENCE_TOLERANCE:
                    draw.ellipse((px - 8, py - 8, px + 8, py + 8), outline="#111827", width=2)
            if len(points) > 1:
                draw.line(points, fill=color, width=3)
            draw.text((left, bottom_metric + 3), f"{x_values.min():g}", fill="#6b7280", font=small_font)
            draw.text((right - 25, bottom_metric + 3), f"{x_values.max():g}", fill="#6b7280", font=small_font)
        if stream_index < len(streams) - 1:
            draw.line((25, y0 + panel_height - 8, width - 25, y0 + panel_height - 8), fill="#e5e7eb", width=1)
    image.save(output, format="PNG")


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    if frame.empty:
        return "No rows."
    values = frame.loc[:, columns].copy()
    for column in values.columns:
        if pd.api.types.is_float_dtype(values[column]):
            values[column] = values[column].map(lambda value: "" if pd.isna(value) else f"{value:.6g}")
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |"
        for row in values.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def _read_checkpoint_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def import_compatible_reference_checkpoints(
    source: Path,
    checkpoint_directory: Path,
    baseline_config: dict[str, Any],
) -> int:
    """Import completed legacy/default stream results after exact config validation."""
    source_config_path = source / "analyzed_configuration.json"
    if not source_config_path.is_file():
        raise FileNotFoundError(
            f"Resume source has no analyzed_configuration.json: {source}"
        )
    source_config = json.loads(source_config_path.read_text(encoding="utf-8"))
    if source_config != baseline_config:
        raise ValueError(
            "Resume source configuration does not exactly match the current "
            "normalized baseline; refusing to mix incomparable ablations."
        )
    sweep = pd.read_csv(source / "ablation_tq_sweep_summary.tsv", sep="\t")
    node = pd.read_csv(
        source / "reference_node_posterior_changes.tsv.gz", sep="\t"
    )
    edge = pd.read_csv(source / "reference_edge_cutoff_flips.tsv.gz", sep="\t")
    path = pd.read_csv(source / "reference_path_rank_changes.tsv", sep="\t")
    imported = 0
    for (stage, stream_id), group in sweep.groupby(["stage", "stream_id"], sort=False):
        prefix = checkpoint_directory / f"{stage}__{stream_id}"
        paths = {
            "sweep": prefix.with_name(prefix.name + "__sweep.tsv"),
            "node": prefix.with_name(prefix.name + "__node_changes.tsv.gz"),
            "edge": prefix.with_name(prefix.name + "__edge_flips.tsv.gz"),
            "path": prefix.with_name(prefix.name + "__path_changes.tsv"),
        }
        if all(item.is_file() for item in paths.values()):
            continue
        group = group.copy()
        group["analysis_context"] = "supplied_configuration"
        group["enabled_in_supplied_configuration"] = True
        node_group = node.loc[node["stream_id"].astype(str).eq(str(stream_id))].copy()
        edge_group = edge.loc[edge["stream_id"].astype(str).eq(str(stream_id))].copy()
        path_group = path.loc[path["stream_id"].astype(str).eq(str(stream_id))].copy()
        for table in (node_group, edge_group, path_group):
            table["analysis_context"] = "supplied_configuration"
        group.to_csv(paths["sweep"], sep="\t", index=False)
        node_group.to_csv(paths["node"], sep="\t", index=False, compression="gzip")
        edge_group.to_csv(paths["edge"], sep="\t", index=False, compression="gzip")
        path_group.to_csv(paths["path"], sep="\t", index=False)
        imported += 1
    return imported


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = parse_args()
    project = args.project_root.resolve()
    registry = load_registry(project)
    if args.configuration:
        supplied = json.loads(args.configuration.read_text(encoding="utf-8"))
        config_source = str(args.configuration.resolve())
    else:
        supplied = default_configuration(registry)
        config_source = "registered_default_configuration"
    config = normalize_configuration(supplied, registry)
    calibration_was_enabled = bool(config["calibration"]["enabled"])
    config["calibration"]["enabled"] = False
    if not config["path"]["enabled"]:
        raise ValueError("End-to-end ablation requires path finding to be enabled.")
    output = (
        args.output_dir.resolve()
        if args.output_dir
        else project / "results/sensitivity_analysis" / f"end_to_end_ablation_{date.today().isoformat()}"
    )
    output.mkdir(parents=True, exist_ok=True)
    checkpoint_directory = output / "checkpoints"
    checkpoint_directory.mkdir(parents=True, exist_ok=True)
    (output / "supplied_configuration.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )

    node_definitions = {item["id"]: item for item in registry["node_streams"]}
    edge_definitions = {item["id"]: item for item in registry["edge_streams"]}
    analysis_config = copy.deepcopy(config)
    exclusive_choices = {
        group: next(
            (
                stream_id
                for stream_id in members
                if analysis_config["edge_streams"][stream_id]["enabled"]
            ),
            "none",
        )
        for group, members in _exclusive_groups(registry).items()
    }
    (output / "analyzed_configuration.json").write_text(
        json.dumps(analysis_config, indent=2), encoding="utf-8"
    )
    if args.resume_compatible_reference:
        imported = import_compatible_reference_checkpoints(
            args.resume_compatible_reference.resolve(),
            checkpoint_directory,
            analysis_config,
        )
        print(
            f"Imported {imported} compatible per-stream checkpoints from "
            f"{args.resume_compatible_reference.resolve()}",
            flush=True,
        )
    active: list[tuple[str, str, dict[str, Any], dict[str, Any]]] = []
    if args.stream_scope == "all":
        active.extend(
            ("node", stream_id, definition, analysis_config["node_streams"][stream_id])
            for stream_id, definition in node_definitions.items()
        )
        active.extend(
            ("edge", stream_id, definition, analysis_config["edge_streams"][stream_id])
            for stream_id, definition in edge_definitions.items()
        )
    else:
        for stream_id, state in analysis_config["node_streams"].items():
            if state["enabled"] and state["weight"] > 0:
                active.append(("node", stream_id, node_definitions[stream_id], state))
        for stream_id, state in analysis_config["edge_streams"].items():
            if state["enabled"] and state["weight"] > 0:
                active.append(("edge", stream_id, edge_definitions[stream_id], state))

    evaluator = Evaluator(project, registry)
    print(
        f"Evaluating {len(active)} streams with scope={args.stream_scope!r}...",
        flush=True,
    )
    reference_full_by_context: dict[str, Snapshot] = {}
    sweep_rows: list[dict[str, Any]] = []
    node_change_tables: list[pd.DataFrame] = []
    edge_change_tables: list[pd.DataFrame] = []
    path_change_tables: list[pd.DataFrame] = []

    for position, (stage, stream_id, definition, state) in enumerate(active, start=1):
        checkpoint_prefix = checkpoint_directory / f"{stage}__{stream_id}"
        checkpoint_paths = {
            "sweep": checkpoint_prefix.with_name(checkpoint_prefix.name + "__sweep.tsv"),
            "partial_sweep": checkpoint_prefix.with_name(
                checkpoint_prefix.name + "__sweep.partial.tsv"
            ),
            "node": checkpoint_prefix.with_name(checkpoint_prefix.name + "__node_changes.tsv.gz"),
            "edge": checkpoint_prefix.with_name(checkpoint_prefix.name + "__edge_flips.tsv.gz"),
            "path": checkpoint_prefix.with_name(checkpoint_prefix.name + "__path_changes.tsv"),
        }
        if all(
            checkpoint_paths[key].is_file()
            for key in ("sweep", "node", "edge", "path")
        ):
            print(
                f"[{position}/{len(active)}] Resuming completed {stage} stream: "
                f"{definition['label']}",
                flush=True,
            )
            checkpoint_sweep = _read_checkpoint_table(checkpoint_paths["sweep"])
            enabled_in_supplied = bool(config[f"{stage}_streams"][stream_id]["enabled"])
            if (
                args.stream_scope == "all"
                and not enabled_in_supplied
                and args.optional_tq_grid == "reference"
                and "is_configured_reference_multiplier" in checkpoint_sweep.columns
            ):
                # A resumable directory may contain richer compact sweeps from an
                # earlier run.  Keep the aggregate output faithful to the policy
                # requested for this run while preserving those checkpoint files.
                reference_mask = (
                    checkpoint_sweep["is_configured_reference_multiplier"]
                    .astype(str)
                    .str.lower()
                    .eq("true")
                )
                checkpoint_sweep = checkpoint_sweep.loc[reference_mask].copy()
            sweep_rows.extend(checkpoint_sweep.to_dict("records"))
            node_change_tables.append(_read_checkpoint_table(checkpoint_paths["node"]))
            edge_change_tables.append(_read_checkpoint_table(checkpoint_paths["edge"]))
            path_change_tables.append(_read_checkpoint_table(checkpoint_paths["path"]))
            continue
        print(f"[{position}/{len(active)}] Ablating {stage} stream: {definition['label']}", flush=True)
        if args.stream_scope == "all":
            full_config, analysis_context = focal_context_configuration(
                analysis_config,
                registry,
                stage=stage,
                stream_id=stream_id,
            )
        else:
            full_config = copy.deepcopy(analysis_config)
            analysis_context = "supplied_configuration"
        if not full_config[f"{stage}_streams"][stream_id]["enabled"]:
            raise RuntimeError(
                f"Internal sensitivity-analysis error: focal stream {stream_id!r} "
                "is not enabled in its full-model context"
            )
        reduced_config = copy.deepcopy(full_config)
        reduced_config[f"{stage}_streams"][stream_id]["enabled"] = False
        focal_state = full_config[f"{stage}_streams"][stream_id]
        reference_multiplier = float(focal_state.get("tq_multiplier", 1.0))
        has_tq = definition.get("normalization", {}).get("user_control", True) and not definition.get("derived")
        enabled_in_supplied = bool(config[f"{stage}_streams"][stream_id]["enabled"])
        if not has_tq:
            grid = [reference_multiplier]
        elif not enabled_in_supplied and args.optional_tq_grid == "reference":
            grid = [reference_multiplier]
        else:
            grid = tq_grid(
                definition,
                reference_multiplier,
                compact=not enabled_in_supplied,
            )
        partial_sweep = (
            _read_checkpoint_table(checkpoint_paths["partial_sweep"])
            if checkpoint_paths["partial_sweep"].is_file()
            else pd.DataFrame()
        )
        completed_multipliers = (
            partial_sweep["full_model_focal_multiplier"].astype(float).tolist()
            if "full_model_focal_multiplier" in partial_sweep
            else []
        )
        pending_grid = [
            multiplier
            for multiplier in grid
            if not any(
                abs(multiplier - completed) <= REFERENCE_TOLERANCE
                for completed in completed_multipliers
            )
        ]
        reference_tables_complete = all(
            checkpoint_paths[key].is_file() for key in ("node", "edge", "path")
        )
        need_reference_snapshot = (
            not reference_tables_complete
            or any(
                abs(multiplier - reference_multiplier) <= REFERENCE_TOLERANCE
                for multiplier in pending_grid
            )
        )
        reference_full = None
        if need_reference_snapshot:
            reference_full = reference_full_by_context.get(analysis_context)
            if reference_full is None:
                print(f"    Evaluating full context: {analysis_context}", flush=True)
                reference_full = evaluator.evaluate(full_config)
                reference_full_by_context[analysis_context] = reference_full
        reduced = evaluator.evaluate(reduced_config)
        if not reference_tables_complete:
            if reference_full is None:
                raise RuntimeError("Reference snapshot is required for audit tables")
            reference_node_change_table(
                reference_full,
                reduced,
                stage,
                stream_id,
                definition["label"],
                analysis_context,
            ).to_csv(
                checkpoint_paths["node"], sep="\t", index=False, compression="gzip"
            )
            reference_edge_change_table(
                reference_full,
                reduced,
                stage,
                stream_id,
                definition["label"],
                analysis_context,
            ).to_csv(
                checkpoint_paths["edge"], sep="\t", index=False, compression="gzip"
            )
            reference_path_table(
                reference_full,
                reduced,
                stage,
                stream_id,
                definition["label"],
                analysis_context,
            ).to_csv(checkpoint_paths["path"], sep="\t", index=False)
        focal_sweep_rows = partial_sweep.to_dict("records")
        reference_row_complete = any(
            abs(completed - reference_multiplier) <= REFERENCE_TOLERANCE
            for completed in completed_multipliers
        )
        if not reference_row_complete:
            if reference_full is None:
                raise RuntimeError("Reference snapshot is required for the reference Tq row")
            focal_sweep_rows.append(
                compare_snapshots(
                    reference_full,
                    reduced,
                    stage=stage,
                    stream_id=stream_id,
                    stream_label=definition["label"],
                    control_label=definition.get("normalization", {}).get(
                        "control_label", "not applicable"
                    ),
                    multiplier=reference_multiplier,
                    reference_multiplier=reference_multiplier,
                    analysis_context=analysis_context,
                    enabled_in_supplied_configuration=enabled_in_supplied,
                )
            )
            pd.DataFrame(focal_sweep_rows).to_csv(
                checkpoint_paths["partial_sweep"], sep="\t", index=False
            )
            pending_grid = [
                multiplier
                for multiplier in pending_grid
                if abs(multiplier - reference_multiplier) > REFERENCE_TOLERANCE
            ]
        for multiplier in pending_grid:
            print(
                f"    {definition.get('normalization', {}).get('control_label', 'scale')} "
                f"{multiplier:g}",
                flush=True,
            )
            if abs(multiplier - reference_multiplier) <= REFERENCE_TOLERANCE:
                if reference_full is None:
                    raise RuntimeError("Reference snapshot missing for reference Tq row")
                full = reference_full
            else:
                varied_config = copy.deepcopy(full_config)
                varied_config[f"{stage}_streams"][stream_id]["tq_multiplier"] = float(multiplier)
                full = evaluator.evaluate(varied_config)
            focal_sweep_rows.append(
                compare_snapshots(
                    full,
                    reduced,
                    stage=stage,
                    stream_id=stream_id,
                    stream_label=definition["label"],
                    control_label=definition.get("normalization", {}).get("control_label", "not applicable"),
                    multiplier=multiplier,
                    reference_multiplier=reference_multiplier,
                    analysis_context=analysis_context,
                    enabled_in_supplied_configuration=enabled_in_supplied,
                )
            )
            pd.DataFrame(focal_sweep_rows).to_csv(
                checkpoint_paths["partial_sweep"], sep="\t", index=False
            )
        focal_sweep = pd.DataFrame(focal_sweep_rows).sort_values(
            "full_model_focal_multiplier"
        )
        focal_sweep.to_csv(checkpoint_paths["sweep"], sep="\t", index=False)
        focal_node_changes = _read_checkpoint_table(checkpoint_paths["node"])
        focal_edge_changes = _read_checkpoint_table(checkpoint_paths["edge"])
        focal_path_changes = _read_checkpoint_table(checkpoint_paths["path"])
        sweep_rows.extend(focal_sweep_rows)
        node_change_tables.append(focal_node_changes)
        edge_change_tables.append(focal_edge_changes)
        path_change_tables.append(focal_path_changes)
        if stage == "node":
            reference_full_by_context.pop(analysis_context, None)
            evaluator.edge_contribution_cache.clear()
            evaluator.downstream_cache.clear()
            gc.collect()

    sweep = pd.DataFrame(sweep_rows)
    reference = sweep.loc[sweep["is_configured_reference_multiplier"]].copy()
    sweep.to_csv(output / "ablation_tq_sweep_summary.tsv", sep="\t", index=False)
    reference.to_csv(output / "reference_tq_ablation_summary.tsv", sep="\t", index=False)
    pd.concat(node_change_tables, ignore_index=True).to_csv(
        output / "reference_node_posterior_changes.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    nonempty_edge_changes = [table for table in edge_change_tables if not table.empty]
    edge_changes = pd.concat(nonempty_edge_changes, ignore_index=True) if nonempty_edge_changes else pd.DataFrame(
        columns=["stage", "stream_id", "stream_label", "analysis_context", "node_a", "node_b", "full_probability", "ablated_probability", "probability_change_full_minus_ablated", "change_type"]
    )
    edge_changes.to_csv(
        output / "reference_edge_cutoff_flips.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    pd.concat(path_change_tables, ignore_index=True).to_csv(
        output / "reference_path_rank_changes.tsv", sep="\t", index=False
    )

    robustness = (
        sweep.groupby(
            ["stage", "stream_id", "stream_label", "analysis_context"],
            as_index=False,
        )
        .agg(
            tq_grid_min=("full_model_focal_multiplier", "min"),
            tq_grid_max=("full_model_focal_multiplier", "max"),
            tq_grid_points=("full_model_focal_multiplier", "count"),
            min_nodes_lost=("selected_graph_nodes_lost", "min"),
            max_nodes_lost=("selected_graph_nodes_lost", "max"),
            min_edges_lost=("supported_edges_lost_total", "min"),
            max_edges_lost=("supported_edges_lost_total", "max"),
            min_path_jaccard=("top_k_path_jaccard", "min"),
            max_path_jaccard=("top_k_path_jaccard", "max"),
            top_path_changed_at_any_tq=("top_path_changed", "max"),
            top_path_changed_at_every_tq=("top_path_changed", "min"),
        )
    )
    robustness.to_csv(output / "tq_robustness_ranges.tsv", sep="\t", index=False)

    draw_reference_impact(reference, output / "reference_tq_end_to_end_ablation.png")
    draw_tq_robustness(sweep, output / "node_stream_tq_aware_ablation.png", "node")
    draw_tq_robustness(sweep, output / "edge_stream_tq_aware_ablation.png", "edge")

    ranked = reference.sort_values(
        ["top_k_path_jaccard", "supported_edges_lost_total", "selected_graph_nodes_lost"],
        ascending=[True, False, False],
    )
    configured_focal = next(
        (
            (stage, stream_id)
            for stage, stream_id, _, _ in active
            if config[f"{stage}_streams"][stream_id]["enabled"]
        ),
        (active[0][0], active[0][1]),
    )
    primary_config, primary_context = focal_context_configuration(
        analysis_config,
        registry,
        stage=configured_focal[0],
        stream_id=configured_focal[1],
    )
    primary_reference_full = reference_full_by_context.get(primary_context)
    if primary_reference_full is None:
        primary_reference_full = evaluator.evaluate(primary_config)
        reference_full_by_context[primary_context] = primary_reference_full
    optional_included = [
        stream_id
        for stage, stream_id, _, _ in active
        if not config[f"{stage}_streams"][stream_id]["enabled"]
    ]
    report = f"""# Tq-aware end-to-end ablation sensitivity analysis

Date: {date.today().isoformat()}

## Question

For each evidence stream in the selected scope (`{args.stream_scope}`), what changes when that stream is removed from a matched Bayesian workflow containing that stream? The analysis measures effects on node posteriors and selection, the rebuilt edge graph, ontology-constrained path propagation, and the ranked top-{analysis_config['path']['top_k']} paths from `{analysis_config['path']['start']}` to `{analysis_config['path']['target']}`.

## Comparison design

1. The **reference comparison** disables one stream at its currently configured Tq/reference multiplier. A stream already enabled in the supplied configuration is removed from that configuration. A stream that was off is first added by itself to the supplied configuration, then removed in the comparator. This covers every registered stream without constructing an artificial model in which every optional assay is forced on simultaneously. All other weights, Tq values, priors, probability cutoffs, directionality rules, and path settings are unchanged.
2. The **Tq-aware comparison** varies only the focal stream's multiplier while comparing with the same ablated model. Streams enabled in the supplied profile retain their registered calibration sweep. Optional streams use policy `{args.optional_tq_grid}`: the default `reference` policy performs their expensive end-to-end comparison at the configured Tq, while the separate all-stream Tq-response analysis supplies dense Tq-vs-BF and Tq-vs-information curves without expanding every optional node graph. The opt-in `compact` policy adds low/high registered bounds and a distinct preferred value. STRING/STITCH use a reference-score multiplier, so their control is labeled `Ref x` rather than a Gaussian Tq.
3. Calibration is not rerun after ablation. Allowing remaining streams to refit would measure compensation by the optimizer rather than the marginal importance of the omitted evidence. The analyzed file records fixed weights and fixed non-focal Tq values.
4. Node-stream ablation is propagated end to end. It can remove graph nodes and consequently every incident edge before edge evidence is integrated. Edge-stream ablation does not alter node selection, but it rebuilds the edge probabilities and paths.
5. Mutually exclusive alternatives are evaluated in matched contexts. For example, optional HPA high-confidence is tested in a full model where it replaces HPA primary; its ablated comparator contains neither HPA alternative. The `analysis_context` and `enabled_in_supplied_configuration` columns record the comparison design explicitly.

The configured output probability cutoffs are strict: a node or edge is supported only when its posterior probability is greater than its cutoff, not equal to it.

## Reference-Tq results

{markdown_table(ranked, ['stage', 'stream_label', 'selected_graph_nodes_lost', 'supported_edges_lost_total', 'supported_edge_jaccard_among_shared_nodes', 'top_k_path_jaccard', 'top_path_changed', 'node_sum_conditional_kl_bits', 'edge_sum_conditional_kl_bits_shared_nodes'])}

## Metric interpretation

- `selected_graph_nodes_lost` is the number of nodes present in the full graph but absent after ablation. The target endpoint is included in both graphs.
- `supported_edges_lost_total` counts unique undirected pairs, once per pair. `full_supported_edges_incident_to_nodes_absent_after_ablation` separates edge loss caused mechanically by node removal from probability-cutoff changes among shared nodes.
- Posterior/probability shifts are signed as full minus ablated; the summary uses absolute shifts.
- `sum_conditional_kl_bits` is `sum_i D_KL(Bernoulli(p_full,i) || Bernoulli(p_ablated,i))`. It measures how distinguishable the full posterior is from the ablated posterior, conditional on all other active streams. It is not an expected mutual information estimate.
- `top_k_path_jaccard` is the intersection divided by the union of exact path routes in the two top-k lists. A value of 1 means identical route membership; 0 means no shared routes. Rank and score changes for individual routes are in `reference_path_rank_changes.tsv`.
- A large effect at only one Tq but a small effect elsewhere is Tq-sensitive. A consistently large effect across `tq_robustness_ranges.tsv` is more robust evidence that the stream is influential.

## Files

- `supplied_configuration.json`: original normalized configuration before sensitivity-scope expansion.
- `analyzed_configuration.json`: supplied baseline configuration used to create the matched contexts.
- `reference_tq_ablation_summary.tsv`: one row per active stream at its configured multiplier.
- `ablation_tq_sweep_summary.tsv`: all focal Tq/reference comparisons.
- `tq_robustness_ranges.tsv`: min/max effects across each stream's multiplier grid.
- `reference_node_posterior_changes.tsv.gz`: per-candidate node posterior changes for every reference ablation.
- `reference_edge_cutoff_flips.tsv.gz`: shared-node edge pairs that cross the graph cutoff.
- `reference_path_rank_changes.tsv`: path membership, rank, and score changes.
- Three PNGs summarize reference effects and Tq robustness.

## Reproducibility notes

- Configuration source: `{config_source}`.
- Stream scope: `{args.stream_scope}`.
- Optional-stream end-to-end Tq policy: `{args.optional_tq_grid}`.
- Calibration was enabled in the supplied configuration: `{calibration_was_enabled}`. It was deliberately not refit during ablation.
- Analyzed node streams: {', '.join(stream_id for stage, stream_id, _, _ in active if stage == 'node')}.
- Analyzed edge streams: {', '.join(stream_id for stage, stream_id, _, _ in active if stage == 'edge')}.
- Optional streams evaluated in matched add-one contexts: {', '.join(optional_included) if optional_included else 'none'}.
- Primary exclusive-group choices: {json.dumps(exclusive_choices, sort_keys=True)}.
- Supplied-model reference graph: {len(primary_reference_full.matrix):,} nodes and {len(primary_reference_full.supported_edges):,} supported unique undirected edges.
- Primary full reference paths found: {len(primary_reference_full.paths):,}; status `{primary_reference_full.path_status}`.
"""
    (output / "README.md").write_text(report, encoding="utf-8")
    summary = {
        "analysis_date": date.today().isoformat(),
        "configuration_source": config_source,
        "stream_scope": args.stream_scope,
        "optional_stream_end_to_end_tq_policy": args.optional_tq_grid,
        "calibration_refit_during_ablation": False,
        "analyzed_node_stream_count": sum(stage == "node" for stage, *_ in active),
        "analyzed_edge_stream_count": sum(stage == "edge" for stage, *_ in active),
        "registered_node_stream_count": len(node_definitions),
        "registered_edge_stream_count": len(edge_definitions),
        "optional_streams_evaluated_in_matched_contexts": optional_included,
        "primary_exclusive_group_choices": exclusive_choices,
        "analysis_contexts": sorted(
            reference["analysis_context"].dropna().astype(str).unique().tolist()
        ),
        "primary_full_graph_nodes": int(len(primary_reference_full.matrix)),
        "primary_full_supported_unique_undirected_edges": int(len(primary_reference_full.supported_edges)),
        "primary_full_paths_found": int(len(primary_reference_full.paths)),
        "primary_full_path_status": primary_reference_full.path_status,
        "reference_top_path": (
            str(primary_reference_full.paths.iloc[0]["path_symbols"])
            if len(primary_reference_full.paths)
            else None
        ),
        "output_directory": str(output),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Analysis complete: {output}", flush=True)


if __name__ == "__main__":
    main()
