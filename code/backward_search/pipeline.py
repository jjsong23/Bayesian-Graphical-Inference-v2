#!/usr/bin/env python3
"""Unified node-selection -> cheap-edge -> backward-search Version 2 pipeline.

This module connects the configurable workbench evidence model to the bounded
frontier search.  Node selection and inexpensive edge evidence are evaluated
with the same registry, Tq/reference controls, weights, priors, and negative-
evidence policies as Version 1.  AlphaPulldown is still requested only for the
shortlisted frontier pairs produced by Version 2.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd


CODE_DIRECTORY = Path(__file__).resolve().parents[1]
PROJECT_DIRECTORY = CODE_DIRECTORY.parent
GUI_DIRECTORY = PROJECT_DIRECTORY / "gui"
EDGE_DIRECTORY = CODE_DIRECTORY / "edge_characterization"
PATH_DIRECTORY = CODE_DIRECTORY / "path_finding"
for search_path in (GUI_DIRECTORY, CODE_DIRECTORY, EDGE_DIRECTORY, PATH_DIRECTORY):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from incremental_edge_cache import ensure_incremental_pairs  # noqa: E402
from ontology_directionality import (  # noqa: E402
    RULE_CATALOG_PATH,
    apply_ontology_directionality,
    load_direction_rule_catalog,
)
from workflow_engine import (  # noqa: E402
    OMNIPATH_DIRECTION_RAW_RELATIVE,
    UNIPROT_RELATIVE,
    UNIVERSE_RELATIVE,
    append_requested_graph_endpoint,
    build_omnipath_direction_evidence as workflow_build_omnipath_direction_evidence,
    calibrate_edge_parameters,
    calibrate_node_parameters,
    combine_edge_factors,
    default_configuration,
    load_registry,
    normalize_configuration,
    select_nodes,
)

from .development import record_structural_import, render_development_report
from .legacy_screen import LEGACY_PROTOCOL_ID, import_legacy_ppi_screen
from .huri import build_huri_search_evidence
from .engine import (
    build_node_sequence_fasta,
    collect_alphapulldown_scores,
    import_structural_scores,
    initialize_run,
    load_run_state,
    step_run,
    utc_now,
)


PIPELINE_SCHEMA_VERSION = 1
ProgressCallback = Callable[[str, float], None]
CancelCallback = Callable[[], bool]


DEFAULT_BACKWARD_OPTIONS: dict[str, Any] = {
    "target": "Aqp2",
    "receptor": "Avpr2",
    "cheap_top_n_per_frontier": 20,
    "beam_width": 5,
    "maximum_depth": 6,
    "always_include_receptor": True,
    "stop_when_receptor_reached": True,
    "development_mode": True,
    "huri": {
        "enabled": False,
        "interactions_file": "runtime/reference/huri/HuRI.tsv",
        "mapping_file": "data/edge_characterization/localization/hpa/v25.1/processed/mouse_human_hpa_mapping.tsv",
        "weight": 1.0,
        "positive_bayes_factor": 5.0,
        "nonreported_scope": "neutral",
        "nonreported_bayes_factor": 0.9,
        "screened_genes_file": "",
    },
    "structural": {
        "enabled": True,
        "metric": "iptm",
        "reference_score": 0.6,
        "weight": 1.0,
        "bayes_factor_floor": 0.1,
        "bayes_factor_ceiling": 10.0,
        "protocol_id": "alphapulldown_biowulf_v2_model1_cycles3_predictions1_compact",
        "compatible_protocol_ids": ["legacy_ppi_screen_5models_predictions1"],
        "num_cycle": 3,
        "num_predictions_per_model": 1,
        "model_names": "model_1_multimer_v3",
        "maximum_concurrent_jobs": 4,
        "legacy_screen_directory": "../ppi_screen",
        "legacy_protocol_id": LEGACY_PROTOCOL_ID,
        "require_legacy_ranked_structure": True,
    },
}


def _write_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, default=str) + "\n", encoding="utf-8")
    temporary.replace(path)


def _pipeline_path(run_dir: Path) -> Path:
    return run_dir / "pipeline_state.json"


def load_pipeline_state(run_dir: Path) -> dict[str, Any]:
    path = _pipeline_path(run_dir.resolve())
    if not path.is_file():
        raise FileNotFoundError(f"pipeline state is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _write_pipeline_state(run_dir: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = utc_now()
    _write_json(_pipeline_path(run_dir), state)


def _check_cancel(cancel_requested: CancelCallback | None) -> None:
    if cancel_requested is not None and cancel_requested():
        raise RuntimeError("Version 2 pipeline cancelled")


def _update(progress: ProgressCallback | None, message: str, fraction: float) -> None:
    if progress:
        progress(message, max(0.0, min(1.0, float(fraction))))


def _number(value: Any, name: str, minimum: float, maximum: float) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be numeric") from error
    if not math.isfinite(result) or not minimum <= result <= maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return result


def normalize_backward_options(value: dict[str, Any] | None) -> dict[str, Any]:
    supplied = value or {}
    options = {**DEFAULT_BACKWARD_OPTIONS, **supplied}
    options["structural"] = {
        **DEFAULT_BACKWARD_OPTIONS["structural"],
        **supplied.get("structural", {}),
    }
    options["huri"] = {
        **DEFAULT_BACKWARD_OPTIONS["huri"],
        **supplied.get("huri", {}),
    }
    options["target"] = str(options["target"]).strip()
    options["receptor"] = str(options["receptor"]).strip()
    if not options["target"] or not options["receptor"]:
        raise ValueError("Version 2 requires a target and receptor")
    if options["target"].casefold() == options["receptor"].casefold():
        raise ValueError("target and receptor must be different")
    for key in ("cheap_top_n_per_frontier", "beam_width", "maximum_depth"):
        options[key] = int(_number(options[key], key, 1, 10_000))
    options["always_include_receptor"] = bool(options["always_include_receptor"])
    options["stop_when_receptor_reached"] = bool(options["stop_when_receptor_reached"])
    options["development_mode"] = bool(options.get("development_mode", True))
    structural = options["structural"]
    structural["enabled"] = bool(structural["enabled"])
    structural["reference_score"] = _number(
        structural["reference_score"], "structural reference score", 1e-6, 1 - 1e-6
    )
    structural["weight"] = _number(structural["weight"], "structural weight", 0.0, 10.0)
    structural["bayes_factor_floor"] = _number(
        structural["bayes_factor_floor"], "structural BF floor", 1e-12, 1.0
    )
    structural["bayes_factor_ceiling"] = _number(
        structural["bayes_factor_ceiling"], "structural BF ceiling", 1.0, 1e12
    )
    structural["num_cycle"] = int(_number(structural["num_cycle"], "num_cycle", 1, 20))
    structural["num_predictions_per_model"] = int(
        _number(structural["num_predictions_per_model"], "predictions per model", 1, 100)
    )
    if structural["num_predictions_per_model"] != 1:
        raise ValueError("Version 2 permits exactly one prediction per AlphaFold model")
    structural["model_names"] = str(structural.get("model_names", "")).strip()
    if not structural["model_names"] or "," in structural["model_names"]:
        raise ValueError("Version 2 requires exactly one AlphaFold model name")
    compatible = structural.get("compatible_protocol_ids", [])
    if isinstance(compatible, str):
        compatible = [item.strip() for item in compatible.split(",") if item.strip()]
    structural["compatible_protocol_ids"] = [
        str(item).strip() for item in compatible if str(item).strip()
    ]
    structural["maximum_concurrent_jobs"] = int(
        _number(structural["maximum_concurrent_jobs"], "concurrent GPU jobs", 1, 1000)
    )
    structural["metric"] = str(structural["metric"]).strip()
    structural["protocol_id"] = str(structural["protocol_id"]).strip()
    structural["legacy_screen_directory"] = str(
        structural.get("legacy_screen_directory", "")
    ).strip()
    structural["legacy_protocol_id"] = str(
        structural.get("legacy_protocol_id", LEGACY_PROTOCOL_ID)
    ).strip() or LEGACY_PROTOCOL_ID
    structural["require_legacy_ranked_structure"] = bool(
        structural.get("require_legacy_ranked_structure", True)
    )
    if not structural["metric"] or not structural["protocol_id"]:
        raise ValueError("structural metric and protocol_id are required")
    huri = options["huri"]
    huri["enabled"] = bool(huri.get("enabled", False))
    huri["interactions_file"] = str(huri.get("interactions_file", "")).strip()
    huri["mapping_file"] = str(huri.get("mapping_file", "")).strip()
    huri["screened_genes_file"] = str(huri.get("screened_genes_file", "")).strip()
    huri["weight"] = _number(huri.get("weight", 1.0), "HuRI weight", 0.0, 10.0)
    huri["positive_bayes_factor"] = _number(
        huri.get("positive_bayes_factor", 5.0), "HuRI positive BF", 1.0000001, 1e12
    )
    huri["nonreported_bayes_factor"] = _number(
        huri.get("nonreported_bayes_factor", 0.9), "HuRI non-report BF", 1e-12, 1.0
    )
    huri["nonreported_scope"] = str(huri.get("nonreported_scope", "neutral")).strip()
    if huri["nonreported_scope"] not in {"neutral", "huri_interactor_set", "screened_gene_file"}:
        raise ValueError("HuRI non-report scope is invalid")
    if huri["enabled"] and (not huri["interactions_file"] or not huri["mapping_file"]):
        raise ValueError("enabled HuRI evidence requires interactions and mapping files")
    if huri["enabled"] and huri["nonreported_scope"] == "screened_gene_file" and not huri["screened_genes_file"]:
        raise ValueError("HuRI screened-gene scope requires a screened genes file")
    return options


def default_pipeline_configuration(project_root: Path | str = PROJECT_DIRECTORY) -> dict[str, Any]:
    registry = load_registry(project_root)
    return {
        "bayesian": default_configuration(registry),
        "backward_search": json.loads(json.dumps(DEFAULT_BACKWARD_OPTIONS)),
    }


def _unique_graph_endpoints(
    selected_nodes: pd.DataFrame,
    node_factors: pd.DataFrame,
    project: Path,
    target: str,
    receptor: str,
) -> tuple[list[str], pd.DataFrame, dict[str, Any]]:
    universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    symbols = selected_nodes["symbol"].astype(str).tolist()
    metadata = selected_nodes.copy()
    target_symbol, metadata, target_resolution = append_requested_graph_endpoint(
        target, symbols, metadata, universe, node_factors, project,
        reason="v2_target_endpoint",
    )
    receptor_symbol, metadata, receptor_resolution = append_requested_graph_endpoint(
        receptor, symbols, metadata, universe, node_factors, project,
        reason="v2_receptor_endpoint",
    )
    if len(symbols) != len(set(symbols)):
        raise ValueError("endpoint resolution created duplicate graph symbols")
    return symbols, metadata, {
        "target_input": target,
        "target_symbol": target_symbol,
        "target_resolution": target_resolution,
        "receptor_input": receptor,
        "receptor_symbol": receptor_symbol,
        "receptor_resolution": receptor_resolution,
    }


def _ensure_incremental_evidence(
    project: Path,
    graph_metadata: pd.DataFrame,
    progress: ProgressCallback | None,
    cancel_requested: CancelCallback | None,
) -> dict[str, Any] | None:
    seed = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    seed_symbols = seed["symbol"].astype(str).tolist()
    if set(graph_metadata["symbol"].astype(str)).issubset(seed_symbols):
        return None
    update = ensure_incremental_pairs(
        project,
        graph_metadata,
        seed_symbols,
        progress=(lambda message, fraction: _update(progress, message, 0.22 + 0.20 * fraction)),
        cancel_check=(lambda: _check_cancel(cancel_requested)),
    )
    return update.as_dict()


def _export_applied_edge_streams(
    run_dir: Path,
    symbols: list[str],
    registry: dict[str, Any],
    bayesian: dict[str, Any],
    contributions: dict[str, np.ndarray],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    factor_dir = run_dir / "cheap_edge_factors"
    factor_dir.mkdir(parents=True, exist_ok=True)
    symbol_array = np.asarray(symbols, dtype=object)
    upper_left, upper_right = np.triu_indices(len(symbols), 1)
    definitions = {item["id"]: item for item in registry["edge_streams"]}
    search_streams: list[dict[str, Any]] = []
    manifest: list[dict[str, Any]] = []
    for stream_id, stream_state in bayesian["edge_streams"].items():
        weight = float(stream_state["weight"])
        if not stream_state["enabled"] or weight <= 0 or stream_id not in contributions:
            continue
        weighted_logs = contributions[stream_id][upper_left, upper_right]
        keep = np.abs(weighted_logs) > 1e-14
        factors = np.exp(np.clip(weighted_logs[keep] / weight, -700.0, 700.0))
        table = pd.DataFrame({
            "node_a": symbol_array[upper_left[keep]],
            "node_b": symbol_array[upper_right[keep]],
            "bayes_factor": factors,
        })
        path = factor_dir / f"{stream_id}.tsv.gz"
        table.to_csv(path, sep="\t", index=False, compression="gzip", float_format="%.12g")
        search_streams.append({
            "id": stream_id,
            "file": str(path.resolve()),
            "factor_column": "bayes_factor",
            "weight": weight,
            # Scope-aware penalties and continuous x=0 values were already
            # materialized by combine_edge_factors. Unlisted pairs are neutral.
            "missing_bayes_factor": 1.0,
        })
        definition = definitions[stream_id]
        manifest.append({
            "id": stream_id,
            "label": definition["label"],
            "weight": weight,
            "tq_or_reference_multiplier": stream_state.get("tq_multiplier"),
            "continuous_negative_evidence": stream_state.get("continuous_negative_evidence", False),
            "non_neutral_pair_count": int(len(table)),
            "factor_file": str(path.resolve()),
            "note": "Factors are unweighted here; the original configured weight is applied by the frontier engine.",
        })
    return search_streams, manifest


def _directionality_outputs(
    project: Path,
    run_dir: Path,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    registry: dict[str, Any],
    bayesian: dict[str, Any],
) -> tuple[Path | None, dict[str, Any] | None]:
    path_config = bayesian["path"]
    if not path_config["ontology_directionality_enabled"]:
        return None, None
    catalog = load_direction_rule_catalog(
        RULE_CATALOG_PATH,
        known_classes=[item["id"] for item in registry["path_ontology_classes"]],
    )
    omnipath = None
    if path_config["omnipath_directionality_enabled"]:
        raw_path = project / OMNIPATH_DIRECTION_RAW_RELATIVE
        raw = pd.read_csv(raw_path, sep="\t", dtype=str, keep_default_na=False)
        omnipath = workflow_build_omnipath_direction_evidence(
            raw, matrix.index.astype(str).tolist()
        )
    _, audit, summary = apply_ontology_directionality(
        matrix,
        metadata,
        catalog,
        audit_probability_cutoff=0.0,
        edge_output_cutoff=float(bayesian["edge_integration"]["output_probability_cutoff"]),
        path_probability_cutoff=0.0,
        omnipath_directions=omnipath,
    )
    audit_path = run_dir / "directionality_audit.tsv.gz"
    audit.to_csv(audit_path, sep="\t", index=False, compression="gzip")
    disallowed: list[dict[str, Any]] = []
    for row in audit.itertuples(index=False):
        if not bool(row.allowed_a_to_b):
            disallowed.append({"source": row.node_a, "target": row.node_b, "allowed": False})
        if not bool(row.allowed_b_to_a):
            disallowed.append({"source": row.node_b, "target": row.node_a, "allowed": False})
    path = run_dir / "disallowed_traversals.tsv"
    pd.DataFrame(disallowed, columns=["source", "target", "allowed"]).to_csv(
        path, sep="\t", index=False
    )
    summary["disallowed_traversal_count"] = len(disallowed)
    summary["audit_file"] = str(audit_path.resolve())
    summary["search_constraint_file"] = str(path.resolve())
    return path, summary


def _build_sequences(project: Path, run_dir: Path, metadata: pd.DataFrame) -> dict[str, Any]:
    node_table = run_dir / "graph_nodes.tsv"
    metadata.to_csv(node_table, sep="\t", index=False)
    fasta = run_dir / "node_sequences.fasta"
    summary = build_node_sequence_fasta(
        node_table,
        project / UNIPROT_RELATIVE,
        fasta,
    )
    mapped = set()
    if fasta.is_file():
        mapped = {
            line[1:].split()[0]
            for line in fasta.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        }
    node_types = (
        metadata["node_type"].astype(str)
        if "node_type" in metadata.columns
        else pd.Series("protein", index=metadata.index, dtype=str)
    )
    required = set(
        metadata.loc[node_types.str.casefold().eq("protein"), "symbol"].astype(str)
    )
    summary["protein_sequence_missing_count"] = len(required - mapped)
    summary["protein_sequence_missing_symbols"] = sorted(required - mapped)
    return summary


def initialize_end_to_end_pipeline(
    supplied_configuration: dict[str, Any] | None,
    run_dir: Path,
    *,
    project_root: Path | str = PROJECT_DIRECTORY,
    progress: ProgressCallback | None = None,
    cancel_requested: CancelCallback | None = None,
) -> dict[str, Any]:
    """Run configurable node selection and prepare a resumable V2 search."""
    project = Path(project_root).resolve()
    run_dir = run_dir.resolve()
    if run_dir.exists():
        raise FileExistsError(f"pipeline run already exists: {run_dir}")
    run_dir.mkdir(parents=True)
    supplied = supplied_configuration or {}
    registry = load_registry(project)
    bayesian = normalize_configuration(supplied.get("bayesian", supplied), registry)
    backward = normalize_backward_options(supplied.get("backward_search"))
    _write_json(run_dir / "submitted_configuration.json", supplied)
    _write_json(run_dir / "normalized_bayesian_configuration.json", bayesian)
    _write_json(run_dir / "normalized_backward_configuration.json", backward)
    try:
        _check_cancel(cancel_requested)
        calibration_summary: dict[str, Any] | None = None
        if bayesian["calibration"]["enabled"]:
            calibration_summary = {"enabled": True, "node": None, "edge": None}
            if bayesian["calibration"]["known_nodes"]:
                _update(progress, "Calibrating node evidence", 0.02)
                node_calibration, node_targets, node_trace = calibrate_node_parameters(
                    project,
                    registry,
                    bayesian,
                    progress=(lambda message, fraction: _update(progress, message, 0.02 + 0.03 * fraction)),
                    check_cancel=(lambda: _check_cancel(cancel_requested)),
                )
                calibration_summary["node"] = node_calibration
                node_targets.to_csv(run_dir / "node_calibration_targets.tsv", sep="\t", index=False)
                node_trace.to_csv(run_dir / "node_calibration_trace.tsv", sep="\t", index=False)
        _update(progress, "Selecting signaling nodes", 0.05)
        node_factors, selected_nodes, node_summary = select_nodes(project, registry, bayesian)
        node_factors.to_csv(
            run_dir / "node_posteriors.tsv.gz", sep="\t", index=False, compression="gzip"
        )
        selected_nodes.to_csv(run_dir / "selected_nodes.tsv", sep="\t", index=False)

        _check_cancel(cancel_requested)
        _update(progress, "Resolving receptor and target", 0.15)
        symbols, graph_metadata, endpoint_summary = _unique_graph_endpoints(
            selected_nodes,
            node_factors,
            project,
            backward["target"],
            backward["receptor"],
        )
        resolved_known_edges: list[tuple[str, str]] = []
        calibration_resolutions: list[dict[str, Any]] = []
        if calibration_summary and bayesian["calibration"]["known_edges"]:
            universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
            seen: set[tuple[str, str]] = set()
            for requested_left, requested_right in bayesian["calibration"]["known_edges"]:
                left, graph_metadata, left_resolution = append_requested_graph_endpoint(
                    requested_left, symbols, graph_metadata, universe, node_factors, project,
                    reason="v2_calibration_edge_endpoint",
                )
                right, graph_metadata, right_resolution = append_requested_graph_endpoint(
                    requested_right, symbols, graph_metadata, universe, node_factors, project,
                    reason="v2_calibration_edge_endpoint",
                )
                pair = tuple(sorted((left, right), key=lambda item: (item.casefold(), item)))
                if pair[0] == pair[1]:
                    raise ValueError(f"known calibration edge resolves to a self-edge: {pair[0]}")
                if pair not in seen:
                    resolved_known_edges.append((left, right))
                    seen.add(pair)
                calibration_resolutions.append({
                    "requested_node_a": requested_left,
                    "requested_node_b": requested_right,
                    "resolved_node_a": left,
                    "resolved_node_b": right,
                    "node_a_resolution": left_resolution,
                    "node_b_resolution": right_resolution,
                })
        graph_metadata.to_csv(run_dir / "graph_nodes.tsv", sep="\t", index=False)

        _check_cancel(cancel_requested)
        _update(progress, "Caching inexpensive evidence for newly selected nodes", 0.20)
        incremental_summary = _ensure_incremental_evidence(
            project, graph_metadata, progress, cancel_requested
        )

        if calibration_summary and resolved_known_edges:
            _update(progress, "Calibrating edge evidence", 0.43)
            edge_calibration, edge_targets, edge_trace = calibrate_edge_parameters(
                project,
                registry,
                bayesian,
                symbols,
                graph_metadata,
                resolved_known_edges,
                progress=(lambda message, fraction: _update(progress, message, 0.43 + 0.04 * fraction)),
                check_cancel=(lambda: _check_cancel(cancel_requested)),
            )
            edge_calibration["endpoint_resolution"] = calibration_resolutions
            calibration_summary["edge"] = edge_calibration
            edge_targets.to_csv(run_dir / "edge_calibration_targets.tsv", sep="\t", index=False)
            edge_trace.to_csv(run_dir / "edge_calibration_trace.tsv", sep="\t", index=False)
        if calibration_summary is not None:
            _write_json(run_dir / "calibration_summary.json", calibration_summary)
            _write_json(run_dir / "calibrated_bayesian_configuration.json", bayesian)

        _check_cancel(cancel_requested)
        _update(progress, "Integrating configurable inexpensive edge evidence", 0.48)
        contributions: dict[str, np.ndarray] = {}
        matrix, edge_summary = combine_edge_factors(
            project,
            registry,
            bayesian,
            symbols,
            graph_metadata=graph_metadata,
            contribution_collector=contributions,
        )
        matrix.to_csv(run_dir / "cheap_edge_posteriors.tsv.gz", sep="\t", compression="gzip", float_format="%.9g")
        streams, stream_manifest = _export_applied_edge_streams(
            run_dir, symbols, registry, bayesian, contributions
        )
        huri_summary: dict[str, Any] | None = None
        huri_config = backward["huri"]
        if huri_config["enabled"]:
            _update(progress, "Mapping HuRI binary interaction evidence", 0.64)
            interactions_path = Path(huri_config["interactions_file"])
            if not interactions_path.is_absolute():
                interactions_path = project / interactions_path
            mapping_path = Path(huri_config["mapping_file"])
            if not mapping_path.is_absolute():
                mapping_path = project / mapping_path
            screened_path = None
            if huri_config["screened_genes_file"]:
                screened_path = Path(huri_config["screened_genes_file"])
                if not screened_path.is_absolute():
                    screened_path = project / screened_path
            huri_stream, huri_summary = build_huri_search_evidence(
                interactions_path,
                mapping_path,
                graph_metadata,
                run_dir / "huri_evidence",
                positive_bayes_factor=huri_config["positive_bayes_factor"],
                nonreported_bayes_factor=huri_config["nonreported_bayes_factor"],
                nonreported_scope=huri_config["nonreported_scope"],
                screened_genes_file=screened_path,
            )
            huri_stream["weight"] = float(huri_config["weight"])
            streams.append(huri_stream)
            stream_manifest.append({
                "id": huri_stream["id"],
                "label": "HuRI binary interaction",
                "weight": float(huri_config["weight"]),
                "non_neutral_pair_count": huri_summary["reported_positive_graph_pair_count"],
                "factor_file": huri_stream["file"],
                "scope_nodes_file": huri_stream["scope_nodes_file"],
                "scoped_missing_bayes_factor": huri_stream["scoped_missing_bayes_factor"],
                "note": huri_summary["scope_note"],
            })
        _write_json(run_dir / "cheap_edge_stream_manifest.json", stream_manifest)

        _check_cancel(cancel_requested)
        _update(progress, "Preparing directionality constraints", 0.68)
        direction_file, direction_summary = _directionality_outputs(
            project, run_dir, matrix, graph_metadata, registry, bayesian
        )

        _check_cancel(cancel_requested)
        _update(progress, "Mapping selected proteins to AlphaPulldown sequences", 0.75)
        sequence_summary = _build_sequences(project, run_dir, graph_metadata)
        target_symbol = endpoint_summary["target_symbol"]
        receptor_symbol = endpoint_summary["receptor_symbol"]
        node_types = graph_metadata.set_index("symbol").get("node_type", pd.Series(dtype=str))
        protein_endpoints = {
            symbol for symbol in (target_symbol, receptor_symbol)
            if str(node_types.get(symbol, "protein")).casefold() == "protein"
        }
        missing_endpoint_sequences = protein_endpoints.intersection(
            sequence_summary["protein_sequence_missing_symbols"]
        )
        if missing_endpoint_sequences:
            raise ValueError(
                "Protein endpoint sequence mapping failed for: "
                + ", ".join(sorted(missing_endpoint_sequences))
            )

        search_configuration = {
            "project_root": str(project),
            "target": target_symbol,
            "receptor": receptor_symbol,
            "node_table": str((run_dir / "graph_nodes.tsv").resolve()),
            "sequence_fasta": str((run_dir / "node_sequences.fasta").resolve()),
            "symbol_column": "symbol",
            "classes_column": "classes",
            "node_type_column": "node_type",
            "allowed_intermediate_classes": list(
                bayesian["path"]["allowed_intermediate_classes"]
            ),
            "prior_probability": float(bayesian["edge_integration"]["prior_probability"]),
            "cheap_top_n_per_frontier": backward["cheap_top_n_per_frontier"],
            "beam_width": backward["beam_width"],
            "maximum_depth": backward["maximum_depth"],
            "always_include_receptor": backward["always_include_receptor"],
            "stop_when_receptor_reached": backward["stop_when_receptor_reached"],
            "development_mode": backward["development_mode"],
            "pair_cache": str((PROJECT_DIRECTORY / "runtime" / "pair_cache.sqlite3").resolve()),
            "cheap_streams": streams,
            "structural": backward["structural"],
        }
        if direction_file is not None:
            search_configuration.update({
                "directionality_file": str(direction_file.resolve()),
                "directionality_mode": "disallowed_rows",
            })
        search_config_path = run_dir / "backward_search_configuration.json"
        _write_json(search_config_path, search_configuration)
        _check_cancel(cancel_requested)
        _update(progress, "Initializing the stepwise backward frontier", 0.88)
        search_run = run_dir / "backward_search"
        search_state = initialize_run(
            search_config_path,
            search_run,
            development_mode=backward["development_mode"],
        )
        legacy_summary: dict[str, Any] | None = None
        legacy_directory_value = backward["structural"].get("legacy_screen_directory", "")
        if legacy_directory_value:
            legacy_directory = Path(legacy_directory_value)
            if not legacy_directory.is_absolute():
                legacy_directory = project / legacy_directory
            if legacy_directory.joinpath("scores").is_dir():
                _update(progress, "Importing valid legacy ppi_screen results", 0.94)
                legacy_summary = import_legacy_ppi_screen(
                    legacy_directory,
                    Path(search_state["pair_cache"]),
                    protocol_id=backward["structural"]["legacy_protocol_id"],
                    metric=backward["structural"]["metric"],
                    require_ranked_structure=backward["structural"]["require_legacy_ranked_structure"],
                    audit_directory=run_dir / "legacy_ppi_screen_import",
                )
        pipeline_state = {
            "schema_version": PIPELINE_SCHEMA_VERSION,
            "created_at": utc_now(),
            "status": search_state["status"],
            "stage": "backward_search",
            "run_directory": str(run_dir),
            "project_root": str(project),
            "code_project_root": str(PROJECT_DIRECTORY),
            "search_run_directory": str(search_run),
            "selected_node_count": int(node_summary["selected_node_count"]),
            "selected_protein_count": int(node_summary["selected_protein_count"]),
            "node_selection_summary": node_summary,
            "cheap_edge_summary": edge_summary,
            "incremental_edge_cache": incremental_summary,
            "endpoint_summary": endpoint_summary,
            "sequence_summary": sequence_summary,
            "directionality_summary": direction_summary,
            "calibration_summary": calibration_summary,
            "cheap_edge_streams": stream_manifest,
            "huri_evidence": huri_summary,
            "search_status": search_state["status"],
            "search_round": search_state["round"],
            "next_search_stage": search_state.get("development_stage", ""),
            "pending_unique_structural_pairs": search_state.get("pending_unique_structural_pairs", 0),
            "solution_count": len(search_state.get("solutions", [])),
            "structural_submission": None,
            "legacy_ppi_screen_import": legacy_summary,
        }
        _write_pipeline_state(run_dir, pipeline_state)
        _update(progress, "Complete pipeline initialized", 1.0)
        return pipeline_state
    except Exception as error:
        failed = {
            "schema_version": PIPELINE_SCHEMA_VERSION,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "status": "failed",
            "stage": "initialization",
            "run_directory": str(run_dir),
            "project_root": str(project),
            "error": str(error),
        }
        _write_json(_pipeline_path(run_dir), failed)
        raise


def refresh_pipeline_state(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    state = load_pipeline_state(run_dir)
    search_dir_value = state.get("search_run_directory")
    if search_dir_value and Path(search_dir_value).joinpath("state.json").is_file():
        search = load_run_state(Path(search_dir_value))
        submitted = state.get("structural_submission")
        submitted_for_round = bool(
            submitted
            and int(submitted.get("round", -1)) == int(search["round"])
            and search["status"] == "waiting_for_structural"
        )
        state.update({
            "status": "structural_submitted" if submitted_for_round else search["status"],
            "search_status": search["status"],
            "search_round": search["round"],
            "next_search_stage": search.get("development_stage", ""),
            "pending_unique_structural_pairs": search.get("pending_unique_structural_pairs", 0),
            "solution_count": len(search.get("solutions", [])),
            "completion_reason": search.get("completion_reason", ""),
            "development_trace": str(Path(search_dir_value) / "development_trace.html"),
        })
        _write_pipeline_state(run_dir, state)
    return state


def step_pipeline(run_dir: Path, *, until_checkpoint: bool = False) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    pipeline = load_pipeline_state(run_dir)
    search_dir = Path(pipeline["search_run_directory"])
    maximum_steps = 100 if until_checkpoint else 1
    for _ in range(maximum_steps):
        before = load_run_state(search_dir)
        if before["status"] in {"complete", "waiting_for_structural"}:
            break
        after = step_run(search_dir)
        if after["status"] in {"complete", "waiting_for_structural"}:
            break
    return refresh_pipeline_state(run_dir)


def collect_pipeline_structural_scores(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    pipeline = load_pipeline_state(run_dir)
    search_dir = Path(pipeline["search_run_directory"])
    score_path = collect_alphapulldown_scores(search_dir)
    search_state = import_structural_scores(search_dir, score_path)
    record_structural_import(search_dir, search_state)
    pipeline = refresh_pipeline_state(run_dir)
    pipeline["last_structural_collection"] = {
        "at": utc_now(),
        "score_file": str(score_path),
    }
    _write_pipeline_state(run_dir, pipeline)
    return pipeline


def submit_pipeline_structural_round(run_dir: Path) -> dict[str, Any]:
    """Submit the prepared round when running inside a Slurm environment."""
    run_dir = run_dir.resolve()
    pipeline = refresh_pipeline_state(run_dir)
    if pipeline["search_status"] != "waiting_for_structural":
        raise ValueError("pipeline is not waiting for structural predictions")
    if os.name == "nt" or shutil.which("sbatch") is None or shutil.which("bash") is None:
        raise RuntimeError(
            "Slurm submission is available only when this GUI is running on Biowulf; "
            "the prepared round can still be copied and submitted there."
        )
    project = Path(pipeline.get("code_project_root", PROJECT_DIRECTORY))
    round_dir = Path(pipeline["search_run_directory"]) / f"round_{int(pipeline['search_round']):03d}"
    command = ["bash", str(project / "biowulf/submit_round.sh"), str(round_dir)]
    environment = os.environ.copy()
    structural_config = load_run_state(Path(pipeline["search_run_directory"]))["configuration"]["structural"]
    environment["GBI_MAX_CONCURRENT_GPU_JOBS"] = str(
        structural_config["maximum_concurrent_jobs"]
    )
    result = subprocess.run(
        command,
        cwd=project,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or "Slurm submission failed")
    submission_path = round_dir / "submission.tsv"
    jobs: list[dict[str, str]] = []
    if submission_path.is_file():
        jobs = pd.read_csv(submission_path, sep="\t", dtype=str).fillna("").to_dict("records")
    pipeline["structural_submission"] = {
        "at": utc_now(),
        "round": int(pipeline["search_round"]),
        "jobs": jobs,
        "stdout": result.stdout.strip(),
        "round_directory": str(round_dir),
    }
    pipeline["status"] = "structural_submitted"
    _write_pipeline_state(run_dir, pipeline)
    return pipeline


def render_pipeline_trace(run_dir: Path) -> Path:
    state = refresh_pipeline_state(run_dir)
    return render_development_report(Path(state["search_run_directory"]))
