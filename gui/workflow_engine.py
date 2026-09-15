#!/usr/bin/env python3
"""Configurable end-to-end Bayesian signaling workflow used by the GUI."""

from __future__ import annotations

import gzip
import copy
import json
import math
import re
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd


GUI_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = GUI_DIR.parent
CODE_DIR = PROJECT_ROOT / "code"
PATH_CODE_DIR = CODE_DIR / "path_finding"
EDGE_CODE_DIR = CODE_DIR / "edge_characterization"
for search_path in (CODE_DIR, PATH_CODE_DIR, EDGE_CODE_DIR):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from build_target_adjacency_vector import (  # noqa: E402
    build_target_adjacency_vector,
    resolve_target,
    safe_name,
)
from find_ranked_paths import (  # noqa: E402
    DEFAULT_SIGNAL_RELAY_CLASSES,
    build_adjacency,
    build_intermediate_eligibility,
    cached_target_extension,
    connected_component_size,
    k_shortest_simple_paths,
    resolve_existing_node,
    serialize_paths,
)
from incremental_edge_cache import (  # noqa: E402
    ensure_incremental_pairs,
    incremental_factor_table,
    incremental_pair_factor_table,
    incremental_node_scope_table,
)
from ontology_directionality import (  # noqa: E402
    RULE_CATALOG_PATH,
    apply_ontology_directionality,
    build_omnipath_direction_evidence,
    build_complete_class_pair_catalog,
    load_direction_rule_catalog,
)
from temporal_path_ranking import (  # noqa: E402
    P_ADJUST_METHODS,
    build_replicate_responses,
    load_site_trajectories,
    rank_paths_replicate,
    replicate_response_table,
    temporal_validation_summary,
    variance_trend_table,
)
from parameter_calibration import (  # noqa: E402
    bounded_powell_positive_calibration,
)


ProgressCallback = Callable[[str, float], None]
CancellationCallback = Callable[[], bool]


class WorkflowCancelled(RuntimeError):
    """Raised when a user requests cooperative workflow cancellation."""


NODE_FACTORS_RELATIVE = Path(
    "results/combined_kinase_phosphoprotein_evidence/combined_all_nodes.tsv"
)
UNIVERSE_RELATIVE = Path("data/node_selection/node_universe_combined_nonzero.tsv")
UNIPROT_RELATIVE = Path(
    "data/edge_characterization/kinase_predictor/"
    "phosphosite_database/raw/uniprot_mouse_reference_proteome.tsv.gz"
)
RUNS_RELATIVE = Path("results/gui_runs")
TEMPORAL_PHOSPHO_RELATIVE = Path("data/phospho_data_original.xlsx")
PROTEIN_RAW_RELATIVE = Path("results/mpkccd_protein_abundance_bayes_factors.tsv")
PC_RAW_RELATIVE = Path("results/pc_median_tpm_bayes_factors.tsv")
PHOSPHOPROTEIN_RAW_RELATIVE = Path(
    "results/phosphoprotein_evidence/node_selection_protein_pc_phosphosite_posterior.tsv"
)
MPKCCD_PROFILES_RELATIVE = Path(
    "data/edge_characterization/localization/processed/node_localization_profiles.tsv"
)
KINASE_PREDICTIONS_RELATIVE = Path(
    "results/edge_characterization/localization_kinase_predictor/"
    "observed_phosphosite_top10_predictions.tsv.gz"
)
STRING_SCORE_MATRIX_RELATIVE = Path(
    "results/edge_characterization/localization_kinase_predictor_string/"
    "string_combined_score_matrix.tsv"
)
HPA_PROFILES_RELATIVE = Path(
    "data/edge_characterization/localization/hpa/v25.1/processed/"
    "node_hpa_localization_profiles.tsv"
)
OMNIPATH_EFFORT_MATRIX_RELATIVE = Path(
    "results/edge_characterization/localization_kinase_predictor_string_hpa_omnipath/"
    "omnipath_curation_effort_matrix.tsv"
)
OMNIPATH_DIRECTION_RAW_RELATIVE = Path(
    "data/edge_characterization/omnipath/2026-07-30/raw/"
    "omnipath_mouse_core_post_translational.tsv"
)
STITCH_UNIVERSE_EDGES_RELATIVE = Path(
    "results/edge_characterization/"
    "localization_kinase_predictor_string_hpa_omnipath_stitch/"
    "stitch_secondary_messenger_edges.tsv"
)
STITCH_ALL_MOUSE_EDGES_RELATIVE = Path(
    "data/edge_characterization/stitch/v5.0/processed/"
    "stitch_secondary_messenger_all_mouse_edges.tsv.gz"
)
STRING_MAPPING_RELATIVE = Path(
    "data/edge_characterization/string/v12.0/processed/node_to_string_mapping.tsv"
)
OBSERVED_PHOSPHOSITES_RELATIVE = Path(
    "data/edge_characterization/kinase_predictor/phosphosite_database/"
    "observed_phosphosites.tsv"
)
STRING_REFERENCE_SCORE = 0.041
STITCH_REFERENCE_SCORE = 0.150
NEUTRAL_LIKELIHOOD = 0.5
FACTOR_EPSILON = 1e-12
PHOSPHOPROTEOMIC_CALIBRATION_STREAMS = {
    "kinase_activity",
    "phosphoprotein_response",
    "pka_ca_ko_phosphoprotein_response",
    "pka_cb_ko_phosphoprotein_response",
    "kinase_predictor",
}
CALIBRATION_WEIGHT_BOUNDS = (0.0, 3.0)
CALIBRATION_TQ_BOUNDS = (0.25, 4.0)
CALIBRATION_PHOSPHO_TQ_BOUNDS = (0.05, 1.0)


def _calibration_tq_bounds(definition: dict[str, Any]) -> tuple[float, float]:
    """Return registry-defined calibration bounds for one scale multiplier."""
    configured = definition.get("normalization", {}).get("calibration_bounds")
    if configured is None:
        return (
            CALIBRATION_PHOSPHO_TQ_BOUNDS
            if definition["id"] in PHOSPHOPROTEOMIC_CALIBRATION_STREAMS
            else CALIBRATION_TQ_BOUNDS
        )
    if not isinstance(configured, list) or len(configured) != 2:
        raise ValueError(
            f"node/edge stream {definition['id']} has invalid calibration bounds"
        )
    lower, upper = map(float, configured)
    if not (math.isfinite(lower) and math.isfinite(upper) and 0 < lower < upper):
        raise ValueError(
            f"node/edge stream {definition['id']} has invalid calibration bounds"
        )
    return lower, upper


def _default_preferred_tq_multiplier(definition: dict[str, Any]) -> float:
    """Return the scientist-editable regularization anchor for one stream."""
    normalization = definition.get("normalization", {})
    if "calibration_preferred_multiplier" in normalization:
        return float(normalization["calibration_preferred_multiplier"])
    return (
        0.1
        if definition["id"] in PHOSPHOPROTEOMIC_CALIBRATION_STREAMS
        else 1.0
    )


@dataclass
class WorkflowResult:
    run_id: str
    output_directory: Path
    summary: dict[str, Any]
    preview: dict[str, Any]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_registry(project_root: Path | str = PROJECT_ROOT) -> dict[str, Any]:
    project = Path(project_root).resolve()
    registry_path = project / "gui/evidence_registry.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    validate_registry(registry, project)
    return registry


def validate_registry(registry: dict[str, Any], project: Path) -> None:
    if registry.get("schema_version") != 1:
        raise ValueError("unsupported evidence-registry schema")
    for group_name in ("node_streams", "edge_streams"):
        streams = registry.get(group_name, [])
        ids = [stream.get("id") for stream in streams]
        if not streams or None in ids or len(ids) != len(set(ids)):
            raise ValueError(f"invalid or duplicate IDs in {group_name}")
    for stream in registry["node_streams"]:
        factor_file = stream.get("factor_file")
        if factor_file and not (project / factor_file).exists():
            raise FileNotFoundError(
                f"node factor table is missing: {project / factor_file}"
            )
    for stream in registry["edge_streams"]:
        if stream.get("derived"):
            continue
        factor_path = project / stream["factor_file"]
        if not factor_path.exists():
            raise FileNotFoundError(f"edge factor table is missing: {factor_path}")
    ontology_classes = registry.get("path_ontology_classes", [])
    ontology_ids = [item.get("id") for item in ontology_classes]
    if (
        not ontology_classes
        or None in ontology_ids
        or len(ontology_ids) != len(set(ontology_ids))
    ):
        raise ValueError("invalid or duplicate IDs in path_ontology_classes")
    missing_defaults = set(DEFAULT_SIGNAL_RELAY_CLASSES).difference(ontology_ids)
    if missing_defaults:
        raise ValueError(
            "path ontology registry is missing default relay classes: "
            + ", ".join(sorted(missing_defaults))
        )
    enabled_by_default = {
        item["id"] for item in ontology_classes if item.get("default_enabled")
    }
    if enabled_by_default != set(DEFAULT_SIGNAL_RELAY_CLASSES):
        raise ValueError(
            "default-enabled path ontology classes must match the audited "
            "DEFAULT_SIGNAL_RELAY_CLASSES policy"
        )
    load_direction_rule_catalog(
        RULE_CATALOG_PATH,
        known_classes=ontology_ids,
    )


def default_configuration(registry: dict[str, Any]) -> dict[str, Any]:
    defaults = registry["defaults"]
    return {
        "calibration": {
            "enabled": False,
            "known_nodes": [],
            "known_edges": [],
            "regularization_strength": 0.1,
            "multistart_count": 2,
        },
        "node_streams": {
            stream["id"]: {
                "enabled": bool(stream["default_enabled"]),
                "weight": float(stream["default_weight"]),
                "continuous_negative_evidence": False,
                **(
                    {
                        "tq_multiplier": float(
                            stream.get("normalization", {}).get(
                                "default_multiplier", 1.0
                            )
                        ),
                        "preferred_tq_multiplier": (
                            _default_preferred_tq_multiplier(stream)
                        ),
                    }
                    if stream.get("normalization", {}).get("user_control", True)
                    else {}
                ),
            }
            for stream in registry["node_streams"]
        },
        "node_integration": {
            "include_second_messengers": bool(defaults["include_second_messengers"]),
            "non_neutral_tolerance": float(defaults["node_non_neutral_tolerance"]),
            "prior_probability": float(defaults["node_prior_probability"]),
            "output_probability_cutoff": float(defaults["node_probability_cutoff"]),
            "penalize_unobserved": bool(defaults["penalize_unobserved_nodes"]),
            "unobserved_bayes_factor": float(
                defaults["unobserved_node_bayes_factor"]
            ),
            "continuous_negative_evidence": bool(
                defaults["continuous_negative_node_evidence"]
            ),
            "continuous_bayes_factor_floor": float(
                defaults["continuous_node_bayes_factor_floor"]
            ),
        },
        "edge_streams": {
            stream["id"]: {
                "enabled": bool(stream["default_enabled"]),
                "weight": float(stream["default_weight"]),
                **(
                    {"continuous_negative_evidence": False}
                    if not stream.get("derived")
                    else {}
                ),
                **(
                    {
                        "tq_multiplier": float(
                            stream.get("normalization", {}).get(
                                "default_multiplier", 1.0
                            )
                        ),
                        "preferred_tq_multiplier": (
                            _default_preferred_tq_multiplier(stream)
                        ),
                    }
                    if stream.get("normalization", {}).get("user_control", True)
                    else {}
                ),
                "parameters": {
                    parameter["id"]: float(parameter["default"])
                    for parameter in stream.get("parameters", [])
                },
            }
            for stream in registry["edge_streams"]
        },
        "edge_integration": {
            "prior_probability": float(defaults["edge_prior_probability"]),
            "output_probability_cutoff": float(defaults["edge_probability_cutoff"]),
            "penalize_unsupported": bool(defaults["penalize_unsupported_edges"]),
            "unsupported_bayes_factor": float(
                defaults["unsupported_edge_bayes_factor"]
            ),
            "continuous_negative_evidence": bool(
                defaults["continuous_negative_edge_evidence"]
            ),
            "continuous_bayes_factor_floor": float(
                defaults["continuous_edge_bayes_factor_floor"]
            ),
        },
        "path": {
            "enabled": True,
            "start": defaults["start_node"],
            "target": defaults["target_node"],
            "top_k": int(defaults["top_k_paths"]),
            "max_hops": int(defaults["maximum_hops"]),
            "minimum_edge_probability": float(
                defaults["path_minimum_edge_probability"]
            ),
            "ontology_directionality_enabled": bool(
                defaults["ontology_directionality_enabled"]
            ),
            "omnipath_directionality_enabled": bool(
                defaults["omnipath_directionality_enabled"]
            ),
            "signaling_intermediates_only": bool(
                defaults["signaling_intermediates_only"]
            ),
            "exclude_multirole_scaffolds": bool(
                defaults["exclude_multirole_scaffolds"]
            ),
            "allowed_intermediate_classes": [
                item["id"]
                for item in registry["path_ontology_classes"]
                if item.get("default_enabled")
            ],
        },
        "temporal_validation": {
            "enabled": bool(defaults["temporal_validation_enabled"]),
            "prior_df": float(defaults["temporal_prior_df"]),
            "alpha": float(defaults["temporal_alpha"]),
            "monte_carlo_draws": int(defaults["temporal_monte_carlo_draws"]),
            "random_seed": int(defaults["temporal_random_seed"]),
            "p_adjust_method": str(defaults["temporal_p_adjust_method"]),
            "minimum_scored_nodes": int(defaults["temporal_minimum_scored_nodes"]),
        },
    }


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return number


def normalize_configuration(
    supplied: dict[str, Any] | None,
    registry: dict[str, Any],
) -> dict[str, Any]:
    config = default_configuration(registry)
    supplied = supplied or {}
    supplied_calibration = supplied.get("calibration", {})
    if not isinstance(supplied_calibration, dict):
        raise ValueError("calibration configuration must be an object")
    legacy_phospho_preference: float | None = None
    if "phosphoproteomic_preferred_tq_multiplier" in supplied_calibration:
        legacy_phospho_preference = _number(
            supplied_calibration["phosphoproteomic_preferred_tq_multiplier"],
            "legacy phosphoproteomic preferred Tq multiplier",
            CALIBRATION_PHOSPHO_TQ_BOUNDS[0],
            CALIBRATION_PHOSPHO_TQ_BOUNDS[1],
        )
    for group in ("node_streams", "edge_streams"):
        incoming = supplied.get(group, {})
        definitions = {item["id"]: item for item in registry[group]}
        for stream_id, state in config[group].items():
            definition = definitions[stream_id]
            if stream_id in incoming:
                state["enabled"] = bool(incoming[stream_id].get("enabled", state["enabled"]))
                state["weight"] = _number(
                    incoming[stream_id].get("weight", state["weight"]),
                    f"{stream_id} weight",
                    0.0,
                    10.0,
                )
                if "continuous_negative_evidence" in state:
                    state["continuous_negative_evidence"] = bool(
                        incoming[stream_id].get(
                            "continuous_negative_evidence",
                            state["continuous_negative_evidence"],
                        )
                    )
                if definition.get("normalization", {}).get("user_control", True):
                    state["tq_multiplier"] = _number(
                        incoming[stream_id].get(
                            "tq_multiplier", state["tq_multiplier"]
                        ),
                        f"{stream_id} Tq multiplier",
                        0.05,
                        20.0,
                    )
                parameter_values = incoming[stream_id].get("parameters", {})
                for parameter in definition.get("parameters", []):
                    parameter_id = parameter["id"]
                    state["parameters"][parameter_id] = _number(
                        parameter_values.get(
                            parameter_id, state["parameters"][parameter_id]
                        ),
                        f"{stream_id} {parameter_id}",
                        float(parameter["minimum"]),
                        float(parameter["maximum"]),
                    )
            if definition.get("normalization", {}).get("user_control", True):
                preferred_value = state["preferred_tq_multiplier"]
                if stream_id in incoming:
                    preferred_value = incoming[stream_id].get(
                        "preferred_tq_multiplier", preferred_value
                    )
                if (
                    legacy_phospho_preference is not None
                    and stream_id in PHOSPHOPROTEOMIC_CALIBRATION_STREAMS
                    and (
                        stream_id not in incoming
                        or "preferred_tq_multiplier" not in incoming[stream_id]
                    )
                ):
                    preferred_value = legacy_phospho_preference
                preferred_bounds = _calibration_tq_bounds(definition)
                state["preferred_tq_multiplier"] = _number(
                    preferred_value,
                    f"{stream_id} preferred Tq/reference multiplier",
                    preferred_bounds[0],
                    preferred_bounds[1],
                )
    config["node_integration"].update(supplied.get("node_integration", {}))
    config["edge_integration"].update(supplied.get("edge_integration", {}))
    config["path"].update(supplied.get("path", {}))
    config["temporal_validation"].update(supplied.get("temporal_validation", {}))

    calibration = config["calibration"]
    calibration.update(
        {
            key: value
            for key, value in supplied_calibration.items()
            if key != "phosphoproteomic_preferred_tq_multiplier"
        }
    )
    calibration["enabled"] = bool(calibration["enabled"])
    calibration["regularization_strength"] = _number(
        calibration["regularization_strength"],
        "calibration regularization strength",
        0.0,
        100.0,
    )
    calibration["multistart_count"] = int(
        _number(
            calibration["multistart_count"],
            "calibration multistart count",
            1,
            3,
        )
    )
    raw_nodes = calibration.get("known_nodes", [])
    if isinstance(raw_nodes, str):
        raw_nodes = re.split(r"[\s,;]+", raw_nodes)
    if not isinstance(raw_nodes, list):
        raise ValueError("known calibration nodes must be a list")
    normalized_nodes: list[str] = []
    seen_nodes: set[str] = set()
    for value in raw_nodes:
        symbol = str(value).strip()
        if symbol and symbol.casefold() not in seen_nodes:
            normalized_nodes.append(symbol)
            seen_nodes.add(symbol.casefold())
    calibration["known_nodes"] = normalized_nodes
    raw_edges = calibration.get("known_edges", [])
    if isinstance(raw_edges, str):
        raw_edges = [line for line in raw_edges.splitlines() if line.strip()]
    if not isinstance(raw_edges, list):
        raise ValueError("known calibration edges must be a list")
    normalized_edges: list[list[str]] = []
    seen_edges: set[tuple[str, str]] = set()
    for item in raw_edges:
        if isinstance(item, str):
            fields = [field.strip() for field in re.split(r"[\t,;]+", item) if field.strip()]
        elif isinstance(item, (list, tuple)):
            fields = [str(field).strip() for field in item if str(field).strip()]
        else:
            raise ValueError("each known edge must be a two-node list or delimited line")
        if len(fields) != 2:
            raise ValueError(
                "each known edge must contain exactly two nodes, for example Prkaca,Aqp2"
            )
        if fields[0].casefold() == fields[1].casefold():
            raise ValueError("known calibration edges cannot be self-edges")
        key = tuple(sorted((fields[0].casefold(), fields[1].casefold())))
        if key not in seen_edges:
            normalized_edges.append(fields)
            seen_edges.add(key)
    calibration["known_edges"] = normalized_edges
    if calibration["enabled"] and not (
        calibration["known_nodes"] or calibration["known_edges"]
    ):
        raise ValueError(
            "calibration is enabled but no known nodes or known edges were supplied"
        )

    node_integration = config["node_integration"]
    node_integration["include_second_messengers"] = bool(
        node_integration["include_second_messengers"]
    )
    node_integration["penalize_unobserved"] = bool(
        node_integration["penalize_unobserved"]
    )
    node_integration["continuous_negative_evidence"] = bool(
        node_integration["continuous_negative_evidence"]
    )
    if node_integration["continuous_negative_evidence"]:
        # The continuous mode supersedes the older fixed nondetection rule.
        node_integration["penalize_unobserved"] = False
    node_integration["unobserved_bayes_factor"] = _number(
        node_integration["unobserved_bayes_factor"],
        "unobserved node Bayes factor",
        1e-6,
        1.0,
    )
    node_integration["continuous_bayes_factor_floor"] = _number(
        node_integration["continuous_bayes_factor_floor"],
        "continuous node Bayes-factor floor",
        1e-12,
        1.0,
    )
    node_integration["non_neutral_tolerance"] = _number(
        node_integration["non_neutral_tolerance"],
        "node non-neutral tolerance",
        0.0,
        0.1,
    )
    node_integration["prior_probability"] = _number(
        node_integration["prior_probability"],
        "node prior probability",
        1e-9,
        1 - 1e-9,
    )
    node_integration["output_probability_cutoff"] = _number(
        node_integration["output_probability_cutoff"],
        "node output cutoff",
        0.0,
        1.0,
    )
    edge_integration = config["edge_integration"]
    edge_integration["penalize_unsupported"] = bool(
        edge_integration["penalize_unsupported"]
    )
    edge_integration["continuous_negative_evidence"] = bool(
        edge_integration["continuous_negative_evidence"]
    )
    if edge_integration["continuous_negative_evidence"]:
        edge_integration["penalize_unsupported"] = False
    edge_integration["unsupported_bayes_factor"] = _number(
        edge_integration["unsupported_bayes_factor"],
        "unsupported edge Bayes factor",
        1e-6,
        1.0,
    )
    edge_integration["continuous_bayes_factor_floor"] = _number(
        edge_integration["continuous_bayes_factor_floor"],
        "continuous edge Bayes-factor floor",
        1e-12,
        1.0,
    )
    edge_integration["prior_probability"] = _number(
        edge_integration["prior_probability"], "edge prior probability", 1e-9, 1 - 1e-9
    )
    edge_integration["output_probability_cutoff"] = _number(
        edge_integration["output_probability_cutoff"],
        "edge output cutoff",
        0.0,
        1.0,
    )
    path = config["path"]
    path["enabled"] = bool(path["enabled"])
    path["ontology_directionality_enabled"] = bool(
        path["ontology_directionality_enabled"]
    )
    path["omnipath_directionality_enabled"] = bool(
        path["omnipath_directionality_enabled"]
    )
    path["start"] = str(path["start"]).strip()
    path["target"] = str(path["target"]).strip()
    path["top_k"] = int(_number(path["top_k"], "top paths", 1, 500))
    path["max_hops"] = int(_number(path["max_hops"], "maximum hops", 1, 12))
    path["minimum_edge_probability"] = _number(
        path["minimum_edge_probability"], "path edge cutoff", 0.0, 1 - 1e-12
    )
    path["signaling_intermediates_only"] = bool(path["signaling_intermediates_only"])
    path["exclude_multirole_scaffolds"] = bool(path["exclude_multirole_scaffolds"])
    selected_classes = path.get("allowed_intermediate_classes")
    if not isinstance(selected_classes, list):
        raise ValueError("allowed intermediate ontology classes must be a list")
    selected_class_ids = {
        str(value).strip() for value in selected_classes if str(value).strip()
    }
    known_class_ids = {
        item["id"] for item in registry["path_ontology_classes"]
    }
    unknown_class_ids = selected_class_ids.difference(known_class_ids)
    if unknown_class_ids:
        raise ValueError(
            "unknown intermediate ontology classes: "
            + ", ".join(sorted(unknown_class_ids))
        )
    path["allowed_intermediate_classes"] = [
        item["id"]
        for item in registry["path_ontology_classes"]
        if item["id"] in selected_class_ids
    ]
    if path["enabled"] and (not path["start"] or not path["target"]):
        raise ValueError("both a starting node and target are required for path finding")

    temporal = config["temporal_validation"]
    temporal["enabled"] = bool(temporal["enabled"])
    temporal["prior_df"] = _number(
        temporal["prior_df"], "temporal prior degrees of freedom", 0.1, 1000.0
    )
    temporal["alpha"] = _number(
        temporal["alpha"], "temporal significance alpha", 1e-6, 1.0
    )
    temporal["monte_carlo_draws"] = int(
        _number(
            temporal["monte_carlo_draws"],
            "temporal Monte-Carlo draws",
            100,
            100000,
        )
    )
    temporal["random_seed"] = int(
        _number(temporal["random_seed"], "temporal random seed", 0, 2**32 - 1)
    )
    temporal["p_adjust_method"] = str(temporal["p_adjust_method"]).strip()
    if temporal["p_adjust_method"] not in P_ADJUST_METHODS:
        raise ValueError(
            "temporal p-adjustment must be one of: " + ", ".join(P_ADJUST_METHODS)
        )
    temporal["minimum_scored_nodes"] = int(
        _number(
            temporal["minimum_scored_nodes"],
            "minimum temporally scored path nodes",
            2,
            20,
        )
    )
    if temporal["enabled"] and not path["enabled"]:
        raise ValueError("temporal path validation requires path finding to be enabled")

    effective_nodes = [
        stream_id
        for stream_id, state in config["node_streams"].items()
        if state["enabled"] and state["weight"] > 0
    ]
    effective_edges = [
        stream_id
        for stream_id, state in config["edge_streams"].items()
        if state["enabled"] and state["weight"] > 0
    ]
    if not effective_nodes:
        raise ValueError("enable at least one node evidence stream with positive weight")
    if not effective_edges:
        raise ValueError("enable at least one edge evidence stream with positive weight")

    edge_defs = {stream["id"]: stream for stream in registry["edge_streams"]}
    if all(edge_defs[stream_id].get("derived") for stream_id in effective_edges):
        raise ValueError(
            "scaffold closure requires at least one non-derived edge stream to "
            "construct its pre-closure graph"
        )
    exclusive: dict[str, list[str]] = {}
    for stream_id in effective_edges:
        group = edge_defs[stream_id].get("exclusive_group")
        if group:
            exclusive.setdefault(group, []).append(stream_id)
    conflicts = {group: ids for group, ids in exclusive.items() if len(ids) > 1}
    if conflicts:
        text = "; ".join(f"{group}: {', '.join(ids)}" for group, ids in conflicts.items())
        raise ValueError(f"mutually exclusive edge streams are enabled: {text}")
    return config


def stable_expit(values: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(values, -700.0, 700.0)))


def probability_distribution_summary(
    values: np.ndarray | pd.Series,
    *,
    prior_probability: float,
    output_cutoff: float,
    bin_count: int = 50,
) -> dict[str, Any]:
    """Summarize posterior probabilities for compact, browser-safe plotting.

    Sending a full edge vector to the GUI would make multi-million-pair runs
    unnecessarily large. Fixed bins spanning the complete probability domain
    preserve a directly comparable distribution while keeping the preview
    payload small and deterministic.
    """
    probabilities = np.asarray(values, dtype=float).reshape(-1)
    probabilities = probabilities[np.isfinite(probabilities)]
    if bin_count < 1:
        raise ValueError("probability histogram bin count must be positive")
    if np.any((probabilities < -1e-12) | (probabilities > 1.0 + 1e-12)):
        raise ValueError("probability distribution contains values outside [0, 1]")
    probabilities = np.clip(probabilities, 0.0, 1.0)
    edges = np.linspace(0.0, 1.0, bin_count + 1)
    counts, edges = np.histogram(probabilities, bins=edges)
    return {
        "hypothesis_count": int(len(probabilities)),
        "minimum": float(probabilities.min()) if len(probabilities) else None,
        "mean": float(probabilities.mean()) if len(probabilities) else None,
        "maximum": float(probabilities.max()) if len(probabilities) else None,
        "prior_probability": float(prior_probability),
        "output_probability_cutoff_exclusive": float(output_cutoff),
        "at_exact_prior_count": int(
            np.isclose(
                probabilities,
                prior_probability,
                rtol=0.0,
                atol=1e-12,
            ).sum()
        ),
        "below_prior_count": int((probabilities < prior_probability).sum()),
        "above_output_cutoff_count": int((probabilities > output_cutoff).sum()),
        "bin_edges": [float(value) for value in edges],
        "bin_counts": [int(value) for value in counts],
        "binning": f"{bin_count} equal-width bins over [0, 1]",
    }


def evidence_factor_distribution_summary(
    values: np.ndarray | pd.Series,
    *,
    distribution_scope: str,
    bin_count: int = 41,
) -> dict[str, Any]:
    """Summarize applied Bayes factors on a neutral-centered log2 scale.

    Evidence streams use incompatible raw units, whereas every integrated
    stream ultimately contributes a positive Bayes factor. Summarizing the
    applied factors therefore gives the GUI a common, auditable distribution:
    BF < 1 refutes, BF = 1 is neutral, and BF > 1 supports. Quantiles allow an
    inspected hypothesis to be placed approximately without serializing a
    full edge-sized vector into the browser payload.
    """
    factors = np.asarray(values, dtype=float).reshape(-1)
    factors = factors[np.isfinite(factors)]
    if bin_count < 1:
        raise ValueError("evidence-factor histogram bin count must be positive")
    if np.any(factors <= 0):
        raise ValueError("evidence-factor distribution contains non-positive values")
    if not len(factors):
        return {
            "scale": "log2_bayes_factor",
            "distribution_scope": distribution_scope,
            "hypothesis_count": 0,
            "minimum": None,
            "median": None,
            "mean": None,
            "maximum": None,
            "refuting_count": 0,
            "neutral_count": 0,
            "supporting_count": 0,
            "bin_edges_log2": [],
            "bin_counts": [],
            "quantile_percentages": [],
            "quantile_bayes_factors": [],
            "value_counts": [],
        }

    log2_factors = np.log2(factors)
    span = max(float(np.max(np.abs(log2_factors))), 0.25)
    edges = np.linspace(-span, span, bin_count + 1)
    counts, edges = np.histogram(log2_factors, bins=edges)
    quantile_percentages = np.linspace(0.0, 100.0, 101)
    quantile_values = np.quantile(factors, quantile_percentages / 100.0)
    neutral = np.isclose(factors, 1.0, rtol=0.0, atol=1e-12)
    sample_size = min(len(factors), 8192)
    sample_indices = np.linspace(0, len(factors) - 1, sample_size, dtype=int)
    sampled_unique = np.unique(factors[sample_indices])
    if len(sampled_unique) <= 64:
        unique_values, unique_counts = np.unique(factors, return_counts=True)
    else:
        unique_values = np.asarray([], dtype=float)
        unique_counts = np.asarray([], dtype=int)
    compact_value_counts = [
        {"bayes_factor": float(value), "count": int(count)}
        for value, count in zip(unique_values, unique_counts, strict=True)
    ] if len(unique_values) <= 64 else []
    return {
        "scale": "log2_bayes_factor",
        "distribution_scope": distribution_scope,
        "hypothesis_count": int(len(factors)),
        "minimum": float(factors.min()),
        "median": float(np.median(factors)),
        "mean": float(factors.mean()),
        "maximum": float(factors.max()),
        "refuting_count": int((factors < 1.0 - 1e-12).sum()),
        "neutral_count": int(neutral.sum()),
        "supporting_count": int((factors > 1.0 + 1e-12).sum()),
        "bin_edges_log2": [float(value) for value in edges],
        "bin_counts": [int(value) for value in counts],
        "quantile_percentages": [float(value) for value in quantile_percentages],
        "quantile_bayes_factors": [float(value) for value in quantile_values],
        "value_counts": compact_value_counts,
        "binning": f"{bin_count} equal-width bins on a symmetric log2(BF) domain",
    }


def complement_minimum_likelihood(
    values: np.ndarray | pd.Series,
    thresholds: np.ndarray | pd.Series | float,
) -> np.ndarray:
    """Apply the project's floor-preserving Gaussian complement kernel."""
    x = np.asarray(values, dtype=float)
    tq = np.asarray(thresholds, dtype=float)
    result = np.full(np.broadcast_shapes(x.shape, tq.shape), NEUTRAL_LIKELIHOOD)
    x_broadcast, tq_broadcast = np.broadcast_arrays(x, tq)
    valid = np.isfinite(x_broadcast) & np.isfinite(tq_broadcast) & (tq_broadcast > 0)
    if np.any(valid):
        z = np.maximum(x_broadcast[valid], 0.0) / tq_broadcast[valid]
        result[valid] = np.maximum(
            NEUTRAL_LIKELIHOOD,
            1.0 - np.exp(-0.5 * np.square(z)),
        )
    return result


def continuous_complement_bayes_factor(
    values: np.ndarray | pd.Series,
    thresholds: np.ndarray | pd.Series | float,
    *,
    minimum_bayes_factor: float,
) -> np.ndarray:
    """Score quantitative evidence without the historical 0.5 likelihood floor.

    The raw likelihood is ``1 - exp(-0.5 * (x / Tq)^2)`` and its neutral
    reference is 0.5, so the returned Bayes factor is raw likelihood / 0.5.
    A true zero would make log-odds updates singular; it is therefore replaced
    by the configured small positive BF floor. Missing values must be converted
    to zero by the caller only inside the source's explicit eligibility scope.
    """
    if not math.isfinite(minimum_bayes_factor) or not 0 < minimum_bayes_factor <= 1:
        raise ValueError("minimum continuous Bayes factor must be in (0, 1]")
    x = np.asarray(values, dtype=float)
    tq = np.asarray(thresholds, dtype=float)
    x_broadcast, tq_broadcast = np.broadcast_arrays(x, tq)
    factors = np.full(x_broadcast.shape, minimum_bayes_factor, dtype=float)
    valid = np.isfinite(x_broadcast) & np.isfinite(tq_broadcast) & (tq_broadcast > 0)
    if np.any(valid):
        z = np.maximum(x_broadcast[valid], 0.0) / tq_broadcast[valid]
        likelihood = 1.0 - np.exp(-0.5 * np.square(z))
        factors[valid] = np.maximum(
            minimum_bayes_factor,
            likelihood / NEUTRAL_LIKELIHOOD,
        )
    return factors


def _stream_multiplier(state: dict[str, Any]) -> float:
    return float(state.get("tq_multiplier", 1.0))


def _aligned_series(
    table: pd.DataFrame,
    key_column: str,
    value_column: str,
    keys: pd.Series,
) -> np.ndarray:
    lookup = table.drop_duplicates(key_column).set_index(key_column)[value_column]
    return pd.to_numeric(keys.map(lookup), errors="coerce").to_numpy(float)


def _node_stream_raw_values_and_thresholds(
    project: Path,
    factors: pd.DataFrame,
    definition: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Return the quantitative x and positive Tq arrays for one node stream."""
    genes = factors["gene_symbol"].astype(str)
    handler = definition["normalization"]["handler"]
    if handler == "protein_abundance":
        raw = pd.read_csv(project / PROTEIN_RAW_RELATIVE, sep="\t")
        values = _aligned_series(
            raw, "gene_symbol", "linear_relative_abundance", genes
        )
        threshold = float(pd.to_numeric(raw["T_q"], errors="coerce").dropna().iloc[0])
        return values, np.full(len(factors), threshold, dtype=float)
    if handler == "pc_transcript":
        raw = pd.read_csv(project / PC_RAW_RELATIVE, sep="\t")
        values = _aligned_series(raw, "gene_symbol", "pc_median_tpm", genes)
        threshold = float(pd.to_numeric(raw["T_q"], errors="coerce").dropna().iloc[0])
        return values, np.full(len(factors), threshold, dtype=float)
    if handler == "kinase_activity":
        values = pd.to_numeric(
            factors["kinase_absolute_lfc"], errors="coerce"
        ).to_numpy(float)
        return values, np.full(len(factors), 0.17, dtype=float)
    if handler == "phosphoprotein_response":
        raw = pd.read_csv(project / PHOSPHOPROTEIN_RAW_RELATIVE, sep="\t")
        values = _aligned_series(raw, "gene_symbol", "max_absolute_lfc", genes)
        thresholds = _aligned_series(raw, "gene_symbol", "site_T_q", genes)
    elif handler == "site_level_phosphoprotein_response":
        values = pd.to_numeric(
            factors[definition["raw_value_column"]], errors="coerce"
        ).to_numpy(float)
        thresholds = pd.to_numeric(
            factors[definition["tq_column"]], errors="coerce"
        ).to_numpy(float)
    elif handler == "collecting_duct_abundance":
        values = pd.to_numeric(
            factors[definition["raw_value_column"]], errors="coerce"
        ).to_numpy(float)
        thresholds = pd.to_numeric(
            factors[definition["tq_column"]], errors="coerce"
        ).to_numpy(float)
    else:
        raise ValueError(f"unsupported node normalization handler: {handler}")

    positive_thresholds = thresholds[np.isfinite(thresholds) & (thresholds > 0)]
    if not len(positive_thresholds):
        raise ValueError(f"node stream {definition['id']} has no stored positive Tq")
    # Current phosphoproteomic and segment streams use a single source-wide Tq.
    # Filling missing candidate rows with that positive reference lets x=0 be
    # scored without pretending the missing row supplied its own threshold.
    threshold_reference = float(np.median(positive_thresholds))
    thresholds = np.where(
        np.isfinite(thresholds) & (thresholds > 0),
        thresholds,
        threshold_reference,
    )
    return values, thresholds


def node_stream_values(
    project: Path,
    factors: pd.DataFrame,
    definition: dict[str, Any],
    state: dict[str, Any],
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> np.ndarray:
    """Return one node stream after applying its independent Tq setting."""
    multiplier = _stream_multiplier(state)
    if continuous_negative:
        values, thresholds = _node_stream_raw_values_and_thresholds(
            project, factors, definition
        )
        eligible, _ = node_stream_observation_masks(factors, definition)
        source_neutral = float(definition["neutral_value"])
        result = np.full(len(factors), source_neutral, dtype=float)
        calculable = eligible & np.isfinite(thresholds) & (thresholds > 0)
        zero_inclusive_values = np.where(np.isfinite(values), values, 0.0)
        result[calculable] = source_neutral * continuous_complement_bayes_factor(
            zero_inclusive_values[calculable],
            thresholds[calculable] * multiplier,
            minimum_bayes_factor=minimum_bayes_factor,
        )
        return result
    if math.isclose(multiplier, 1.0, rel_tol=0.0, abs_tol=1e-15):
        return pd.to_numeric(
            factors[definition["column"]], errors="coerce"
        ).fillna(definition["neutral_value"]).to_numpy(float)

    handler = definition["normalization"]["handler"]
    genes = factors["gene_symbol"].astype(str)
    if handler == "protein_abundance":
        raw = pd.read_csv(project / PROTEIN_RAW_RELATIVE, sep="\t")
        values = _aligned_series(
            raw, "gene_symbol", "linear_relative_abundance", genes
        )
        base_tq = float(pd.to_numeric(raw["T_q"], errors="coerce").dropna().iloc[0])
        return complement_minimum_likelihood(values, base_tq * multiplier)
    if handler == "pc_transcript":
        raw = pd.read_csv(project / PC_RAW_RELATIVE, sep="\t")
        values = _aligned_series(raw, "gene_symbol", "pc_median_tpm", genes)
        base_tq = float(pd.to_numeric(raw["T_q"], errors="coerce").dropna().iloc[0])
        observed_positive = np.isfinite(values) & (values > 0)
        result = np.full(len(factors), NEUTRAL_LIKELIHOOD)
        result[observed_positive] = complement_minimum_likelihood(
            values[observed_positive], base_tq * multiplier
        )
        return result
    if handler == "kinase_activity":
        values = pd.to_numeric(
            factors["kinase_absolute_lfc"], errors="coerce"
        ).to_numpy(float)
        observed = factors["kinase_evidence_observed"].fillna(False).astype(bool).to_numpy()
        result = np.ones(len(factors), dtype=float)
        result[observed] = complement_minimum_likelihood(
            values[observed], 0.17 * multiplier
        ) / NEUTRAL_LIKELIHOOD
        return result
    if handler == "phosphoprotein_response":
        raw = pd.read_csv(project / PHOSPHOPROTEIN_RAW_RELATIVE, sep="\t")
        values = _aligned_series(raw, "gene_symbol", "max_absolute_lfc", genes)
        thresholds = _aligned_series(raw, "gene_symbol", "site_T_q", genes)
        observed = np.isfinite(values) & np.isfinite(thresholds) & (thresholds > 0)
        result = np.ones(len(factors), dtype=float)
        result[observed] = complement_minimum_likelihood(
            values[observed], thresholds[observed] * multiplier
        ) / NEUTRAL_LIKELIHOOD
        return result
    if handler == "site_level_phosphoprotein_response":
        values = pd.to_numeric(
            factors[definition["raw_value_column"]], errors="coerce"
        ).to_numpy(float)
        thresholds = pd.to_numeric(
            factors[definition["tq_column"]], errors="coerce"
        ).to_numpy(float)
        observed = _boolean_mask(factors[definition["observed_column"]])
        observed &= np.isfinite(values) & np.isfinite(thresholds) & (thresholds > 0)
        result = np.ones(len(factors), dtype=float)
        result[observed] = complement_minimum_likelihood(
            values[observed], thresholds[observed] * multiplier
        ) / NEUTRAL_LIKELIHOOD
        return result
    if handler == "collecting_duct_abundance":
        values = pd.to_numeric(
            factors[definition["raw_value_column"]], errors="coerce"
        ).to_numpy(float)
        observed_series = factors[definition["observed_column"]]
        if pd.api.types.is_bool_dtype(observed_series):
            observed = observed_series.fillna(False).to_numpy(bool).copy()
        else:
            observed = (
                observed_series.fillna("")
                .astype(str)
                .str.strip()
                .str.casefold()
                .isin({"true", "1", "yes"})
                .to_numpy(bool)
            )
        base_tq_values = pd.to_numeric(
            factors[definition["tq_column"]], errors="coerce"
        ).dropna()
        if base_tq_values.empty:
            raise ValueError(
                f"node stream {definition['id']} has no stored positive Tq"
            )
        base_tq = float(base_tq_values.iloc[0])
        observed &= np.isfinite(values) & (values > 0)
        result = np.full(len(factors), NEUTRAL_LIKELIHOOD)
        result[observed] = complement_minimum_likelihood(
            values[observed], base_tq * multiplier
        )
        return result
    raise ValueError(f"unsupported node normalization handler: {handler}")


def _boolean_mask(values: pd.Series) -> np.ndarray:
    """Interpret stored boolean-like observation flags without truthy strings."""
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).to_numpy(bool).copy()
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"true", "1", "yes"})
        .to_numpy(bool)
    )


def node_stream_negative_evidence_policy(
    definition: dict[str, Any],
) -> dict[str, Any] | None:
    """Return an explicit or conventional observation policy for a node stream."""
    policy = definition.get("negative_evidence")
    if not policy and definition.get("observed_column"):
        policy = {
            "observed_column": definition["observed_column"],
            "eligibility": "all_candidates",
            "absence_definition": (
                "No positive segment-specific abundance mapped to the candidate."
            ),
        }
    return policy


def node_stream_observation_masks(
    factors: pd.DataFrame,
    definition: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (eligible, observed) masks for optional negative node evidence.

    Eligibility is explicit: most abundance streams can assess every protein
    candidate, whereas kinase activity can assess only kinase-annotated nodes.
    Nodes outside a stream's scope remain neutral.
    """
    policy = node_stream_negative_evidence_policy(definition)
    if not policy:
        return np.zeros(len(factors), dtype=bool), np.zeros(len(factors), dtype=bool)

    observed_column = str(policy.get("observed_column", "")).strip()
    if observed_column not in factors.columns:
        raise ValueError(
            f"node stream {definition['id']} lacks observation column "
            f"{observed_column!r}"
        )
    observed = _boolean_mask(factors[observed_column])

    eligibility = str(policy.get("eligibility", "all_candidates"))
    if eligibility == "all_candidates":
        eligible = np.ones(len(factors), dtype=bool)
    elif eligibility == "ontology_class":
        ontology_class = str(policy.get("ontology_class", "")).strip()
        if not ontology_class:
            raise ValueError(
                f"node stream {definition['id']} has no negative-evidence ontology class"
            )
        eligible = (
            factors["node_classes"]
            .fillna("")
            .astype(str)
            .str.split(";")
            .apply(lambda classes: ontology_class in classes)
            .to_numpy(bool)
        )
    else:
        raise ValueError(
            f"node stream {definition['id']} has unsupported negative-evidence "
            f"eligibility {eligibility!r}"
        )
    observed &= eligible
    return eligible, observed


def load_node_factor_catalog(
    project: Path,
    registry: dict[str, Any],
) -> pd.DataFrame:
    """Load the core node catalog and merge each distinct auxiliary table once."""
    factors = pd.read_csv(project / NODE_FACTORS_RELATIVE, sep="\t")
    loaded_paths: set[Path] = set()
    for definition in registry["node_streams"]:
        relative = definition.get("factor_file")
        if not relative:
            continue
        path = (project / relative).resolve()
        if path in loaded_paths:
            continue
        auxiliary = pd.read_csv(path, sep="\t")
        if "gene_symbol" not in auxiliary.columns:
            raise ValueError(f"node factor table lacks gene_symbol: {path}")
        if auxiliary["gene_symbol"].duplicated().any():
            raise ValueError(f"node factor table has duplicate gene symbols: {path}")
        collisions = set(factors.columns).intersection(auxiliary.columns) - {
            "gene_symbol"
        }
        if collisions:
            raise ValueError(
                f"node factor table has conflicting columns {sorted(collisions)}: {path}"
            )
        factors = factors.merge(
            auxiliary,
            on="gene_symbol",
            how="left",
            validate="one_to_one",
        )
        loaded_paths.add(path)
    return factors


def select_nodes(
    project: Path,
    registry: dict[str, Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    factors = load_node_factor_catalog(project, registry)
    universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    stream_defs = {stream["id"]: stream for stream in registry["node_streams"]}
    tolerance = config["node_integration"]["non_neutral_tolerance"]
    prior_probability = float(config["node_integration"]["prior_probability"])
    output_cutoff = float(config["node_integration"]["output_probability_cutoff"])
    penalize_unobserved = bool(
        config["node_integration"]["penalize_unobserved"]
    )
    unobserved_bayes_factor = float(
        config["node_integration"]["unobserved_bayes_factor"]
    )
    global_continuous_negative = bool(
        config["node_integration"]["continuous_negative_evidence"]
    )
    continuous_bf_floor = float(
        config["node_integration"]["continuous_bayes_factor_floor"]
    )
    prior_log_odds = math.log(prior_probability / (1.0 - prior_probability))
    log_odds = np.full(len(factors), prior_log_odds, dtype=float)
    if "initial_prior" in factors.columns:
        factors["catalog_historical_initial_prior"] = factors["initial_prior"]
    factors["initial_prior"] = prior_probability
    any_non_neutral = np.zeros(len(factors), dtype=bool)
    any_positive_support = np.zeros(len(factors), dtype=bool)
    active_streams: list[dict[str, Any]] = []
    stream_audit_columns: dict[str, np.ndarray] = {}
    for stream_id, state in config["node_streams"].items():
        if not state["enabled"] or state["weight"] <= 0:
            continue
        definition = stream_defs[stream_id]
        stream_continuous_negative = (
            global_continuous_negative
            or bool(state.get("continuous_negative_evidence", False))
        )
        values = node_stream_values(
            project,
            factors,
            definition,
            state,
            continuous_negative=stream_continuous_negative,
            minimum_bayes_factor=continuous_bf_floor,
        )
        if (values <= 0).any():
            raise ValueError(f"node stream {stream_id} contains a non-positive factor")
        neutral_value = float(definition["neutral_value"])
        if not math.isfinite(neutral_value) or neutral_value <= 0.0:
            raise ValueError(f"node stream {stream_id} has an invalid neutral value")
        source_bayes_factors = values / neutral_value
        bayes_factors = source_bayes_factors.copy()
        eligible, observed = node_stream_observation_masks(factors, definition)
        unobserved_eligible = eligible & ~observed
        penalty_applied = (
            (penalize_unobserved and not stream_continuous_negative)
            & unobserved_eligible
        )
        bayes_factors[penalty_applied] = unobserved_bayes_factor
        continuous_negative_applied = (
            stream_continuous_negative
            & eligible
            & (bayes_factors < 1.0 - tolerance)
        )
        weight = float(state["weight"])
        weighted_log_bf = weight * np.log(bayes_factors)
        log_odds += weighted_log_bf
        non_neutral = np.abs(bayes_factors - 1.0) > tolerance
        positive_support = bayes_factors > 1.0 + tolerance
        any_non_neutral |= non_neutral
        any_positive_support |= positive_support
        stream_audit_columns.update(
            {
                f"gui_{stream_id}_non_neutral": non_neutral,
                f"gui_{stream_id}_factor": values,
                f"gui_{stream_id}_source_bayes_factor": source_bayes_factors,
                f"gui_{stream_id}_bayes_factor": bayes_factors,
                f"gui_{stream_id}_weighted_log_bayes_factor": weighted_log_bf,
                f"gui_{stream_id}_negative_evidence_eligible": eligible,
                f"gui_{stream_id}_observed": observed,
                f"gui_{stream_id}_unobserved_penalty_applied": penalty_applied,
                f"gui_{stream_id}_continuous_negative_evidence_applied": (
                    continuous_negative_applied
                ),
            }
        )
        active_streams.append(
            {
                "id": stream_id,
                "label": definition["label"],
                "weight": weight,
                "tq_multiplier": _stream_multiplier(state),
                "normalization_reference": definition["normalization"]["reference"],
                "dependence_group": definition.get("dependence_group"),
                "non_neutral_candidates": int(non_neutral.sum()),
                "positive_support_candidates": int(positive_support.sum()),
                "negative_evidence_eligible_candidates": int(eligible.sum()),
                "observed_eligible_candidates": int(observed.sum()),
                "unobserved_eligible_candidates": int(unobserved_eligible.sum()),
                "unobserved_penalties_applied": int(penalty_applied.sum()),
                "continuous_negative_candidates": int(
                    continuous_negative_applied.sum()
                ),
                "continuous_negative_evidence": stream_continuous_negative,
                "unobserved_bayes_factor": (
                    unobserved_bayes_factor if penalize_unobserved else None
                ),
                "negative_evidence_scope": node_stream_negative_evidence_policy(
                    definition
                ),
                "source_neutral_value": neutral_value,
                "factor_distribution": evidence_factor_distribution_summary(
                    bayes_factors,
                    distribution_scope="all modeled protein candidates",
                ),
            }
        )
    if stream_audit_columns:
        factors = pd.concat(
            [factors, pd.DataFrame(stream_audit_columns, index=factors.index)],
            axis=1,
        )
    posterior = stable_expit(log_odds)
    selected_probability = posterior > output_cutoff
    factors["gui_initial_prior_probability"] = prior_probability
    factors["gui_posterior"] = posterior
    factors["gui_rank"] = (
        pd.Series(posterior).rank(method="min", ascending=False).astype(int).to_numpy()
    )
    factors["gui_any_selected_stream_non_neutral"] = any_non_neutral
    factors["gui_any_selected_stream_positive_support"] = any_positive_support
    factors["gui_above_node_probability_cutoff"] = selected_probability

    seed_symbols = set(universe["symbol"])
    selected_protein_symbols = set(
        factors.loc[selected_probability, "gene_symbol"].astype(str)
    )
    selected_seed_symbols = selected_protein_symbols.intersection(seed_symbols)
    added_symbols = selected_protein_symbols.difference(seed_symbols)
    factors["gui_available_in_edge_catalog"] = factors["gene_symbol"].isin(seed_symbols)
    factors["gui_selected_in_graph"] = selected_probability
    include_messengers = config["node_integration"]["include_second_messengers"]
    selected_rows = universe[universe["symbol"].isin(selected_seed_symbols)].copy()
    selected_rows["gui_selection_reason"] = "posterior_above_node_cutoff"
    if added_symbols:
        liberal = pd.read_csv(
            project / "data/node_selection/mouse_signaling_nodes_liberal.tsv",
            sep="\t",
            dtype=str,
        ).fillna("")
        liberal = liberal.drop_duplicates("symbol").set_index("symbol")
        dynamic_source = factors.loc[
            factors["gene_symbol"].isin(added_symbols)
        ].copy()
        dynamic_source = dynamic_source.sort_values(
            ["gui_rank", "gene_symbol"], kind="stable"
        )
        dynamic_rows = pd.DataFrame(
            {
                "symbol": dynamic_source["gene_symbol"].astype(str),
                "display_symbol": dynamic_source["gene_symbol"].astype(str),
                "name": dynamic_source["node_name"].fillna("").astype(str),
                "classes": dynamic_source["node_classes"].fillna("").astype(str),
                "node_type": "protein",
                "stable_id": "",
                "scope_tier": "incrementally_characterized",
                "selection_basis": "binary_node_posterior_above_cutoff",
                "bayesian_score_status": "scored",
                "included_by_curated_rule": False,
                "source_url": "",
                "scope_note": "Added beyond the immutable 891-node seed catalog.",
            }
        )
        dynamic_rows["selected_uniprot"] = dynamic_rows["symbol"].map(
            liberal["uniprot"]
        ).fillna("")
        dynamic_rows["gui_selection_reason"] = "posterior_above_node_cutoff_incremental"
        selected_rows = pd.concat([selected_rows, dynamic_rows], ignore_index=True)
    if include_messengers:
        messengers = universe[universe["node_type"] == "molecule"].copy()
        messengers["gui_selection_reason"] = "curated_second_messenger"
        selected_rows = pd.concat([selected_rows, messengers], ignore_index=True)
    current_node_scores = factors[
        [
            "gene_symbol",
            "gui_initial_prior_probability",
            "gui_posterior",
            "gui_rank",
            "gui_any_selected_stream_non_neutral",
            "gui_any_selected_stream_positive_support",
            "gui_above_node_probability_cutoff",
        ]
    ].rename(columns={"gene_symbol": "symbol"})
    selected_rows = selected_rows.merge(
        current_node_scores,
        on="symbol",
        how="left",
        validate="many_to_one",
    )
    order = {symbol: index for index, symbol in enumerate(universe["symbol"])}
    selected_rows["gui_universe_index"] = selected_rows["symbol"].map(order)
    dynamic_order = {
        symbol: len(order) + index
        for index, symbol in enumerate(
            sorted(added_symbols, key=lambda value: (value.casefold(), value))
        )
    }
    selected_rows["gui_universe_index"] = selected_rows["gui_universe_index"].fillna(
        selected_rows["symbol"].map(dynamic_order)
    )
    selected_rows = selected_rows.sort_values("gui_universe_index").reset_index(drop=True)
    summary = {
        "candidate_protein_count": int(len(factors)),
        "selected_protein_count": int(len(selected_protein_symbols)),
        "selected_seed_protein_count": int(len(selected_seed_symbols)),
        "incrementally_added_protein_count": int(len(added_symbols)),
        "non_neutral_candidates_outside_edge_catalog": int(len(added_symbols)),
        "curated_second_messenger_count": int(
            (selected_rows["node_type"] == "molecule").sum()
        ),
        "selected_node_count": int(len(selected_rows)),
        "active_streams": active_streams,
        "node_prior_probability": prior_probability,
        "node_output_probability_cutoff": output_cutoff,
        "penalize_unobserved": penalize_unobserved,
        "unobserved_bayes_factor": unobserved_bayes_factor,
        "continuous_negative_evidence": global_continuous_negative,
        "continuous_negative_streams": [
            stream["id"]
            for stream in active_streams
            if stream["continuous_negative_evidence"]
        ],
        "continuous_bayes_factor_floor": continuous_bf_floor,
        "total_unobserved_penalties_applied": int(
            sum(stream["unobserved_penalties_applied"] for stream in active_streams)
        ),
        "total_continuous_negative_applications": int(
            sum(stream["continuous_negative_candidates"] for stream in active_streams)
        ),
        "candidates_below_prior": int((posterior < prior_probability).sum()),
        "posterior_minimum": float(posterior.min()),
        "posterior_mean": float(posterior.mean()),
        "posterior_maximum": float(posterior.max()),
        "neutral_posterior_count": int(
            np.isclose(posterior, prior_probability, rtol=0.0, atol=tolerance).sum()
        ),
        "probability_distribution": probability_distribution_summary(
            posterior,
            prior_probability=prior_probability,
            output_cutoff=output_cutoff,
        ),
        "posterior_probabilities_are_independent": True,
        "selection_rule": (
            "Each protein is an independent present-versus-absent hypothesis. Every "
            "protein begins at the configured Bernoulli prior (default 0.5). Source "
            "scores are divided by their neutral values to obtain BF=1 at neutrality; "
            "weighted Bayes factors multiply prior odds. The fixed-absence mode gives "
            "each eligible nondetection one configured BF<1. Alternatively, continuous "
            "negative evidence removes the 0.5 likelihood floor and scores both weak "
            "detections and x=0 nondetections with the same Tq kernel, bounded only by "
            "the numerical BF floor. Out-of-scope nodes remain at BF=1. A protein is selected when "
            "its posterior is strictly above the configured node cutoff (default 0.5). "
            "Posteriors are not normalized across proteins."
        ),
        "catalog_constraint": (
            "The validated 891-node graph is the immutable seed. Newly non-neutral "
            "proteins are inserted into the graph after all new unordered pairs are "
            "characterized and written to the persistent incremental edge cache."
        ),
    }
    return factors, selected_rows, summary


def _calibration_parameter_specs(
    registry: dict[str, Any],
    config: dict[str, Any],
    group: str,
) -> list[dict[str, Any]]:
    """Describe bounded primary-stream weights and scale multipliers."""
    definitions = {item["id"]: item for item in registry[group]}
    specs: list[dict[str, Any]] = []
    for stream_id, state in config[group].items():
        definition = definitions[stream_id]
        if not state["enabled"] or definition.get("derived"):
            continue
        specs.append(
            {
                "key": f"{stream_id}:weight",
                "stream_id": stream_id,
                "parameter": "weight",
                "current": float(state["weight"]),
                "preferred": 1.0,
                "bounds": CALIBRATION_WEIGHT_BOUNDS,
                "scale": 0.5,
                "transform": "linear",
            }
        )
        if definition.get("normalization", {}).get("user_control", True):
            specs.append(
                {
                    "key": f"{stream_id}:tq_multiplier",
                    "stream_id": stream_id,
                    "parameter": "tq_multiplier",
                    "current": _stream_multiplier(state),
                    "preferred": float(state["preferred_tq_multiplier"]),
                    "bounds": _calibration_tq_bounds(definition),
                    "scale": math.log(2.0),
                    "transform": "log",
                }
            )
    return specs


def _apply_calibrated_parameters(
    config: dict[str, Any],
    group: str,
    parameter_rows: list[dict[str, Any]],
) -> None:
    for row in parameter_rows:
        stream_id = str(row["stream_id"])
        parameter = str(row["parameter"])
        config[group][stream_id][parameter] = float(row["fitted"])


def _prepare_node_calibration_context(
    project: Path,
    factors: pd.DataFrame,
    definition: dict[str, Any],
    target_indices: np.ndarray,
) -> dict[str, np.ndarray]:
    values, thresholds = _node_stream_raw_values_and_thresholds(
        project, factors, definition
    )
    eligible, observed = node_stream_observation_masks(factors, definition)
    return {
        "values": values[target_indices],
        "thresholds": thresholds[target_indices],
        "eligible": eligible[target_indices],
        "observed": observed[target_indices],
    }


def calibrate_node_parameters(
    project: Path,
    registry: dict[str, Any],
    config: dict[str, Any],
    *,
    progress: Callable[[str, float], None] | None = None,
    check_cancel: Callable[[], None] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Fit node stream settings against supplied known-present candidates."""
    requested = config["calibration"]["known_nodes"]
    factors = load_node_factor_catalog(project, registry)
    symbol_lookup = {
        symbol.casefold(): (index, symbol)
        for index, symbol in enumerate(factors["gene_symbol"].astype(str))
    }
    missing = [value for value in requested if value.casefold() not in symbol_lookup]
    if missing:
        raise ValueError(
            "known calibration nodes are not modeled signaling candidates: "
            + ", ".join(missing)
        )
    resolved = [symbol_lookup[value.casefold()][1] for value in requested]
    target_indices = np.asarray(
        [symbol_lookup[value.casefold()][0] for value in requested], dtype=int
    )
    definitions = {item["id"]: item for item in registry["node_streams"]}
    contexts = {
        stream_id: _prepare_node_calibration_context(
            project, factors, definitions[stream_id], target_indices
        )
        for stream_id, state in config["node_streams"].items()
        if state["enabled"]
    }
    specs = _calibration_parameter_specs(registry, config, "node_streams")
    prior = float(config["node_integration"]["prior_probability"])
    prior_log_odds = math.log(prior / (1.0 - prior))
    penalize = bool(config["node_integration"]["penalize_unobserved"])
    missing_factor = float(config["node_integration"]["unobserved_bayes_factor"])
    global_continuous_negative = bool(
        config["node_integration"]["continuous_negative_evidence"]
    )
    continuous_floor = float(
        config["node_integration"]["continuous_bayes_factor_floor"]
    )

    def probabilities(settings: dict[str, float]) -> np.ndarray:
        log_odds = np.full(len(target_indices), prior_log_odds, dtype=float)
        for stream_id, context in contexts.items():
            multiplier = settings[f"{stream_id}:tq_multiplier"]
            weight = settings[f"{stream_id}:weight"]
            stream_continuous_negative = (
                global_continuous_negative
                or bool(
                    config["node_streams"][stream_id].get(
                        "continuous_negative_evidence", False
                    )
                )
            )
            bayes_factors = np.ones(len(target_indices), dtype=float)
            observed = context["observed"].astype(bool)
            if stream_continuous_negative:
                calculable = (
                    context["eligible"].astype(bool)
                    & np.isfinite(context["thresholds"])
                    & (context["thresholds"] > 0)
                )
                zero_inclusive = np.where(
                    np.isfinite(context["values"]), context["values"], 0.0
                )
                bayes_factors[calculable] = continuous_complement_bayes_factor(
                    zero_inclusive[calculable],
                    context["thresholds"][calculable] * multiplier,
                    minimum_bayes_factor=continuous_floor,
                )
            else:
                calculable = (
                    observed
                    & np.isfinite(context["values"])
                    & np.isfinite(context["thresholds"])
                    & (context["thresholds"] > 0)
                )
                bayes_factors[calculable] = complement_minimum_likelihood(
                    context["values"][calculable],
                    context["thresholds"][calculable] * multiplier,
                ) / NEUTRAL_LIKELIHOOD
            if penalize and not stream_continuous_negative:
                bayes_factors[context["eligible"].astype(bool) & ~observed] = missing_factor
            log_odds += weight * np.log(bayes_factors)
        return stable_expit(log_odds)

    summary, trace = bounded_powell_positive_calibration(
        specs,
        probabilities,
        target_count=len(resolved),
        regularization_strength=config["calibration"]["regularization_strength"],
        multistart_count=config["calibration"]["multistart_count"],
        progress=progress,
        check_cancel=check_cancel,
        stage_label="node evidence",
    )
    _apply_calibrated_parameters(config, "node_streams", summary["parameters"])
    target_table = pd.DataFrame(
        {
            "requested_node": requested,
            "resolved_symbol": resolved,
            "initial_probability": summary.pop("initial_probabilities"),
            "fitted_probability": summary.pop("final_probabilities"),
        }
    )
    summary["known_nodes"] = resolved
    summary["scope_note"] = (
        "Only known-present node controls entered the fit loss. All other candidates "
        "remained unlabeled and were evaluated only after fitting."
    )
    return summary, target_table, trace


def resolve_target_for_graph(
    target: str,
    universe: pd.DataFrame,
    project: Path,
) -> tuple[str, bool, dict[str, Any] | None]:
    symbols_ci = {symbol.casefold(): symbol for symbol in universe["symbol"].astype(str)}
    if target.casefold() in symbols_ci:
        return symbols_ci[target.casefold()], False, None
    symbol, entry, method = resolve_target(target, project / UNIPROT_RELATIVE)
    return symbol, True, {
        "target_mapping_method": method,
        "target_uniprot": str(entry["Entry"]),
        "target_protein_name": str(entry["Protein names"]),
    }


def append_requested_graph_endpoint(
    query: str,
    graph_symbols: list[str],
    graph_metadata: pd.DataFrame,
    universe: pd.DataFrame,
    node_factors: pd.DataFrame,
    project: Path,
    *,
    reason: str,
) -> tuple[str, pd.DataFrame, dict[str, Any] | None]:
    """Resolve and append an endpoint needed by paths or calibration."""
    selected_ci = {symbol.casefold(): symbol for symbol in graph_symbols}
    if query.casefold() in selected_ci:
        return selected_ci[query.casefold()], graph_metadata, None
    universe_ci = {
        symbol.casefold(): symbol for symbol in universe["symbol"].astype(str)
    }
    resolution: dict[str, Any] | None = None
    if query.casefold() in universe_ci:
        symbol = universe_ci[query.casefold()]
        row = universe.loc[universe["symbol"].astype(str).eq(symbol)].iloc[0].to_dict()
    else:
        symbol, _, resolution = resolve_target_for_graph(query, universe, project)
        candidate = node_factors.loc[
            node_factors["gene_symbol"].astype(str).str.casefold().eq(symbol.casefold())
        ]
        row = {
            "symbol": symbol,
            "display_symbol": symbol,
            "name": (
                str(candidate.iloc[0]["node_name"])
                if len(candidate)
                else (resolution or {}).get("target_protein_name", "")
            ),
            "classes": (
                str(candidate.iloc[0]["node_classes"])
                if len(candidate)
                else "external_target"
            ),
            "node_type": "protein",
            "selected_uniprot": (resolution or {}).get("target_uniprot", ""),
        }
    row["gui_selection_reason"] = reason
    graph_symbols.append(symbol)
    graph_metadata = pd.concat(
        [graph_metadata, pd.DataFrame([row])], ignore_index=True, sort=False
    ).fillna("")
    return symbol, graph_metadata, resolution


def _read_aligned_matrix(
    path: Path,
    symbols: list[str],
) -> np.ndarray:
    source = pd.read_csv(path, sep="\t", index_col=0)
    if source.index.has_duplicates or source.columns.duplicated().any():
        raise ValueError(f"matrix has duplicate labels: {path}")
    # Protein-only source matrices legitimately omit curated small molecules and
    # newly requested nodes. Their missing raw score is zero; source-specific
    # eligibility still decides whether that zero is negative or out of scope.
    frame = source.reindex(index=symbols, columns=symbols).fillna(0.0)
    return frame.to_numpy(float)


def _upper_factor_table(
    symbols: list[str],
    factors: np.ndarray,
    *,
    include: np.ndarray | None = None,
    retain_included_neutral: bool = False,
) -> pd.DataFrame:
    left, right = np.triu_indices(len(symbols), 1)
    values = np.asarray(factors, dtype=float)[left, right]
    keep = np.abs(values - 1.0) > FACTOR_EPSILON
    if include is not None:
        included = np.asarray(include, dtype=bool)[left, right]
        keep = included if retain_included_neutral else keep & included
    symbol_array = np.asarray(symbols)
    return pd.DataFrame(
        {
            "node_a": symbol_array[left[keep]],
            "node_b": symbol_array[right[keep]],
            "bayes_factor": values[keep],
        }
    )


def _localization_edge_factors(
    project: Path,
    symbols: list[str],
    multiplier: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    profiles = pd.read_csv(project / MPKCCD_PROFILES_RELATIVE, sep="\t").set_index("symbol")
    profiles = profiles.reindex(symbols)
    columns = ["1K", "4K", "17K", "200Kp", "200Ks"]
    observed = profiles["localization_observed"].fillna(False).astype(bool).to_numpy().copy()
    raw = profiles[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    thresholds = pd.to_numeric(profiles["localization_Tq"], errors="coerce").to_numpy(float)
    observed &= np.isfinite(raw).all(axis=1) & np.isfinite(thresholds) & (thresholds > 0)
    likelihood = np.full((len(symbols), len(symbols)), NEUTRAL_LIKELIHOOD)
    indices = np.flatnonzero(observed)
    if len(indices):
        dot = raw[indices] @ raw[indices].T
        tq = thresholds[indices] * multiplier
        if continuous_negative:
            directed_bf = continuous_complement_bayes_factor(
                dot,
                tq[np.newaxis, :],
                minimum_bayes_factor=minimum_bayes_factor,
            )
            symmetric = (directed_bf + directed_bf.T) / 2.0
        else:
            directed = complement_minimum_likelihood(dot, tq[np.newaxis, :])
            symmetric = (directed + directed.T) / 2.0
        likelihood[np.ix_(indices, indices)] = symmetric
    np.fill_diagonal(likelihood, NEUTRAL_LIKELIHOOD)
    if continuous_negative:
        factors = np.ones_like(likelihood)
        factors[np.ix_(indices, indices)] = likelihood[np.ix_(indices, indices)]
        eligible = np.zeros_like(factors, dtype=bool)
        eligible[np.ix_(indices, indices)] = True
        np.fill_diagonal(eligible, False)
        return _upper_factor_table(
            symbols,
            factors,
            include=eligible,
            retain_included_neutral=True,
        )
    return _upper_factor_table(symbols, likelihood / NEUTRAL_LIKELIHOOD)


def _kinase_predictor_edge_factors(
    project: Path,
    symbols: list[str],
    multiplier: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    predictions = pd.read_csv(project / KINASE_PREDICTIONS_RELATIVE, sep="\t")
    used = predictions["used_for_undirected_edge"]
    if used.dtype != bool:
        used = used.astype(str).str.casefold().eq("true")
    predictions = predictions.loc[used].copy()
    symbol_set = set(symbols)
    predictions = predictions.loc[
        predictions["undirected_node_a"].isin(symbol_set)
        & predictions["undirected_node_b"].isin(symbol_set)
    ]
    raw_score = pd.to_numeric(predictions["raw_score"], errors="raise").to_numpy(float)
    tq = pd.to_numeric(predictions["site_Tq_q75"], errors="raise").to_numpy(float)
    if continuous_negative:
        site_bf = continuous_complement_bayes_factor(
            np.maximum(raw_score, 0.0),
            tq * multiplier,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    else:
        site_bf = complement_minimum_likelihood(
            np.maximum(raw_score, 0.0), tq * multiplier
        ) / NEUTRAL_LIKELIHOOD
    predictions["gui_log_bayes_factor"] = np.log(site_bf)
    grouped = predictions.groupby(
        ["undirected_node_a", "undirected_node_b"], sort=False
    )["gui_log_bayes_factor"].sum().reset_index()
    grouped["bayes_factor"] = np.exp(
        np.clip(grouped["gui_log_bayes_factor"].to_numpy(float), -700.0, 700.0)
    )
    if continuous_negative:
        grouped["bayes_factor"] = np.maximum(
            grouped["bayes_factor"].to_numpy(float), minimum_bayes_factor
        )
    grouped = grouped.loc[
        np.abs(grouped["bayes_factor"] - 1.0) > FACTOR_EPSILON
    ].rename(
        columns={"undirected_node_a": "node_a", "undirected_node_b": "node_b"}
    )
    return grouped[["node_a", "node_b", "bayes_factor"]]


def _string_edge_factors(
    project: Path,
    symbols: list[str],
    multiplier: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    score = _read_aligned_matrix(project / STRING_SCORE_MATRIX_RELATIVE, symbols)
    reference = STRING_REFERENCE_SCORE * multiplier
    if not 0 < reference < 1:
        raise ValueError("STRING reference score after scaling must be between 0 and 1")
    observed = score > 0
    safe_score = np.clip(score, 1e-12, 1 - 1e-12)
    score_odds = safe_score / (1.0 - safe_score)
    reference_odds = reference / (1.0 - reference)
    factors = np.ones_like(score)
    factors[observed] = score_odds[observed] / reference_odds
    if continuous_negative:
        factors[observed] = np.maximum(factors[observed], minimum_bayes_factor)
    return _upper_factor_table(
        symbols,
        factors,
        include=observed,
        retain_included_neutral=continuous_negative,
    )


def _hpa_binary_profiles(
    profiles: pd.DataFrame,
    high_confidence: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    location_columns = [column for column in profiles.columns if column.startswith("location__")]
    if high_confidence:
        binary = np.zeros((len(profiles), len(location_columns)), dtype=float)
        column_index = {column.removeprefix("location__"): i for i, column in enumerate(location_columns)}
        for row_index, value in enumerate(profiles["hpa_high_confidence_locations"].fillna("")):
            for location in str(value).split(";"):
                safe = re.sub(r"[^A-Za-z0-9]+", "_", location).strip("_").lower()
                if safe in column_index:
                    binary[row_index, column_index[safe]] = 1.0
        observed = profiles["hpa_high_confidence_profile_observed"].fillna(False).astype(bool).to_numpy().copy()
        threshold_column = "hpa_high_confidence_Tq_q75"
    else:
        binary = profiles[location_columns].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
        observed = profiles["hpa_profile_observed"].fillna(False).astype(bool).to_numpy().copy()
        threshold_column = "hpa_Tq_q75"
    thresholds = pd.to_numeric(profiles[threshold_column], errors="coerce").to_numpy(float)
    observed &= binary.sum(axis=1) > 0
    return binary, observed, thresholds


def _hpa_edge_factors(
    project: Path,
    symbols: list[str],
    multiplier: float,
    *,
    high_confidence: bool,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    profiles = pd.read_csv(project / HPA_PROFILES_RELATIVE, sep="\t").set_index("symbol").reindex(symbols)
    binary, observed, thresholds = _hpa_binary_profiles(profiles, high_confidence)
    likelihood = np.full((len(symbols), len(symbols)), NEUTRAL_LIKELIHOOD)
    indices = np.flatnonzero(observed)
    if len(indices):
        selected = binary[indices]
        unit = selected / np.linalg.norm(selected, axis=1, keepdims=True)
        similarity = unit @ unit.T
        tq = thresholds[indices] * multiplier
        directed = np.full_like(
            similarity,
            1.0 if continuous_negative else NEUTRAL_LIKELIHOOD,
        )
        valid_tq = np.isfinite(tq) & (tq > 0)
        if np.any(valid_tq):
            if continuous_negative:
                directed[valid_tq] = continuous_complement_bayes_factor(
                    similarity[valid_tq],
                    tq[valid_tq, np.newaxis],
                    minimum_bayes_factor=minimum_bayes_factor,
                )
            else:
                directed[valid_tq] = complement_minimum_likelihood(
                    similarity[valid_tq], tq[valid_tq, np.newaxis]
                )
        likelihood[np.ix_(indices, indices)] = (directed + directed.T) / 2.0
    np.fill_diagonal(likelihood, NEUTRAL_LIKELIHOOD)
    if continuous_negative:
        factors = np.ones_like(likelihood)
        factors[np.ix_(indices, indices)] = likelihood[np.ix_(indices, indices)]
        eligible = np.zeros_like(factors, dtype=bool)
        eligible[np.ix_(indices, indices)] = True
        np.fill_diagonal(eligible, False)
        return _upper_factor_table(
            symbols,
            factors,
            include=eligible,
            retain_included_neutral=True,
        )
    return _upper_factor_table(symbols, likelihood / NEUTRAL_LIKELIHOOD)


def _omnipath_edge_factors(
    project: Path,
    symbols: list[str],
    multiplier: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    effort = _read_aligned_matrix(project / OMNIPATH_EFFORT_MATRIX_RELATIVE, symbols)
    observed = effort > 0
    factors = np.ones_like(effort)
    if continuous_negative:
        factors[observed] = continuous_complement_bayes_factor(
            effort[observed],
            6.0 * multiplier,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    else:
        support = 1.0 - np.exp(-0.5 * np.square(effort / (6.0 * multiplier)))
        likelihood = NEUTRAL_LIKELIHOOD + (1.0 - NEUTRAL_LIKELIHOOD) * support
        factors[observed] = likelihood[observed] / NEUTRAL_LIKELIHOOD
    return _upper_factor_table(
        symbols,
        factors,
        include=observed,
        retain_included_neutral=continuous_negative,
    )


def _stitch_bayes_factors(
    scores: np.ndarray,
    reference: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> np.ndarray:
    """Convert STITCH confidence scores to positive-only evidence factors."""
    if not 0 < reference < 1:
        raise ValueError("STITCH reference score after scaling must be between 0 and 1")
    scores = np.asarray(scores, dtype=float)
    if (~np.isfinite(scores) | (scores <= 0) | (scores >= 1)).any():
        raise ValueError("STITCH scores must be finite probabilities between 0 and 1")
    score_odds = scores / (1.0 - scores)
    reference_odds = reference / (1.0 - reference)
    factors = score_odds / reference_odds
    if continuous_negative:
        return np.maximum(minimum_bayes_factor, factors)
    return np.maximum(1.0, factors)


def _stitch_edge_factors(
    project: Path,
    symbols: list[str],
    multiplier: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    evidence = pd.read_csv(project / STITCH_UNIVERSE_EDGES_RELATIVE, sep="\t")
    symbol_set = set(symbols)
    evidence = evidence.loc[
        evidence["node_a"].isin(symbol_set) & evidence["node_b"].isin(symbol_set)
    ].copy()
    reference = STITCH_REFERENCE_SCORE * multiplier
    evidence["bayes_factor"] = _stitch_bayes_factors(
        pd.to_numeric(evidence["stitch_score"], errors="raise").to_numpy(float),
        reference,
        continuous_negative=continuous_negative,
        minimum_bayes_factor=minimum_bayes_factor,
    )
    if not continuous_negative:
        evidence = evidence.loc[
            np.abs(evidence["bayes_factor"] - 1.0) > FACTOR_EPSILON
        ]
    return evidence[["node_a", "node_b", "bayes_factor"]]


def edge_stream_factor_table(
    project: Path,
    definition: dict[str, Any],
    state: dict[str, Any],
    symbols: list[str],
    *,
    include_negative: bool = False,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    """Load default BFs or recompute one stream with its requested scale."""
    if definition.get("derived"):
        raise ValueError(
            f"derived stream {definition['id']} must be calculated from the "
            "pre-closure graph, not loaded as an independent factor table"
        )
    multiplier = _stream_multiplier(state)
    handler = definition["normalization"]["handler"]
    use_raw_negative_string = (
        (include_negative or continuous_negative) and handler == "string_v12"
    )
    if (
        math.isclose(multiplier, 1.0, rel_tol=0.0, abs_tol=1e-15)
        and not use_raw_negative_string
        and not continuous_negative
    ):
        table = pd.read_csv(
            project / definition["factor_file"], sep="\t", compression="gzip"
        )
        return table[["node_a", "node_b", definition["factor_column"]]].rename(
            columns={definition["factor_column"]: "bayes_factor"}
        )
    if handler == "mpkccd_localization":
        return _localization_edge_factors(
            project,
            symbols,
            multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "kinase_predictor":
        return _kinase_predictor_edge_factors(
            project,
            symbols,
            multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "string_v12":
        return _string_edge_factors(
            project,
            symbols,
            multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "hpa_primary":
        return _hpa_edge_factors(
            project,
            symbols,
            multiplier,
            high_confidence=False,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "hpa_high_confidence":
        return _hpa_edge_factors(
            project,
            symbols,
            multiplier,
            high_confidence=True,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "omnipath_core":
        return _omnipath_edge_factors(
            project,
            symbols,
            multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "stitch_secondary_messenger":
        return _stitch_edge_factors(
            project,
            symbols,
            multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    raise ValueError(f"unsupported edge normalization handler: {handler}")


def _metadata_for_graph(
    project: Path,
    graph_symbols: list[str],
    graph_metadata: pd.DataFrame | None,
) -> pd.DataFrame:
    """Return one metadata row per graph symbol without inventing classifications."""
    base = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    frames = [base]
    if graph_metadata is not None:
        frames.append(graph_metadata.copy().fillna(""))
    metadata = pd.concat(frames, ignore_index=True, sort=False)
    if "symbol" not in metadata.columns:
        raise ValueError("graph metadata must contain a symbol column")
    metadata["symbol"] = metadata["symbol"].astype(str)
    metadata = metadata.drop_duplicates("symbol", keep="last").set_index("symbol")
    return metadata.reindex(graph_symbols).fillna("")


def _observed_series(values: pd.Series) -> pd.Series:
    """Interpret booleans without treating the string ``False`` as true."""
    if values.dtype == bool:
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.casefold().isin({"true", "1", "yes"})


def edge_stream_eligibility_matrix(
    project: Path,
    definition: dict[str, Any],
    graph_symbols: list[str],
    *,
    graph_metadata: pd.DataFrame | None = None,
) -> np.ndarray:
    """Identify pairs a source could assess for optional negative evidence.

    Eligibility is deliberately source-specific.  A missing relationship can
    lower the odds only when the source contains the necessary measurements or
    identifiers for both endpoints.  This prevents missing coverage from being
    mistaken for evidence against an edge.
    """
    policy = definition.get("negative_evidence", {})
    eligibility = policy.get("eligibility")
    if not eligibility:
        return np.zeros((len(graph_symbols), len(graph_symbols)), dtype=bool)

    symbols = pd.Index(graph_symbols, dtype=str)
    metadata = _metadata_for_graph(project, graph_symbols, graph_metadata)
    protein = (
        metadata.get("node_type", pd.Series("", index=metadata.index))
        .astype(str)
        .str.casefold()
        .eq("protein")
        .to_numpy(bool)
    )
    classes = metadata.get("classes", pd.Series("", index=metadata.index)).astype(str)
    kinase = classes.map(
        lambda value: "kinase"
        in {token.strip() for token in value.split(";") if token.strip()}
    ).to_numpy(bool)

    node_flag = np.zeros(len(symbols), dtype=bool)
    other_flag: np.ndarray | None = None
    if eligibility == "both_mpkccd_profiles":
        profiles = pd.read_csv(project / MPKCCD_PROFILES_RELATIVE, sep="\t").set_index("symbol")
        aligned = profiles.reindex(symbols)
        observed = _observed_series(aligned["localization_observed"])
        tq = pd.to_numeric(aligned["localization_Tq"], errors="coerce")
        node_flag = (observed & tq.gt(0) & tq.notna()).to_numpy(bool)
    elif eligibility in {"both_hpa_primary_profiles", "both_hpa_high_profiles"}:
        profiles = pd.read_csv(project / HPA_PROFILES_RELATIVE, sep="\t").set_index("symbol")
        aligned = profiles.reindex(symbols)
        if eligibility == "both_hpa_high_profiles":
            observed_column = "hpa_high_confidence_profile_observed"
            tq_column = "hpa_high_confidence_Tq_q75"
        else:
            observed_column = "hpa_profile_observed"
            tq_column = "hpa_Tq_q75"
        observed = _observed_series(aligned[observed_column])
        tq = pd.to_numeric(aligned[tq_column], errors="coerce")
        node_flag = (observed & tq.gt(0) & tq.notna()).to_numpy(bool)
    elif eligibility == "kinase_to_scorable_phosphoprotein":
        sites = pd.read_csv(project / OBSERVED_PHOSPHOSITES_RELATIVE, sep="\t")
        scorable = _observed_series(sites["kinasepredictor_scorable"])
        target_symbols = set(sites.loc[scorable, "symbol"].astype(str))
        node_flag = symbols.isin(target_symbols).astype(bool)
        other_flag = kinase
    elif eligibility == "both_string_mapped_proteins":
        mapping = pd.read_csv(project / STRING_MAPPING_RELATIVE, sep="\t", dtype=str).fillna("")
        mapped_symbols = set(
            mapping.loc[
                mapping["mapping_status"].eq("mapped")
                & mapping["string_id"].str.strip().ne(""),
                "symbol",
            ].astype(str)
        )
        node_flag = protein & symbols.isin(mapped_symbols).astype(bool)
    elif eligibility == "all_protein_pairs":
        node_flag = protein
    elif eligibility == "mapped_messenger_protein_pairs":
        mapping = pd.read_csv(project / STRING_MAPPING_RELATIVE, sep="\t", dtype=str).fillna("")
        mapped_symbols = set(
            mapping.loc[
                mapping["mapping_status"].eq("mapped")
                & mapping["string_id"].str.strip().ne(""),
                "symbol",
            ].astype(str)
        )
        stitch = pd.read_csv(
            project / STITCH_ALL_MOUSE_EDGES_RELATIVE,
            sep="\t",
            compression="gzip",
            usecols=["messenger_node"],
        )
        node_flag = symbols.isin(set(stitch["messenger_node"].astype(str))).astype(bool)
        other_flag = protein & symbols.isin(mapped_symbols).astype(bool)
    else:
        raise ValueError(f"unsupported edge negative-evidence eligibility: {eligibility}")

    node_flag = np.asarray(node_flag, dtype=bool).copy()
    if other_flag is not None:
        other_flag = np.asarray(other_flag, dtype=bool).copy()

    # Overlay scope information for dynamically added nodes.  KinasePredictor
    # targets remain conservative because the current cache does not persist a
    # no-hit scorable-site flag; such nodes are left neutral rather than guessed.
    seed_symbols = set(
        pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", usecols=["symbol"])[
            "symbol"
        ].astype(str)
    )
    if any(symbol not in seed_symbols for symbol in graph_symbols) and eligibility in {
        "both_mpkccd_profiles",
        "both_hpa_primary_profiles",
        "both_hpa_high_profiles",
        "both_string_mapped_proteins",
        "mapped_messenger_protein_pairs",
    }:
        cached = incremental_node_scope_table(project, graph_symbols)
        if not cached.empty:
            cached = cached.set_index("symbol").reindex(symbols)
            if eligibility == "both_mpkccd_profiles":
                dynamic_flag = cached["localization_observed"].fillna(False).astype(bool).to_numpy()
                node_flag |= dynamic_flag
            elif eligibility == "both_hpa_primary_profiles":
                dynamic_flag = cached["hpa_primary_observed"].fillna(False).astype(bool).to_numpy()
                node_flag |= dynamic_flag
            elif eligibility == "both_hpa_high_profiles":
                dynamic_flag = cached["hpa_high_observed"].fillna(False).astype(bool).to_numpy()
                node_flag |= dynamic_flag
            else:
                dynamic_mapped = cached["string_mapped"].fillna(False).astype(bool).to_numpy()
                if eligibility == "both_string_mapped_proteins":
                    node_flag |= protein & dynamic_mapped
                else:
                    assert other_flag is not None
                    other_flag |= protein & dynamic_mapped

    if other_flag is None:
        eligible = np.logical_and.outer(node_flag, node_flag)
    else:
        eligible = np.logical_or(
            np.logical_and.outer(node_flag, other_flag),
            np.logical_and.outer(other_flag, node_flag),
        )
    np.fill_diagonal(eligible, False)
    return eligible


def _canonical_undirected_pair(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((str(left), str(right)), key=lambda value: (value.casefold(), value)))


def _raw_seed_edge_factor_table(
    project: Path,
    definition: dict[str, Any],
    symbols: list[str],
    multiplier: float,
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> pd.DataFrame:
    """Rescore a small set of seed endpoints without loading full BF catalogs."""
    handler = definition["normalization"]["handler"]
    if len(symbols) < 2:
        return pd.DataFrame(columns=["node_a", "node_b", "bayes_factor"])
    if handler == "mpkccd_localization":
        return _localization_edge_factors(
            project, symbols, multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "kinase_predictor":
        return _kinase_predictor_edge_factors(
            project, symbols, multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "string_v12":
        return _string_edge_factors(
            project, symbols, multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "hpa_primary":
        return _hpa_edge_factors(
            project, symbols, multiplier, high_confidence=False,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "hpa_high_confidence":
        return _hpa_edge_factors(
            project, symbols, multiplier, high_confidence=True,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "omnipath_core":
        return _omnipath_edge_factors(
            project, symbols, multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    if handler == "stitch_secondary_messenger":
        return _stitch_edge_factors(
            project, symbols, multiplier,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
    raise ValueError(f"unsupported edge calibration handler: {handler}")


def calibrate_edge_parameters(
    project: Path,
    registry: dict[str, Any],
    config: dict[str, Any],
    graph_symbols: list[str],
    graph_metadata: pd.DataFrame,
    known_edges: list[tuple[str, str]],
    *,
    progress: Callable[[str, float], None] | None = None,
    check_cancel: Callable[[], None] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    """Fit primary edge-stream settings against known undirected relationships."""
    if not known_edges:
        raise ValueError("edge calibration requires at least one known edge")
    target_symbols = list(
        dict.fromkeys(symbol for pair in known_edges for symbol in pair)
    )
    target_index = {symbol: index for index, symbol in enumerate(target_symbols)}
    pair_keys = [_canonical_undirected_pair(*pair) for pair in known_edges]
    seed_set = set(
        pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", usecols=["symbol"])[
            "symbol"
        ].astype(str)
    )
    seed_targets = [symbol for symbol in target_symbols if symbol in seed_set]
    incremental_targets = [
        pair for pair in known_edges if pair[0] not in seed_set or pair[1] not in seed_set
    ]
    definitions = {item["id"]: item for item in registry["edge_streams"]}
    active = {
        stream_id: definitions[stream_id]
        for stream_id, state in config["edge_streams"].items()
        if state["enabled"] and not definitions[stream_id].get("derived")
    }
    metadata_subset = graph_metadata.loc[
        graph_metadata["symbol"].astype(str).isin(target_symbols)
    ].copy()
    eligibility: dict[str, np.ndarray] = {}
    for stream_id, definition in active.items():
        matrix = edge_stream_eligibility_matrix(
            project,
            definition,
            target_symbols,
            graph_metadata=metadata_subset,
        )
        eligibility[stream_id] = np.asarray(
            [
                matrix[target_index[left], target_index[right]]
                for left, right in known_edges
            ],
            dtype=bool,
        )
    factor_cache: dict[tuple[str, float], np.ndarray] = {}
    penalize = bool(config["edge_integration"]["penalize_unsupported"])
    unsupported_factor = float(
        config["edge_integration"]["unsupported_bayes_factor"]
    )
    global_continuous_negative = bool(
        config["edge_integration"]["continuous_negative_evidence"]
    )
    continuous_floor = float(
        config["edge_integration"]["continuous_bayes_factor_floor"]
    )

    def stream_factors(stream_id: str, multiplier: float) -> np.ndarray:
        cache_key = (stream_id, round(float(multiplier), 12))
        cached = factor_cache.get(cache_key)
        if cached is not None:
            return cached
        definition = active[stream_id]
        stream_continuous_negative = (
            global_continuous_negative
            or bool(
                config["edge_streams"][stream_id].get(
                    "continuous_negative_evidence", False
                )
            )
        )
        tables = [
            _raw_seed_edge_factor_table(
                project,
                definition,
                seed_targets,
                multiplier,
                continuous_negative=stream_continuous_negative,
                minimum_bayes_factor=continuous_floor,
            )
        ]
        if incremental_targets:
            tables.append(
                incremental_pair_factor_table(
                    project,
                    definition["normalization"]["handler"],
                    multiplier,
                    incremental_targets,
                    continuous_negative=stream_continuous_negative,
                    minimum_bayes_factor=continuous_floor,
                )
            )
        lookup: dict[tuple[str, str], float] = {}
        for table in tables:
            for row in table.itertuples(index=False):
                key = _canonical_undirected_pair(row.node_a, row.node_b)
                factor = float(row.bayes_factor)
                lookup[key] = lookup.get(key, 1.0) * factor
        factors = np.ones(len(known_edges), dtype=float)
        observed = np.zeros(len(known_edges), dtype=bool)
        for index, key in enumerate(pair_keys):
            if key in lookup:
                factors[index] = lookup[key]
                observed[index] = True
        if stream_continuous_negative and definition.get("negative_evidence"):
            factors[eligibility[stream_id] & ~observed] = continuous_floor
        elif penalize and definition.get("negative_evidence"):
            factors[eligibility[stream_id] & ~observed] = unsupported_factor
        if (~np.isfinite(factors) | (factors <= 0)).any():
            raise ValueError(
                f"edge calibration stream {stream_id} produced invalid factors"
            )
        factor_cache[cache_key] = factors
        return factors

    specs = _calibration_parameter_specs(registry, config, "edge_streams")
    prior = float(config["edge_integration"]["prior_probability"])
    prior_log_odds = math.log(prior / (1.0 - prior))

    def probabilities(settings: dict[str, float]) -> np.ndarray:
        log_odds = np.full(len(known_edges), prior_log_odds, dtype=float)
        for stream_id in active:
            multiplier = settings[f"{stream_id}:tq_multiplier"]
            weight = settings[f"{stream_id}:weight"]
            log_odds += weight * np.log(stream_factors(stream_id, multiplier))
        return stable_expit(log_odds)

    summary, trace = bounded_powell_positive_calibration(
        specs,
        probabilities,
        target_count=len(known_edges),
        regularization_strength=config["calibration"]["regularization_strength"],
        multistart_count=config["calibration"]["multistart_count"],
        progress=progress,
        check_cancel=check_cancel,
        stage_label="edge evidence",
    )
    _apply_calibrated_parameters(config, "edge_streams", summary["parameters"])
    target_table = pd.DataFrame(
        {
            "node_a": [pair[0] for pair in known_edges],
            "node_b": [pair[1] for pair in known_edges],
            "initial_probability": summary.pop("initial_probabilities"),
            "fitted_probability": summary.pop("final_probabilities"),
        }
    )
    summary["known_edges"] = [list(pair) for pair in known_edges]
    summary["derived_streams_excluded"] = [
        stream_id
        for stream_id, state in config["edge_streams"].items()
        if state["enabled"] and definitions[stream_id].get("derived")
    ]
    summary["scope_note"] = (
        "Only known-present undirected edge controls entered the fit loss. Unknown "
        "pairs were unlabeled. Derived scaffold closure was held fixed and excluded "
        "from optimization."
    )
    return summary, target_table, trace


def scaffold_triadic_closure_factors(
    project: Path,
    graph_symbols: list[str],
    preclosure_probabilities: np.ndarray,
    state: dict[str, Any],
    *,
    graph_metadata: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Assign one fixed factor to every protein pair sharing a strong scaffold.

    Only protein nodes can be endpoints or common scaffolds. A common node must
    contain the exact ``adaptor_scaffold`` class token. A pair qualifies when
    both of its protein-scaffold probabilities are strictly above the selected
    anchor cutoff for at least one shared scaffold. Every qualifying pair gets
    the same configured support likelihood relative to the neutral likelihood;
    scaffold degree and the number of shared scaffolds do not change the factor.
    The operation is one pass, so inferred closure edges never become anchors.
    """
    symbols = list(graph_symbols)
    metadata = _metadata_for_graph(project, symbols, graph_metadata)
    node_types = metadata.get("node_type", pd.Series("", index=metadata.index))
    classes = metadata.get("classes", pd.Series("", index=metadata.index))
    protein_indices = np.flatnonzero(
        node_types.astype(str).str.casefold().eq("protein").to_numpy()
    )
    scaffold_indices = np.asarray(
        [
            index
            for index in protein_indices
            if "adaptor_scaffold"
            in {
                token.strip()
                for token in str(classes.iloc[index]).split(";")
                if token.strip()
            }
        ],
        dtype=int,
    )
    columns = [
        "node_a",
        "node_b",
        "preclosure_probability",
        "bayes_factor",
        "closure_support_likelihood",
        "supporting_scaffold_count",
        "supporting_scaffolds",
        "anchor_probability_cutoff",
        "closure_rule",
    ]
    anchor_cutoff = float(state.get("parameters", {}).get(
        "anchor_probability_cutoff", 0.9
    ))
    closure_likelihood = float(state.get("parameters", {}).get(
        "closure_likelihood", 0.9
    ))
    if not 0.5 <= closure_likelihood < 1.0:
        raise ValueError("closure likelihood must be at least 0.5 and below 1.0")
    closure_factor = closure_likelihood / NEUTRAL_LIKELIHOOD
    rule_text = (
        "Both protein-scaffold probabilities must be strictly above the anchor "
        "cutoff for at least one shared adaptor_scaffold; every qualifying pair "
        "receives the same factor in one non-recursive pass."
    )
    empty_summary = {
        "scaffold_node_count": int(len(scaffold_indices)),
        "usable_scaffold_count": 0,
        "anchor_probability_cutoff_exclusive": anchor_cutoff,
        "closure_support_likelihood": closure_likelihood,
        "closure_bayes_factor": closure_factor,
        "anchored_protein_scaffold_associations": 0,
        "pairs_sharing_at_least_one_scaffold": 0,
        "pairs_with_non_neutral_closure_factor": 0,
        "degree_adjustment": "none",
        "empirical_tq": None,
        "closure_rule": rule_text,
    }
    if len(protein_indices) < 2 or not len(scaffold_indices):
        return pd.DataFrame(columns=columns), empty_summary

    association = preclosure_probabilities[np.ix_(protein_indices, scaffold_indices)]
    anchors = association > anchor_cutoff
    partner_counts = anchors.sum(axis=0)
    usable = np.flatnonzero(partner_counts >= 2)
    if not len(usable):
        return pd.DataFrame(columns=columns), empty_summary

    usable_anchors = anchors[:, usable].astype(np.uint16)
    support_count = usable_anchors @ usable_anchors.T
    n_proteins = len(protein_indices)
    upper_indices = np.triu_indices(n_proteins, 1)
    qualifying = support_count[upper_indices] > 0
    rows_i = upper_indices[0][qualifying]
    rows_j = upper_indices[1][qualifying]
    protein_symbols = np.asarray(symbols, dtype=object)[protein_indices]
    scaffold_symbols = np.asarray(symbols, dtype=object)[scaffold_indices]
    supporting_lists: list[str] = []
    for left, right in zip(rows_i, rows_j, strict=True):
        support = np.flatnonzero(anchors[left] & anchors[right])
        supporting_lists.append(";".join(scaffold_symbols[support].tolist()))
    audit = pd.DataFrame(
        {
            "node_a": protein_symbols[rows_i],
            "node_b": protein_symbols[rows_j],
            "preclosure_probability": preclosure_probabilities[
                protein_indices[rows_i], protein_indices[rows_j]
            ],
            "bayes_factor": np.full(len(rows_i), closure_factor),
            "closure_support_likelihood": np.full(
                len(rows_i), closure_likelihood
            ),
            "supporting_scaffold_count": support_count[rows_i, rows_j],
            "supporting_scaffolds": supporting_lists,
            "anchor_probability_cutoff": anchor_cutoff,
            "closure_rule": rule_text,
        }
    )
    audit = audit.sort_values(
        ["supporting_scaffold_count", "node_a", "node_b"],
        ascending=[False, True, True],
        kind="stable",
    ).reset_index(drop=True)
    summary = {
        **empty_summary,
        "usable_scaffold_count": int(len(usable)),
        "anchored_protein_scaffold_associations": int(anchors.sum()),
        "pairs_sharing_at_least_one_scaffold": int(len(rows_i)),
        "pairs_with_non_neutral_closure_factor": int(
            len(rows_i) if closure_factor > 1.0 + FACTOR_EPSILON else 0
        ),
        "maximum_closure_factor": closure_factor,
    }
    return audit, summary


def combine_edge_factors(
    project: Path,
    registry: dict[str, Any],
    config: dict[str, Any],
    graph_symbols: list[str],
    *,
    graph_metadata: pd.DataFrame | None = None,
    audit_collector: dict[str, pd.DataFrame] | None = None,
    contribution_collector: dict[str, np.ndarray] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    prior = config["edge_integration"]["prior_probability"]
    penalize_unsupported = config["edge_integration"]["penalize_unsupported"]
    unsupported_factor = config["edge_integration"]["unsupported_bayes_factor"]
    global_continuous_negative = bool(
        config["edge_integration"]["continuous_negative_evidence"]
    )
    continuous_floor = float(
        config["edge_integration"]["continuous_bayes_factor_floor"]
    )
    base_log_odds = math.log(prior / (1.0 - prior))
    log_odds = np.full((len(graph_symbols), len(graph_symbols)), base_log_odds, dtype=float)
    unique_pair_indices = np.triu_indices(len(graph_symbols), 1)
    index = {symbol: position for position, symbol in enumerate(graph_symbols)}
    seed_universe = pd.read_csv(
        project / UNIVERSE_RELATIVE, sep="\t", dtype=str
    ).fillna("")
    seed_set = set(seed_universe["symbol"].astype(str))
    seed_graph_symbols = [symbol for symbol in graph_symbols if symbol in seed_set]
    has_incremental_nodes = len(seed_graph_symbols) != len(graph_symbols)
    incremental_pairs: list[tuple[str, str]] = []
    if has_incremental_nodes:
        incremental_set = set(graph_symbols).difference(seed_set)
        incremental_pairs = [
            (left, right)
            for left_index, left in enumerate(graph_symbols)
            for right in graph_symbols[left_index + 1 :]
            if left in incremental_set or right in incremental_set
        ]
    stream_defs = {stream["id"]: stream for stream in registry["edge_streams"]}
    active_streams: list[dict[str, Any]] = []
    total_penalty_applications = 0
    derived_streams: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    for stream_id, state in config["edge_streams"].items():
        if not state["enabled"] or state["weight"] <= 0:
            continue
        definition = stream_defs[stream_id]
        if definition.get("derived"):
            derived_streams.append((stream_id, definition, state))
            continue
        stream_continuous_negative = (
            global_continuous_negative
            or bool(state.get("continuous_negative_evidence", False))
        )
        seed_table = edge_stream_factor_table(
            project,
            definition,
            state,
            seed_graph_symbols,
            include_negative=penalize_unsupported,
            continuous_negative=stream_continuous_negative,
            minimum_bayes_factor=continuous_floor,
        )
        tables = [seed_table]
        if has_incremental_nodes:
            incremental_loader = (
                incremental_pair_factor_table
                if len(incremental_pairs) <= 500_000
                else incremental_factor_table
            )
            tables.append(
                incremental_loader(
                    project,
                    definition["normalization"]["handler"],
                    _stream_multiplier(state),
                    (
                        incremental_pairs
                        if incremental_loader is incremental_pair_factor_table
                        else graph_symbols
                    ),
                    continuous_negative=stream_continuous_negative,
                    minimum_bayes_factor=continuous_floor,
                )
            )
        table = pd.concat(tables, ignore_index=True)
        left = table["node_a"].map(index)
        right = table["node_b"].map(index)
        mask = left.notna() & right.notna()
        values = pd.to_numeric(
            table.loc[mask, "bayes_factor"], errors="raise"
        ).to_numpy(float)
        if (values <= 0).any():
            raise ValueError(f"edge stream {stream_id} contains non-positive factors")
        weighted_logs = float(state["weight"]) * np.log(values)
        left_index = left[mask].astype(int).to_numpy()
        right_index = right[mask].astype(int).to_numpy()
        stream_log_odds = np.zeros_like(log_odds)
        np.add.at(stream_log_odds, (left_index, right_index), weighted_logs)
        np.add.at(stream_log_odds, (right_index, left_index), weighted_logs)
        observed_pairs = np.zeros_like(log_odds, dtype=bool)
        observed_pairs[left_index, right_index] = True
        observed_pairs[right_index, left_index] = True
        eligible_pairs = edge_stream_eligibility_matrix(
            project,
            definition,
            graph_symbols,
            graph_metadata=graph_metadata,
        )
        eligible_upper = np.triu(eligible_pairs, 1)
        unsupported_upper = eligible_upper & ~observed_pairs
        unsupported_left, unsupported_right = np.where(unsupported_upper)
        penalties_applied = 0
        negative_fill_factor: float | None = None
        if stream_continuous_negative and definition.get("negative_evidence"):
            negative_fill_factor = continuous_floor
        elif penalize_unsupported and definition.get("negative_evidence"):
            negative_fill_factor = unsupported_factor
        if negative_fill_factor is not None:
            penalty_log = float(state["weight"]) * math.log(negative_fill_factor)
            np.add.at(
                stream_log_odds,
                (unsupported_left, unsupported_right),
                penalty_log,
            )
            np.add.at(
                stream_log_odds,
                (unsupported_right, unsupported_left),
                penalty_log,
            )
            penalties_applied = int(len(unsupported_left))
            total_penalty_applications += penalties_applied
        log_odds += stream_log_odds
        if contribution_collector is not None:
            contribution_collector[stream_id] = stream_log_odds
        active_streams.append(
            {
                "id": stream_id,
                "label": definition["label"],
                "weight": float(state["weight"]),
                "tq_multiplier": _stream_multiplier(state),
                "normalization_reference": definition["normalization"]["reference"],
                "supported_pairs_in_selected_graph": int(mask.sum()),
                "eligible_pairs_for_negative_evidence": int(eligible_upper.sum()),
                "unsupported_eligible_pairs": int(unsupported_upper.sum()),
                "negative_penalties_applied": penalties_applied,
                "continuous_zero_imputations": (
                    penalties_applied if stream_continuous_negative else 0
                ),
                "continuous_negative_evidence": stream_continuous_negative,
                "explicit_source_factors_below_one": int((values < 1.0).sum()),
                "unsupported_pair_bayes_factor": (
                    negative_fill_factor if penalties_applied else None
                ),
                "negative_evidence_policy": definition.get("negative_evidence"),
                "factor_distribution": evidence_factor_distribution_summary(
                    np.exp(
                        np.clip(
                            stream_log_odds[unique_pair_indices] / float(state["weight"]),
                            -700.0,
                            700.0,
                        )
                    ),
                    distribution_scope="all unique undirected graph pairs",
                ),
            }
        )
    preclosure_probabilities = stable_expit(log_odds)
    np.fill_diagonal(preclosure_probabilities, 0.0)
    for stream_id, definition, state in derived_streams:
        handler = definition["normalization"]["handler"]
        if handler != "scaffold_triadic_closure":
            raise ValueError(f"unsupported derived edge handler: {handler}")
        table, closure_summary = scaffold_triadic_closure_factors(
            project,
            graph_symbols,
            preclosure_probabilities,
            state,
            graph_metadata=graph_metadata,
        )
        left_index = table["node_a"].map(index).astype(int).to_numpy()
        right_index = table["node_b"].map(index).astype(int).to_numpy()
        values = table["bayes_factor"].to_numpy(float)
        weighted_logs = float(state["weight"]) * np.log(values)
        np.add.at(log_odds, (left_index, right_index), weighted_logs)
        np.add.at(log_odds, (right_index, left_index), weighted_logs)
        derived_log_odds = np.zeros_like(log_odds)
        np.add.at(
            derived_log_odds,
            (left_index, right_index),
            weighted_logs,
        )
        np.add.at(
            derived_log_odds,
            (right_index, left_index),
            weighted_logs,
        )
        if contribution_collector is not None:
            contribution_collector[stream_id] = derived_log_odds
        if len(table):
            pair_pre = np.clip(
                table["preclosure_probability"].to_numpy(float),
                FACTOR_EPSILON,
                1.0 - FACTOR_EPSILON,
            )
            pair_posteriors = stable_expit(
                np.log(pair_pre / (1.0 - pair_pre)) + weighted_logs
            )
        else:
            pair_posteriors = np.asarray([], dtype=float)
        table["postclosure_probability"] = pair_posteriors
        cutoff = config["edge_integration"]["output_probability_cutoff"]
        closure_summary["newly_above_output_cutoff"] = int(
            (
                (table["preclosure_probability"] <= cutoff)
                & (table["postclosure_probability"] > cutoff)
            ).sum()
        )
        if audit_collector is not None:
            audit_collector[stream_id] = table
        active_streams.append(
            {
                "id": stream_id,
                "label": definition["label"],
                "weight": float(state["weight"]),
                "normalization": "binary_fixed_likelihood",
                "normalization_reference": definition["normalization"]["reference"],
                "supported_pairs_in_selected_graph": int(len(table)),
                "derived": True,
                "parameters": dict(state.get("parameters", {})),
                "derivation_summary": closure_summary,
                "factor_distribution": evidence_factor_distribution_summary(
                    np.exp(
                        np.clip(
                            derived_log_odds[unique_pair_indices]
                            / float(state["weight"]),
                            -700.0,
                            700.0,
                        )
                    ),
                    distribution_scope="all unique undirected graph pairs",
                ),
            }
        )
    probabilities = stable_expit(log_odds)
    np.fill_diagonal(probabilities, 0.0)
    matrix = pd.DataFrame(probabilities, index=graph_symbols, columns=graph_symbols)
    matrix.index.name = "symbol"
    upper = probabilities[unique_pair_indices]
    cutoff = config["edge_integration"]["output_probability_cutoff"]
    summary = {
        "edge_prior_probability": prior,
        "penalize_unsupported": penalize_unsupported,
        "unsupported_edge_bayes_factor": unsupported_factor,
        "continuous_negative_evidence": global_continuous_negative,
        "continuous_negative_streams": [
            stream["id"]
            for stream in active_streams
            if stream.get("continuous_negative_evidence")
        ],
        "continuous_bayes_factor_floor": continuous_floor,
        "negative_penalty_applications": total_penalty_applications,
        "unique_pair_count": int(len(upper)),
        "pairs_above_output_cutoff": int((upper > cutoff).sum()),
        "pairs_below_prior": int((upper < prior - FACTOR_EPSILON).sum()),
        "pairs_at_exact_prior": int(np.isclose(upper, prior, atol=1e-12, rtol=0).sum()),
        "output_probability_cutoff_exclusive": cutoff,
        "probability_distribution": probability_distribution_summary(
            upper,
            prior_probability=prior,
            output_cutoff=cutoff,
        ),
        "active_streams": active_streams,
        "integration_rule": (
            "Edge prior odds are multiplied by each selected Bayes factor raised "
            "to its user-specified weight. Fixed-absence mode gives a source-eligible "
            "pair with no record one configured BF below 1. Alternatively, continuous "
            "negative mode removes the 0.5 likelihood floor, retains quantitative "
            "factors below 1, and scores eligible no-record pairs as x=0 using the "
            "numerical BF floor. Out-of-scope pairs remain at BF=1. When both negative "
            "modes are off, missing sparse-table entries receive BF=1. "
            "If enabled, scaffold closure is derived once from the pre-closure graph "
            "and appended without recursive feedback or absence penalties."
        ),
    }
    return matrix, summary


def load_or_build_external_target_vector(
    project: Path,
    target_input: str,
    target_symbol: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    cached = cached_target_extension(project, target_input, target_symbol)
    if cached is None:
        result = build_target_adjacency_vector(
            target_input,
            project_root=project,
            write_outputs=True,
            write_extended_matrix=True,
        )
        return result.vector.copy(), result.summary
    _, matrix_path, summary = cached
    vector_path = matrix_path.parent / "target_adjacency_vector.tsv"
    return pd.read_csv(vector_path, sep="\t"), summary


def _truth_array(values: pd.Series) -> np.ndarray:
    if values.dtype == bool:
        return values.to_numpy(bool)
    return values.fillna(False).astype(str).str.casefold().eq("true").to_numpy()


def external_target_stream_values(
    project: Path,
    target_symbol: str,
    vector_by_symbol: pd.DataFrame,
    symbols: pd.Index,
    definition: dict[str, Any],
    state: dict[str, Any],
    *,
    continuous_negative: bool = False,
    minimum_bayes_factor: float = 1e-6,
) -> np.ndarray | None:
    multiplier = _stream_multiplier(state)
    handler = definition["normalization"]["handler"]
    if handler == "stitch_secondary_messenger":
        result = np.ones(len(symbols), dtype=float)
        if "string_target_id" not in vector_by_symbol.columns:
            return result
        target_ids = (
            vector_by_symbol["string_target_id"].fillna("").astype(str).str.strip()
        )
        target_ids = target_ids.loc[target_ids.ne("")].unique()
        if len(target_ids) != 1:
            return result
        evidence = pd.read_csv(
            project / STITCH_ALL_MOUSE_EDGES_RELATIVE,
            sep="\t",
            compression="gzip",
        )
        evidence = evidence.loc[evidence["protein"].eq(target_ids[0])].copy()
        if evidence.empty:
            return result
        reference = STITCH_REFERENCE_SCORE * multiplier
        evidence["bayes_factor"] = _stitch_bayes_factors(
            pd.to_numeric(evidence["stitch_score"], errors="raise").to_numpy(float),
            reference,
            continuous_negative=continuous_negative,
            minimum_bayes_factor=minimum_bayes_factor,
        )
        by_messenger = evidence.groupby("messenger_node")["bayes_factor"].max()
        mapped = symbols.to_series().map(by_messenger)
        found = mapped.notna().to_numpy()
        result[found] = mapped[found].to_numpy(float)
        return result
    column = definition.get("target_factor_column")
    if (
        math.isclose(multiplier, 1.0, rel_tol=0.0, abs_tol=1e-15)
        and not continuous_negative
    ):
        if not column:
            return None
        return pd.to_numeric(
            vector_by_symbol.reindex(symbols)[column], errors="coerce"
        ).fillna(1.0).to_numpy(float)

    aligned = vector_by_symbol.reindex(symbols)
    result = np.ones(len(aligned), dtype=float)
    if handler == "mpkccd_localization":
        observed = _truth_array(aligned["mpkccd_localization_observed_for_both"])
        dot = pd.to_numeric(aligned["mpkccd_dot_product"], errors="coerce").to_numpy(float)
        target_tq = pd.to_numeric(aligned["mpkccd_target_Tq"], errors="coerce").to_numpy(float)
        node_tq = pd.to_numeric(aligned["mpkccd_node_Tq"], errors="coerce").to_numpy(float)
        valid = observed & np.isfinite(dot) & np.isfinite(target_tq) & np.isfinite(node_tq)
        if np.any(valid):
            if continuous_negative:
                target_factor = continuous_complement_bayes_factor(
                    dot[valid],
                    target_tq[valid] * multiplier,
                    minimum_bayes_factor=minimum_bayes_factor,
                )
                node_factor = continuous_complement_bayes_factor(
                    dot[valid],
                    node_tq[valid] * multiplier,
                    minimum_bayes_factor=minimum_bayes_factor,
                )
                result[valid] = (target_factor + node_factor) / 2.0
            else:
                target_factor = complement_minimum_likelihood(
                    dot[valid], target_tq[valid] * multiplier
                )
                node_factor = complement_minimum_likelihood(
                    dot[valid], node_tq[valid] * multiplier
                )
                result[valid] = ((target_factor + node_factor) / 2.0) / NEUTRAL_LIKELIHOOD
        return result
    if handler == "kinase_predictor":
        prediction_path = (
            project / "results/path_finding/target_extensions" / safe_name(target_symbol)
            / "target_kinase_predictions.tsv"
        )
        if not prediction_path.exists():
            return result
        predictions = pd.read_csv(prediction_path, sep="\t")
        raw_score = pd.to_numeric(predictions["raw_score"], errors="raise").to_numpy(float)
        tq = pd.to_numeric(predictions["site_Tq_q75"], errors="raise").to_numpy(float)
        if continuous_negative:
            site_bf = continuous_complement_bayes_factor(
                np.maximum(raw_score, 0.0),
                tq * multiplier,
                minimum_bayes_factor=minimum_bayes_factor,
            )
        else:
            site_bf = complement_minimum_likelihood(
                np.maximum(raw_score, 0.0), tq * multiplier
            ) / NEUTRAL_LIKELIHOOD
        predictions["gui_log_bf"] = np.log(site_bf)
        grouped = predictions.groupby("kinase_node")["gui_log_bf"].sum()
        mapped = aligned.index.to_series().map(grouped)
        found = mapped.notna().to_numpy()
        result[found] = np.exp(np.clip(mapped[found].to_numpy(float), -700.0, 700.0))
        if continuous_negative:
            result[found] = np.maximum(result[found], minimum_bayes_factor)
        return result
    if handler == "string_v12":
        score = pd.to_numeric(aligned["string_combined_score"], errors="coerce").to_numpy(float)
        reference = STRING_REFERENCE_SCORE * multiplier
        if not 0 < reference < 1:
            raise ValueError("STRING reference score after scaling must be between 0 and 1")
        observed = np.isfinite(score) & (score > 0) & (score < 1)
        result[observed] = (
            score[observed] / (1.0 - score[observed])
        ) / (reference / (1.0 - reference))
        if continuous_negative:
            result[observed] = np.maximum(
                result[observed], minimum_bayes_factor
            )
        return result
    if handler == "hpa_primary":
        observed = _truth_array(aligned["hpa_profiles_observed_for_both"])
        similarity = pd.to_numeric(aligned["hpa_cosine_similarity"], errors="coerce").to_numpy(float)
        target_tq = pd.to_numeric(aligned["hpa_target_Tq"], errors="coerce").to_numpy(float)
        node_tq = pd.to_numeric(aligned["hpa_node_Tq"], errors="coerce").to_numpy(float)
        valid = observed & np.isfinite(similarity) & np.isfinite(target_tq) & np.isfinite(node_tq)
        if np.any(valid):
            if continuous_negative:
                target_factor = continuous_complement_bayes_factor(
                    similarity[valid],
                    target_tq[valid] * multiplier,
                    minimum_bayes_factor=minimum_bayes_factor,
                )
                node_factor = continuous_complement_bayes_factor(
                    similarity[valid],
                    node_tq[valid] * multiplier,
                    minimum_bayes_factor=minimum_bayes_factor,
                )
                result[valid] = (target_factor + node_factor) / 2.0
            else:
                target_factor = complement_minimum_likelihood(
                    similarity[valid], target_tq[valid] * multiplier
                )
                node_factor = complement_minimum_likelihood(
                    similarity[valid], node_tq[valid] * multiplier
                )
                result[valid] = ((target_factor + node_factor) / 2.0) / NEUTRAL_LIKELIHOOD
        return result
    if handler == "hpa_high_confidence":
        return None
    if handler == "omnipath_core":
        effort = pd.to_numeric(aligned["omnipath_curation_effort"], errors="coerce").to_numpy(float)
        observed = np.isfinite(effort) & (effort > 0)
        if continuous_negative:
            result[observed] = continuous_complement_bayes_factor(
                effort[observed],
                6.0 * multiplier,
                minimum_bayes_factor=minimum_bayes_factor,
            )
        else:
            support = 1.0 - np.exp(-0.5 * np.square(effort[observed] / (6.0 * multiplier)))
            likelihood = NEUTRAL_LIKELIHOOD + (1.0 - NEUTRAL_LIKELIHOOD) * support
            result[observed] = likelihood / NEUTRAL_LIKELIHOOD
        return result
    raise ValueError(f"unsupported target normalization handler: {handler}")


def external_target_stream_eligibility(
    project: Path,
    target_symbol: str,
    vector_by_symbol: pd.DataFrame,
    symbols: pd.Index,
    definition: dict[str, Any],
) -> np.ndarray:
    """Return source scope for edges between one external target and the graph."""
    aligned = vector_by_symbol.reindex(symbols)
    handler = definition["normalization"]["handler"]
    if handler == "mpkccd_localization":
        return _truth_array(aligned["mpkccd_localization_observed_for_both"])
    if handler == "kinase_predictor":
        prediction_path = (
            project / "results/path_finding/target_extensions" / safe_name(target_symbol)
            / "target_kinase_predictions.tsv"
        )
        if not prediction_path.exists() or pd.read_csv(prediction_path, sep="\t").empty:
            return np.zeros(len(symbols), dtype=bool)
        classes = aligned["classes"].fillna("").astype(str)
        return classes.map(
            lambda value: "kinase"
            in {token.strip() for token in value.split(";") if token.strip()}
        ).to_numpy(bool)
    if handler == "string_v12":
        target_mapped = bool(
            vector_by_symbol["string_target_id"]
            .fillna("")
            .astype(str)
            .str.strip()
            .ne("")
            .any()
        )
        mapping = pd.read_csv(
            project / STRING_MAPPING_RELATIVE, sep="\t", dtype=str
        ).fillna("")
        mapped_symbols = set(
            mapping.loc[
                mapping["mapping_status"].eq("mapped")
                & mapping["string_id"].str.strip().ne(""),
                "symbol",
            ].astype(str)
        )
        node_mapped = symbols.isin(mapped_symbols).copy()
        cached = incremental_node_scope_table(project, symbols.astype(str).tolist())
        if not cached.empty:
            cached_mapped = set(
                cached.loc[cached["string_mapped"].astype(bool), "symbol"].astype(str)
            )
            node_mapped |= symbols.isin(cached_mapped)
        protein = aligned["node_type"].fillna("").astype(str).str.casefold().eq("protein")
        return target_mapped & node_mapped & protein.to_numpy(bool)
    if handler == "hpa_primary":
        return _truth_array(aligned["hpa_profiles_observed_for_both"])
    if handler == "hpa_high_confidence":
        return np.zeros(len(symbols), dtype=bool)
    if handler == "omnipath_core":
        return aligned["node_type"].fillna("").astype(str).str.casefold().eq("protein").to_numpy(bool)
    if handler == "stitch_secondary_messenger":
        target_mapped = bool(
            vector_by_symbol["string_target_id"]
            .fillna("")
            .astype(str)
            .str.strip()
            .ne("")
            .any()
        )
        evidence = pd.read_csv(
            project / STITCH_ALL_MOUSE_EDGES_RELATIVE,
            sep="\t",
            compression="gzip",
            usecols=["messenger_node"],
        )
        messengers = set(evidence["messenger_node"].astype(str))
        return target_mapped & symbols.isin(messengers)
    return np.zeros(len(symbols), dtype=bool)


def append_external_target(
    project: Path,
    matrix: pd.DataFrame,
    target_symbol: str,
    target_vector: pd.DataFrame,
    registry: dict[str, Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, list[str]]:
    prior = config["edge_integration"]["prior_probability"]
    penalize_unsupported = config["edge_integration"]["penalize_unsupported"]
    unsupported_factor = config["edge_integration"]["unsupported_bayes_factor"]
    global_continuous_negative = bool(
        config["edge_integration"]["continuous_negative_evidence"]
    )
    continuous_floor = float(
        config["edge_integration"]["continuous_bayes_factor_floor"]
    )
    target_log_odds = np.full(len(matrix), math.log(prior / (1.0 - prior)), dtype=float)
    vector_by_symbol = target_vector.set_index("symbol")
    warnings: list[str] = []
    stream_defs = {stream["id"]: stream for stream in registry["edge_streams"]}
    for stream_id, state in config["edge_streams"].items():
        if not state["enabled"] or state["weight"] <= 0:
            continue
        definition = stream_defs[stream_id]
        if definition.get("derived"):
            warnings.append(
                f"{definition['label']} is derived from the internal graph and was "
                "not reapplied to the external target."
            )
            continue
        stream_continuous_negative = (
            global_continuous_negative
            or bool(state.get("continuous_negative_evidence", False))
        )
        values = external_target_stream_values(
            project,
            target_symbol,
            vector_by_symbol,
            matrix.index,
            definition,
            state,
            continuous_negative=stream_continuous_negative,
            minimum_bayes_factor=continuous_floor,
        )
        if values is None:
            warnings.append(
                f"{definition['label']} has no external-target factor implementation; "
                "it was neutral for target edges."
            )
            continue
        if (values <= 0).any() or not np.isfinite(values).all():
            raise ValueError(f"external-target stream {stream_id} produced invalid factors")
        if (
            (stream_continuous_negative or penalize_unsupported)
            and definition.get("negative_evidence")
        ):
            eligible = external_target_stream_eligibility(
                project,
                target_symbol,
                vector_by_symbol,
                matrix.index,
                definition,
            )
            unsupported = eligible & np.isclose(
                values, 1.0, atol=FACTOR_EPSILON, rtol=0.0
            )
            values = values.copy()
            fill_factor = (
                continuous_floor
                if stream_continuous_negative
                else unsupported_factor
            )
            values[unsupported] = fill_factor
            if np.any(unsupported):
                warnings.append(
                    f"{definition['label']} applied "
                    f"{'continuous x=0' if stream_continuous_negative else 'unsupported-pair'} "
                    f"BF {fill_factor:g} to {int(unsupported.sum())} eligible "
                    "external-target edges."
                )
        target_log_odds += float(state["weight"]) * np.log(values)
    target_probabilities = stable_expit(target_log_odds)
    symbols = [*matrix.index.astype(str), target_symbol]
    values = np.zeros((len(symbols), len(symbols)), dtype=float)
    values[:-1, :-1] = matrix.to_numpy(float)
    values[-1, :-1] = target_probabilities
    values[:-1, -1] = target_probabilities
    extended = pd.DataFrame(values, index=symbols, columns=symbols)
    extended.index.name = "symbol"
    return extended, warnings


def supported_edge_table(matrix: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    symbols = matrix.index.astype(str).to_numpy()
    values = matrix.to_numpy(float)
    left, right = np.triu_indices(len(symbols), 1)
    probabilities = values[left, right]
    keep = probabilities > cutoff
    result = pd.DataFrame(
        {
            "node_a": symbols[left[keep]],
            "node_b": symbols[right[keep]],
            "edge_probability": probabilities[keep],
        }
    )
    return result.sort_values("edge_probability", ascending=False).reset_index(drop=True)


def run_paths(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    config: dict[str, Any],
    start_input: str,
    target_symbol: str,
    *,
    directed_graph: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    symbols = matrix.index.astype(str).tolist()
    resolved = resolve_existing_node(
        start_input,
        metadata,
        symbols,
        allow_curated_alias=True,
    )
    if resolved is None:
        raise ValueError(
            f"Starting node {start_input!r} was not selected by the chosen node streams"
        )
    start_symbol, resolution_method = resolved
    if target_symbol not in symbols:
        raise ValueError("target is absent from the constructed graph")
    if start_symbol == target_symbol:
        raise ValueError("start and target resolve to the same node")
    path_config = config["path"]
    eligibility = build_intermediate_eligibility(
        metadata,
        symbols,
        start_symbol,
        target_symbol,
        signaling_intermediates_only=path_config["signaling_intermediates_only"],
        relay_classes=path_config["allowed_intermediate_classes"],
        exclude_any_scaffold=path_config["exclude_multirole_scaffolds"],
    )
    permitted = frozenset(
        eligibility.loc[eligibility["permitted_in_search"], "matrix_index"].astype(int)
    )
    eligible = frozenset(
        eligibility.loc[eligibility["allowed_as_intermediate"], "matrix_index"].astype(int)
    )
    adjacency, retained_edges = build_adjacency(
        matrix,
        path_config["minimum_edge_probability"],
        directed=directed_graph,
    )
    symbol_index = {symbol: index for index, symbol in enumerate(symbols)}
    start_index = symbol_index[start_symbol]
    target_index = symbol_index[target_symbol]
    ranked = k_shortest_simple_paths(
        adjacency,
        matrix.to_numpy(float),
        start_index,
        target_index,
        top_k=path_config["top_k"],
        max_hops=path_config["max_hops"],
        permitted_nodes=permitted,
    )
    paths, edges = serialize_paths(ranked, matrix, metadata)
    validation = {
        "all_paths_simple": all(len(item.nodes) == len(set(item.nodes)) for item in ranked),
        "all_paths_within_hop_limit": all(
            len(item.nodes) - 1 <= path_config["max_hops"] for item in ranked
        ),
        "all_internal_nodes_eligible": all(
            all(node in eligible for node in item.nodes[1:-1]) for item in ranked
        ),
        "all_edges_above_path_cutoff": bool(
            edges.empty
            or (edges["edge_probability"] > path_config["minimum_edge_probability"]).all()
        ),
    }
    if not all(validation.values()):
        raise AssertionError(f"path validation failed: {validation}")
    probability_values = matrix.to_numpy(dtype=float)
    upper_triangle = np.triu_indices(len(matrix), k=1)
    retained_unique_edges = int(
        (
            np.maximum(probability_values, probability_values.T)[upper_triangle]
            > path_config["minimum_edge_probability"]
        ).sum()
    )
    start_reachable_nodes = connected_component_size(
        adjacency, start_index, permitted
    )
    summary = {
        "start_symbol": start_symbol,
        "start_resolution_method": resolution_method,
        "target_symbol": target_symbol,
        "paths_found": int(len(paths)),
        "top_k_requested": int(path_config["top_k"]),
        "maximum_hops": int(path_config["max_hops"]),
        "minimum_edge_probability_exclusive": float(
            path_config["minimum_edge_probability"]
        ),
        "signaling_intermediates_only": bool(
            path_config["signaling_intermediates_only"]
        ),
        "exclude_multirole_scaffolds": bool(
            path_config["exclude_multirole_scaffolds"]
        ),
        "allowed_intermediate_classes": list(
            path_config["allowed_intermediate_classes"]
        ),
        "eligible_intermediate_nodes": int(
            eligibility["allowed_as_intermediate"].sum()
        ),
        "retained_unique_edges_before_intermediate_filter": retained_unique_edges,
        "retained_transition_count_before_intermediate_filter": int(retained_edges),
        "path_graph_is_directed": bool(directed_graph),
        "start_reachable_node_count": start_reachable_nodes,
        "start_component_size": (
            None if directed_graph else start_reachable_nodes
        ),
        "validation": validation,
    }
    return paths, edges, eligibility, summary


def _records(frame: pd.DataFrame, columns: list[str], limit: int) -> list[dict[str, Any]]:
    if frame.empty:
        return []
    result = frame.loc[:, columns].head(limit).copy()
    return json.loads(result.to_json(orient="records"))


PATH_NETWORK_VISUALIZATION_LIMIT = 50


def build_path_network_payload(
    path_rows: pd.DataFrame,
    path_edges: pd.DataFrame,
    graph_metadata: pd.DataFrame,
    propagation_matrix: pd.DataFrame,
    *,
    directed_graph: bool,
    path_limit: int = PATH_NETWORK_VISUALIZATION_LIMIT,
) -> dict[str, Any]:
    """Build a compact, auditable union graph from the highest-ranked paths.

    Repeated edges are represented once and retain the complete list of path
    ranks that use them.  A displayed arrow means the propagation matrix truly
    permits only one traversal direction; unresolved edges stay arrowless even
    though each serialized path necessarily traverses them in one direction.
    """
    available_path_count = int(len(path_rows))
    empty = {
        "available_path_count": available_path_count,
        "visualized_path_count": 0,
        "path_limit": int(path_limit),
        "directed_graph": bool(directed_graph),
        "start_symbol": None,
        "target_symbol": None,
        "nodes": [],
        "edges": [],
    }
    if path_rows.empty or path_edges.empty or path_limit < 1:
        return empty

    ranks = (
        pd.to_numeric(path_rows["rank"], errors="coerce")
        .dropna()
        .astype(int)
        .sort_values()
        .drop_duplicates()
        .head(path_limit)
        .tolist()
    )
    if not ranks:
        return empty
    rank_set = set(ranks)
    edges = path_edges.loc[
        pd.to_numeric(path_edges["path_rank"], errors="coerce").isin(rank_set)
    ].copy()
    if edges.empty:
        return empty
    edges["path_rank"] = pd.to_numeric(edges["path_rank"], errors="raise").astype(int)
    edges["step"] = pd.to_numeric(edges["step"], errors="raise").astype(int)
    edges["edge_probability"] = pd.to_numeric(
        edges["edge_probability"], errors="raise"
    ).astype(float)

    metadata = graph_metadata.copy().fillna("")
    metadata["symbol"] = metadata["symbol"].astype(str)
    metadata = metadata.drop_duplicates("symbol", keep="last").set_index("symbol")
    symbols = propagation_matrix.index.astype(str).tolist()
    symbol_position = {symbol: index for index, symbol in enumerate(symbols)}

    path_nodes: dict[str, dict[str, Any]] = {}
    start_symbol: str | None = None
    target_symbol: str | None = None
    for rank, group in edges.groupby("path_rank", sort=True):
        ordered = group.sort_values("step")
        sequence = [str(ordered.iloc[0]["source_symbol"]), *ordered["target_symbol"].astype(str)]
        if start_symbol is None:
            start_symbol = sequence[0]
            target_symbol = sequence[-1]
        denominator = max(len(sequence) - 1, 1)
        for index, symbol in enumerate(sequence):
            record = path_nodes.setdefault(
                symbol,
                {"path_ranks": set(), "path_positions": []},
            )
            record["path_ranks"].add(int(rank))
            record["path_positions"].append(
                {"rank": int(rank), "position": float(index / denominator)}
            )

    node_rows: list[dict[str, Any]] = []
    for symbol, record in path_nodes.items():
        meta = metadata.loc[symbol] if symbol in metadata.index else pd.Series(dtype=object)
        posterior_value = pd.to_numeric(
            pd.Series([meta.get("gui_posterior", "")]), errors="coerce"
        ).iloc[0]
        posterior = None if pd.isna(posterior_value) else float(posterior_value)
        positions = record["path_positions"]
        path_ranks = sorted(int(value) for value in record["path_ranks"])
        node_rows.append(
            {
                "id": symbol,
                "label": str(meta.get("display_symbol", "") or symbol),
                "name": str(meta.get("name", "")),
                "classes": str(meta.get("classes", "")),
                "node_type": str(meta.get("node_type", "")),
                "posterior_probability": posterior,
                "posterior_available": posterior is not None,
                "evidence_basis": (
                    "bayesian_node_posterior"
                    if posterior is not None
                    else "curated_or_external_node"
                ),
                "is_start": symbol == start_symbol,
                "is_target": symbol == target_symbol,
                "best_path_rank": min(path_ranks),
                "path_count": len(path_ranks),
                "path_ranks": path_ranks,
                "path_positions": positions,
                "mean_path_position": float(
                    np.mean([item["position"] for item in positions])
                ),
            }
        )

    pair_records: dict[tuple[str, str], dict[str, Any]] = {}
    for row in edges.itertuples(index=False):
        source = str(row.source_symbol)
        target = str(row.target_symbol)
        if source not in symbol_position or target not in symbol_position:
            continue
        canonical = tuple(sorted((source, target), key=symbol_position.__getitem__))
        record = pair_records.setdefault(
            canonical,
            {
                "probability": float(row.edge_probability),
                "path_ranks": set(),
                "path_traversals": set(),
            },
        )
        record["probability"] = max(
            float(record["probability"]), float(row.edge_probability)
        )
        record["path_ranks"].add(int(row.path_rank))
        record["path_traversals"].add(f"{source}->{target}")

    edge_rows: list[dict[str, Any]] = []
    for (node_a, node_b), record in pair_records.items():
        forward = float(propagation_matrix.loc[node_a, node_b])
        reverse = float(propagation_matrix.loc[node_b, node_a])
        uniquely_directed = bool(directed_graph and ((forward > 0) ^ (reverse > 0)))
        if uniquely_directed and reverse > 0:
            source, target = node_b, node_a
        else:
            source, target = node_a, node_b
        path_ranks = sorted(int(value) for value in record["path_ranks"])
        edge_rows.append(
            {
                "id": f"{node_a}--{node_b}",
                "node_a": node_a,
                "node_b": node_b,
                "source": source,
                "target": target,
                "edge_probability": float(record["probability"]),
                "directionality": (
                    "uniquely_directed"
                    if uniquely_directed
                    else "unresolved_bidirectional"
                    if directed_graph
                    else "undirected"
                ),
                "path_count": len(path_ranks),
                "best_path_rank": min(path_ranks),
                "path_ranks": path_ranks,
                "path_traversals": sorted(record["path_traversals"]),
            }
        )

    node_rows.sort(key=lambda row: (row["mean_path_position"], row["id"]))
    edge_rows.sort(key=lambda row: (row["best_path_rank"], row["id"]))
    return {
        **empty,
        "visualized_path_count": len(ranks),
        "start_symbol": start_symbol,
        "target_symbol": target_symbol,
        "nodes": node_rows,
        "edges": edge_rows,
    }


def run_workflow(
    supplied_configuration: dict[str, Any] | None = None,
    *,
    project_root: Path | str = PROJECT_ROOT,
    run_id: str | None = None,
    progress: ProgressCallback | None = None,
    cancel_requested: CancellationCallback | None = None,
) -> WorkflowResult:
    project = Path(project_root).resolve()
    registry = load_registry(project)
    config = normalize_configuration(supplied_configuration, registry)
    submitted_config = copy.deepcopy(config)
    run_id = run_id or (
        datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8]
    )
    output = project / RUNS_RELATIVE / safe_name(run_id)
    output.mkdir(parents=True, exist_ok=False)

    def check_cancel() -> None:
        if cancel_requested is not None and cancel_requested():
            raise WorkflowCancelled("Analysis cancelled by user")

    def update(message: str, fraction: float) -> None:
        check_cancel()
        if progress:
            progress(message, fraction)

    try:
        calibration_summary: dict[str, Any] | None = None
        node_calibration_targets = pd.DataFrame()
        node_calibration_trace = pd.DataFrame()
        edge_calibration_targets = pd.DataFrame()
        edge_calibration_trace = pd.DataFrame()
        if config["calibration"]["enabled"]:
            calibration_summary = {
                "enabled": True,
                "model": "positive-control target fitting with one quadratic regularizer",
                "unknown_hypotheses_are_negative_labels": False,
                "node": None,
                "edge": None,
            }
            if config["calibration"]["known_nodes"]:
                update("Preparing node calibration", 0.02)
                (
                    calibration_summary["node"],
                    node_calibration_targets,
                    node_calibration_trace,
                ) = calibrate_node_parameters(
                    project,
                    registry,
                    config,
                    progress=(
                        lambda message, fraction: update(
                            message, 0.02 + 0.12 * float(fraction)
                        )
                    ),
                    check_cancel=check_cancel,
                )

        update("Selecting nodes", 0.16 if calibration_summary else 0.08)
        node_factors, selected_nodes, node_summary = select_nodes(
            project, registry, config
        )
        warnings: list[str] = []
        universe = pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
        selected_symbols = selected_nodes["symbol"].astype(str).tolist()

        target_symbol = ""
        target_external = False
        target_resolution: dict[str, Any] | None = None
        graph_symbols = selected_symbols.copy()
        graph_metadata = selected_nodes.copy()
        if config["path"]["enabled"]:
            target_symbol, graph_metadata, target_resolution = append_requested_graph_endpoint(
                config["path"]["target"],
                graph_symbols,
                graph_metadata,
                universe,
                node_factors,
                project,
                reason="path_endpoint_incremental",
            )
            target_external = target_symbol not in set(universe["symbol"].astype(str))

        resolved_known_edges: list[tuple[str, str]] = []
        calibration_endpoint_resolutions: list[dict[str, Any]] = []
        if calibration_summary and config["calibration"]["known_edges"]:
            seen_resolved: set[tuple[str, str]] = set()
            for requested_left, requested_right in config["calibration"]["known_edges"]:
                left, graph_metadata, left_resolution = append_requested_graph_endpoint(
                    requested_left,
                    graph_symbols,
                    graph_metadata,
                    universe,
                    node_factors,
                    project,
                    reason="calibration_edge_endpoint",
                )
                right, graph_metadata, right_resolution = append_requested_graph_endpoint(
                    requested_right,
                    graph_symbols,
                    graph_metadata,
                    universe,
                    node_factors,
                    project,
                    reason="calibration_edge_endpoint",
                )
                if left.casefold() == right.casefold():
                    raise ValueError(
                        f"known edge {requested_left},{requested_right} resolves to a self-edge"
                    )
                canonical = _canonical_undirected_pair(left, right)
                if canonical not in seen_resolved:
                    resolved_known_edges.append((left, right))
                    seen_resolved.add(canonical)
                calibration_endpoint_resolutions.append(
                    {
                        "requested_node_a": requested_left,
                        "requested_node_b": requested_right,
                        "resolved_node_a": left,
                        "resolved_node_b": right,
                        "node_a_resolution": left_resolution,
                        "node_b_resolution": right_resolution,
                    }
                )
            forced_endpoints = sorted(
                {
                    symbol
                    for pair in resolved_known_edges
                    for symbol in pair
                    if symbol not in set(selected_symbols)
                },
                key=lambda value: (value.casefold(), value),
            )
            calibration_summary["edge_endpoints_added_beyond_node_selection"] = (
                forced_endpoints
            )
            if forced_endpoints:
                warnings.append(
                    "Known-edge calibration forced these endpoints into the graph "
                    "even though node selection did not select them: "
                    + ", ".join(forced_endpoints)
                    + "."
                )

        seed_symbols = universe["symbol"].astype(str).tolist()
        cache_summary: dict[str, Any] | None = None
        if any(symbol not in set(seed_symbols) for symbol in graph_symbols):
            update("Characterizing and caching new edge pairs", 0.22)
            cache_update = ensure_incremental_pairs(
                project,
                graph_metadata,
                seed_symbols,
                progress=(
                    (lambda message, fraction: update(message, fraction))
                    if progress
                    else None
                ),
                cancel_check=check_cancel,
            )
            cache_summary = cache_update.as_dict()

        if calibration_summary and resolved_known_edges:
            update("Preparing edge calibration", 0.48)
            (
                calibration_summary["edge"],
                edge_calibration_targets,
                edge_calibration_trace,
            ) = calibrate_edge_parameters(
                project,
                registry,
                config,
                graph_symbols,
                graph_metadata,
                resolved_known_edges,
                progress=(
                    lambda message, fraction: update(
                        message, 0.48 + 0.12 * float(fraction)
                    )
                ),
                check_cancel=check_cancel,
            )
            calibration_summary["edge"]["endpoint_resolution"] = (
                calibration_endpoint_resolutions
            )

        update("Integrating edge evidence", 0.64 if calibration_summary else 0.62)
        derived_audits: dict[str, pd.DataFrame] = {}
        matrix, edge_summary = combine_edge_factors(
            project,
            registry,
            config,
            graph_symbols,
            graph_metadata=graph_metadata,
            audit_collector=derived_audits,
        )
        if "scaffold_triadic_closure" in derived_audits:
            warnings.append(
                "Scaffold-mediated closure is derived from the selected pre-closure "
                "edge graph. It is dependent proximity/co-complex evidence, not "
                "independent proof of a direct binary PPI."
            )
        cutoff = config["edge_integration"]["output_probability_cutoff"]
        supported = supported_edge_table(matrix, cutoff)
        propagation_matrix = matrix
        directionality_audit = pd.DataFrame()
        directionality_class_catalog = pd.DataFrame()
        directionality_summary: dict[str, Any] | None = None
        directionality_catalog: dict[str, Any] | None = None
        omnipath_direction_evidence = pd.DataFrame()
        if (
            config["path"]["enabled"]
            and config["path"]["ontology_directionality_enabled"]
        ):
            update("Applying ontology directionality", 0.69)
            directionality_catalog = load_direction_rule_catalog(
                RULE_CATALOG_PATH,
                known_classes=[
                    item["id"] for item in registry["path_ontology_classes"]
                ],
            )
            directionality_class_catalog = build_complete_class_pair_catalog(
                directionality_catalog,
                [item["id"] for item in registry["path_ontology_classes"]],
            )
            if config["path"]["omnipath_directionality_enabled"]:
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
                    graph_symbols,
                )
            orientation_cutoff = min(
                float(cutoff),
                float(config["path"]["minimum_edge_probability"]),
            )
            (
                propagation_matrix,
                directionality_audit,
                directionality_summary,
            ) = apply_ontology_directionality(
                matrix,
                graph_metadata,
                directionality_catalog,
                audit_probability_cutoff=orientation_cutoff,
                edge_output_cutoff=float(cutoff),
                path_probability_cutoff=float(
                    config["path"]["minimum_edge_probability"]
                ),
                omnipath_directions=(
                    omnipath_direction_evidence
                    if config["path"]["omnipath_directionality_enabled"]
                    else None
                ),
            )
            directionality_summary["omnipath_direction_source"] = (
                str(OMNIPATH_DIRECTION_RAW_RELATIVE)
                if config["path"]["omnipath_directionality_enabled"]
                else None
            )
        path_rows = pd.DataFrame()
        path_edges = pd.DataFrame()
        eligibility = pd.DataFrame()
        path_summary: dict[str, Any] | None = None
        temporal_path_rows = pd.DataFrame()
        temporal_gene_rows = pd.DataFrame()
        temporal_trend_rows = pd.DataFrame()
        temporal_summary: dict[str, Any] | None = None
        if config["path"]["enabled"]:
            update("Ranking paths", 0.74)
            path_rows, path_edges, eligibility, path_summary = run_paths(
                propagation_matrix,
                graph_metadata,
                config,
                config["path"]["start"],
                target_symbol,
                directed_graph=directionality_summary is not None,
            )
            if path_summary is not None:
                path_summary["ontology_directionality"] = directionality_summary

        path_network = build_path_network_payload(
            path_rows,
            path_edges,
            graph_metadata,
            propagation_matrix,
            directed_graph=directionality_summary is not None,
        )
        if path_summary is not None:
            path_summary["merged_path_visualization"] = {
                "available_path_count": path_network["available_path_count"],
                "visualized_path_count": path_network["visualized_path_count"],
                "path_limit": path_network["path_limit"],
                "merged_node_count": len(path_network["nodes"]),
                "merged_edge_count": len(path_network["edges"]),
                "uniquely_directed_edge_count": sum(
                    edge["directionality"] == "uniquely_directed"
                    for edge in path_network["edges"]
                ),
                "artifact": "top_path_network.json",
            }

        if config["temporal_validation"]["enabled"]:
            temporal_config = config["temporal_validation"]
            if path_rows.empty:
                temporal_summary = {
                    "enabled": True,
                    "executed": False,
                    "reason": "No Bayesian paths were available for temporal validation.",
                    "primary_bayesian_rank_changed": False,
                }
                warnings.append(
                    "Temporal validation was enabled, but no supported Bayesian path "
                    "was available to score."
                )
            else:
                temporal_source = project / TEMPORAL_PHOSPHO_RELATIVE
                update("Loading dDAVP temporal phosphoproteomics", 0.76)
                temporal_sites = load_site_trajectories(
                    temporal_source,
                    progress=(
                        lambda message, fraction: update(
                            message, 0.76 + 0.04 * float(fraction)
                        )
                    ),
                    cancel_check=lambda: (check_cancel() or False),
                )
                temporal_responses, temporal_trend = build_replicate_responses(
                    temporal_sites,
                    prior_df=temporal_config["prior_df"],
                    alpha=temporal_config["alpha"],
                    n_boot=temporal_config["monte_carlo_draws"],
                    seed=temporal_config["random_seed"],
                    p_adjust_method=temporal_config["p_adjust_method"],
                    progress=(
                        lambda message, fraction: update(
                            message, 0.80 + 0.08 * float(fraction)
                        )
                    ),
                    cancel_check=lambda: (check_cancel() or False),
                )
                update("Scoring temporal path order", 0.89)
                temporal_path_rows = rank_paths_replicate(
                    path_rows,
                    temporal_responses,
                    minimum_scored_nodes=temporal_config["minimum_scored_nodes"],
                )
                temporal_gene_rows = replicate_response_table(temporal_responses)
                temporal_trend_rows = variance_trend_table(temporal_trend)
                temporal_summary = {
                    **temporal_validation_summary(
                        temporal_path_rows,
                        temporal_responses,
                        site_count=len(temporal_sites),
                        prior_df=temporal_config["prior_df"],
                        alpha=temporal_config["alpha"],
                        n_boot=temporal_config["monte_carlo_draws"],
                        seed=temporal_config["random_seed"],
                        p_adjust_method=temporal_config["p_adjust_method"],
                        minimum_scored_nodes=temporal_config[
                            "minimum_scored_nodes"
                        ],
                    ),
                    "executed": True,
                    "source_file": str(TEMPORAL_PHOSPHO_RELATIVE).replace("\\", "/"),
                }
                if temporal_summary["temporally_informative_path_count"] == 0:
                    warnings.append(
                        "Temporal validation found no path with enough significant, "
                        "measured nodes for the configured temporal evidence rank."
                    )
                elif temporal_summary["temporally_informative_path_count"] < len(path_rows):
                    warnings.append(
                        "Temporal ordering is only interpretable for "
                        f"{temporal_summary['temporally_informative_path_count']} of "
                        f"{len(path_rows)} paths at the configured coverage gate."
                    )

        update("Writing reproducible outputs", 0.92)
        config_path = output / "configuration.json"
        submitted_config_path = output / "submitted_configuration.json"
        calibration_summary_path = output / "calibration_summary.json"
        node_calibration_parameters_path = output / "node_calibrated_parameters.tsv"
        node_calibration_targets_path = output / "node_calibration_targets.tsv"
        node_calibration_trace_path = output / "node_calibration_optimizer_trace.tsv.gz"
        edge_calibration_parameters_path = output / "edge_calibrated_parameters.tsv"
        edge_calibration_targets_path = output / "edge_calibration_targets.tsv"
        edge_calibration_trace_path = output / "edge_calibration_optimizer_trace.tsv.gz"
        node_factors_path = output / "node_posteriors.tsv.gz"
        selected_nodes_path = output / "selected_nodes.tsv"
        matrix_path = output / "edge_adjacency_matrix.tsv"
        supported_path = output / "supported_edges.tsv.gz"
        propagation_matrix_path = output / "propagation_adjacency_matrix.tsv"
        directionality_audit_path = output / "ontology_directionality_audit.tsv.gz"
        directionality_class_catalog_path = output / "ontology_class_pair_catalog.tsv"
        directionality_rules_path = output / "ontology_direction_rules.json"
        omnipath_direction_evidence_path = output / "omnipath_direction_evidence.tsv.gz"
        temporal_paths_path = output / "ranked_paths_temporal.tsv"
        temporal_genes_path = output / "temporal_gene_responses.tsv.gz"
        temporal_trend_path = output / "temporal_variance_trend.tsv"
        temporal_summary_path = output / "temporal_validation_summary.json"
        path_network_path = output / "top_path_network.json"
        summary_path = output / "analysis_summary.json"
        check_cancel()
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        check_cancel()
        node_factors.to_csv(node_factors_path, sep="\t", index=False, compression="gzip")
        check_cancel()
        selected_nodes.to_csv(selected_nodes_path, sep="\t", index=False)
        check_cancel()
        matrix.to_csv(matrix_path, sep="\t", float_format="%.9g")
        check_cancel()
        supported.to_csv(supported_path, sep="\t", index=False, compression="gzip")
        files = [
            config_path,
            node_factors_path,
            selected_nodes_path,
            matrix_path,
            supported_path,
        ]
        if calibration_summary is not None:
            warnings.append(
                "Calibration used known-present positive controls only. It improves "
                "their fitted probabilities but does not by itself demonstrate "
                "specificity, population calibration, or performance on unknown pairs."
            )
            submitted_config_path.write_text(
                json.dumps(submitted_config, indent=2) + "\n", encoding="utf-8"
            )
            calibration_summary_path.write_text(
                json.dumps(calibration_summary, indent=2, default=str) + "\n",
                encoding="utf-8",
            )
            files.extend([submitted_config_path, calibration_summary_path])
            if calibration_summary.get("node") is not None:
                pd.DataFrame(calibration_summary["node"]["parameters"]).to_csv(
                    node_calibration_parameters_path, sep="\t", index=False
                )
                node_calibration_targets.to_csv(
                    node_calibration_targets_path, sep="\t", index=False
                )
                node_calibration_trace.to_csv(
                    node_calibration_trace_path,
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
                files.extend(
                    [
                        node_calibration_parameters_path,
                        node_calibration_targets_path,
                        node_calibration_trace_path,
                    ]
                )
            if calibration_summary.get("edge") is not None:
                pd.DataFrame(calibration_summary["edge"]["parameters"]).to_csv(
                    edge_calibration_parameters_path, sep="\t", index=False
                )
                edge_calibration_targets.to_csv(
                    edge_calibration_targets_path, sep="\t", index=False
                )
                edge_calibration_trace.to_csv(
                    edge_calibration_trace_path,
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
                files.extend(
                    [
                        edge_calibration_parameters_path,
                        edge_calibration_targets_path,
                        edge_calibration_trace_path,
                    ]
                )
        if directionality_summary is not None and directionality_catalog is not None:
            check_cancel()
            propagation_matrix.to_csv(
                propagation_matrix_path,
                sep="\t",
                float_format="%.9g",
            )
            check_cancel()
            directionality_audit.to_csv(
                directionality_audit_path,
                sep="\t",
                index=False,
                compression="gzip",
            )
            directionality_class_catalog.to_csv(
                directionality_class_catalog_path,
                sep="\t",
                index=False,
            )
            serialized_catalog = {
                key: value
                for key, value in directionality_catalog.items()
                if key != "catalog_path"
            }
            directionality_rules_path.write_text(
                json.dumps(serialized_catalog, indent=2) + "\n",
                encoding="utf-8",
            )
            files.extend(
                [
                    propagation_matrix_path,
                    directionality_audit_path,
                    directionality_class_catalog_path,
                    directionality_rules_path,
                ]
            )
            if config["path"]["omnipath_directionality_enabled"]:
                omnipath_direction_evidence.to_csv(
                    omnipath_direction_evidence_path,
                    sep="\t",
                    index=False,
                    compression="gzip",
                )
                files.append(omnipath_direction_evidence_path)
        for stream_id, audit in derived_audits.items():
            check_cancel()
            audit_path = output / f"{safe_name(stream_id)}_audit.tsv.gz"
            audit.to_csv(audit_path, sep="\t", index=False, compression="gzip")
            files.append(audit_path)
        if config["path"]["enabled"]:
            check_cancel()
            paths_path = output / "ranked_paths.tsv"
            path_edges_path = output / "ranked_path_edges.tsv"
            eligibility_path = output / "intermediate_node_eligibility.tsv"
            path_rows.to_csv(paths_path, sep="\t", index=False, float_format="%.12g")
            path_edges.to_csv(path_edges_path, sep="\t", index=False, float_format="%.12g")
            eligibility.to_csv(eligibility_path, sep="\t", index=False)
            files.extend([paths_path, path_edges_path, eligibility_path])
            path_network_path.write_text(
                json.dumps(path_network, indent=2) + "\n", encoding="utf-8"
            )
            files.append(path_network_path)
        if temporal_summary is not None:
            check_cancel()
            temporal_summary_path.write_text(
                json.dumps(temporal_summary, indent=2, default=str) + "\n",
                encoding="utf-8",
            )
            files.append(temporal_summary_path)
            if temporal_summary.get("executed"):
                temporal_path_rows.to_csv(
                    temporal_paths_path,
                    sep="\t",
                    index=False,
                    float_format="%.12g",
                )
                temporal_gene_rows.to_csv(
                    temporal_genes_path,
                    sep="\t",
                    index=False,
                    compression="gzip",
                    float_format="%.12g",
                )
                temporal_trend_rows.to_csv(
                    temporal_trend_path,
                    sep="\t",
                    index=False,
                    float_format="%.12g",
                )
                files.extend(
                    [temporal_paths_path, temporal_genes_path, temporal_trend_path]
                )

        summary = {
            "schema_version": 1,
            "run_id": run_id,
            "generated_at": utc_now(),
            "calibration": calibration_summary,
            "node_selection": node_summary,
            "edge_characterization": {
                **edge_summary,
                "incremental_pair_cache": cache_summary,
                "matrix_node_count_after_target_extension": int(len(matrix)),
                "reported_supported_edge_count": int(len(supported)),
            },
            "path_finding": path_summary,
            "temporal_validation": temporal_summary,
            "ontology_directionality": directionality_summary,
            "external_target": (
                {
                    "symbol": target_symbol,
                    "resolution": target_resolution,
                    "source_summary": cache_summary,
                }
                if target_external
                else None
            ),
            "warnings": warnings,
            "outputs": [path.name for path in files],
        }
        summary_path.write_text(json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")
        preview = {
            "metrics": {
                "selected_nodes": node_summary["selected_node_count"],
                "supported_edges": int(len(supported)),
                "ranked_paths": int(len(path_rows)),
            },
            "top_nodes": _records(
                node_factors.sort_values("gui_posterior", ascending=False),
                ["gene_symbol", "gui_posterior", "gui_rank"],
                10,
            ),
            "top_paths": _records(
                (
                    temporal_path_rows
                    if not temporal_path_rows.empty
                    else path_rows
                ),
                [
                    "rank",
                    "hop_count",
                    "path_probability_product",
                    "path_symbols",
                    *(
                        [
                            "temporal_n_scored",
                            "temporal_kendall_tau_mean",
                            "temporal_kendall_tau_low",
                            "temporal_kendall_tau_high",
                            "temporal_order_informative",
                            "temporal_evidence_rank",
                        ]
                        if not temporal_path_rows.empty
                        else []
                    ),
                ],
                config["path"]["top_k"],
            )
            if not path_rows.empty
            else [],
            "path_network": path_network,
            "probability_distributions": {
                "nodes": node_summary["probability_distribution"],
                "edges": edge_summary["probability_distribution"],
            },
            "directionality": directionality_summary,
            "temporal_validation": temporal_summary,
            "calibration": calibration_summary,
            "warnings": warnings,
            "files": [path.name for path in [*files, summary_path]],
        }
        update("Complete", 1.0)
        return WorkflowResult(run_id, output, summary, preview)
    except Exception:
        # Keep the run directory and its configuration when possible so a failed
        # scientific run can be audited without overwriting previous results.
        failure_config = output / "configuration.json"
        if not failure_config.exists():
            failure_config.write_text(
                json.dumps(config, indent=2) + "\n", encoding="utf-8"
            )
        raise


if __name__ == "__main__":
    result = run_workflow()
    print(json.dumps(result.preview, indent=2))
    print(f"Output directory: {result.output_directory}")
