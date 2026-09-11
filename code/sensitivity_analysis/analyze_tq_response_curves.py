#!/usr/bin/env python3
"""Trace Bayes-factor and information response curves over every stream's Tq range.

This is a companion to ``analyze_evidence_information.py``.  It deliberately
uses the current positive-or-neutral scoring mode, the registry's default
weights, and one stream at a time.  The x-axis is a scale multiplier because
several streams store a distinct empirical Tq for every site or protein.

STRING and STITCH do not use the Gaussian-complement Tq.  Their same GUI
control scales the reference probability used in the odds-ratio Bayes factor;
their rows are therefore labelled ``Ref x`` rather than ``Tq x``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
GUI_DIR = PROJECT_ROOT / "gui"
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))
if str(SCRIPT_PATH.parent) not in sys.path:
    sys.path.insert(0, str(SCRIPT_PATH.parent))

from workflow_engine import (  # noqa: E402
    UNIVERSE_RELATIVE,
    _calibration_tq_bounds,
    default_configuration,
    load_node_factor_catalog,
    load_registry,
)

from analyze_evidence_information import (  # noqa: E402
    edge_upper_bundle,
    evidence_metrics,
    node_stream_bundle,
)


DEFAULT_GRID_POINTS = 41
MODE = "current_positive_or_neutral"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Project root containing gui/evidence_registry.json.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory (default: results/sensitivity_analysis/"
            "evidence_information_2026-08-24/tq_response_curves)."
        ),
    )
    parser.add_argument(
        "--grid-points",
        type=int,
        default=DEFAULT_GRID_POINTS,
        help="Number of log-spaced multiplier values per stream (minimum 9).",
    )
    return parser.parse_args()


def multiplier_grid(definition: dict[str, Any], point_count: int) -> np.ndarray:
    """Return a log-spaced GUI-calibration grid with multiplier 1 included."""
    if point_count < 9:
        raise ValueError("--grid-points must be at least 9")
    lower, upper = _calibration_tq_bounds(definition)
    values = np.geomspace(lower, upper, point_count)
    if lower <= 1.0 <= upper:
        values = np.append(values, 1.0)
    return np.unique(np.round(values, 12))


def response_row(
    definition: dict[str, Any],
    state: dict[str, Any],
    multiplier: float,
    bayes_factors: np.ndarray,
    eligible: np.ndarray,
    *,
    stage: str,
) -> dict[str, Any]:
    weight = float(state["weight"])
    metrics = evidence_metrics(
        definition["id"],
        definition["label"],
        stage,
        bayes_factors,
        eligible,
        enabled_by_default=bool(state["enabled"]),
        weight=weight,
        mode=MODE,
    )
    mean_signed_bits = float(metrics["mean_signed_bits_per_eligible"])
    geometric_mean_bf = (
        float(2.0**mean_signed_bits) if math.isfinite(mean_signed_bits) else math.nan
    )
    lower, upper = _calibration_tq_bounds(definition)
    return {
        "stage": stage,
        "stream_id": definition["id"],
        "stream_label": definition["label"],
        "control_label": definition["normalization"]["control_label"],
        "normalization_handler": definition["normalization"]["handler"],
        "multiplier": float(multiplier),
        "calibration_lower_bound": float(lower),
        "calibration_upper_bound": float(upper),
        "enabled_by_default": bool(state["enabled"]),
        "weight": weight,
        "eligible_count": int(metrics["eligible_count"]),
        "nonneutral_count": int(metrics["nonneutral_count"]),
        # This is the geometric mean of BF**weight over eligible hypotheses.
        # It is stable where an arithmetic mean or product would be dominated by
        # a few very large factors or overflow numerically.
        "geometric_mean_effective_bayes_factor": geometric_mean_bf,
        "total_signed_log2_bf_bits": float(metrics["signed_net_bits"]),
        "mean_signed_log2_bf_bits_per_eligible": mean_signed_bits,
        "total_absolute_log2_bf_bits": float(metrics["total_absolute_bits"]),
        "mean_absolute_log2_bf_bits_per_eligible": float(
            metrics["mean_absolute_bits_per_eligible"]
        ),
        "total_information_kl_bits": float(
            metrics["total_posterior_kl_bits_from_prior_0_5"]
        ),
        "mean_information_kl_bits_per_eligible": float(
            metrics["mean_posterior_kl_bits_per_eligible"]
        ),
    }


def maxima_table(curves: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (_, stream_id), group in curves.groupby(["stage", "stream_id"], sort=False):
        group = group.sort_values("multiplier")
        bf_index = group["geometric_mean_effective_bayes_factor"].idxmax()
        info_index = group["total_information_kl_bits"].idxmax()
        bf_row = group.loc[bf_index]
        info_row = group.loc[info_index]
        rows.append(
            {
                "stage": bf_row["stage"],
                "stream_id": stream_id,
                "stream_label": bf_row["stream_label"],
                "control_label": bf_row["control_label"],
                "enabled_by_default": bool(bf_row["enabled_by_default"]),
                "weight": float(bf_row["weight"]),
                "calibration_lower_bound": float(bf_row["calibration_lower_bound"]),
                "calibration_upper_bound": float(bf_row["calibration_upper_bound"]),
                "bf_maximizing_multiplier_on_grid": float(bf_row["multiplier"]),
                "maximum_geometric_mean_effective_bf": float(
                    bf_row["geometric_mean_effective_bayes_factor"]
                ),
                "information_maximizing_multiplier_on_grid": float(
                    info_row["multiplier"]
                ),
                "maximum_total_information_kl_bits": float(
                    info_row["total_information_kl_bits"]
                ),
                "bf_maximum_at_search_boundary": bool(
                    math.isclose(
                        float(bf_row["multiplier"]),
                        float(bf_row["calibration_lower_bound"]),
                        rel_tol=1e-9,
                        abs_tol=1e-12,
                    )
                    or math.isclose(
                        float(bf_row["multiplier"]),
                        float(bf_row["calibration_upper_bound"]),
                        rel_tol=1e-9,
                        abs_tol=1e-12,
                    )
                ),
                "information_maximum_at_search_boundary": bool(
                    math.isclose(
                        float(info_row["multiplier"]),
                        float(info_row["calibration_lower_bound"]),
                        rel_tol=1e-9,
                        abs_tol=1e-12,
                    )
                    or math.isclose(
                        float(info_row["multiplier"]),
                        float(info_row["calibration_upper_bound"]),
                        rel_tol=1e-9,
                        abs_tol=1e-12,
                    )
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = parse_args()
    project = args.project_root.resolve()
    output = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else project
        / "results/sensitivity_analysis/evidence_information_2026-08-24/tq_response_curves"
    )
    output.mkdir(parents=True, exist_ok=True)

    registry = load_registry(project)
    config = default_configuration(registry)
    node_definitions = {item["id"]: item for item in registry["node_streams"]}
    edge_definitions = {
        item["id"]: item for item in registry["edge_streams"] if not item.get("derived")
    }
    node_factors = load_node_factor_catalog(project, registry)
    symbols = (
        pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str)["symbol"]
        .astype(str)
        .tolist()
    )
    upper_i, upper_j = np.triu_indices(len(symbols), 1)

    rows: list[dict[str, Any]] = []
    for stream_id, definition in node_definitions.items():
        state = config["node_streams"][stream_id]
        grid = multiplier_grid(definition, args.grid_points)
        for multiplier in grid:
            varied_state = dict(state)
            varied_state["tq_multiplier"] = float(multiplier)
            bundle = node_stream_bundle(
                project, node_factors, definition, varied_state, continuous=False
            )
            rows.append(
                response_row(
                    definition,
                    state,
                    multiplier,
                    bundle["bayes_factors"],
                    bundle["eligible"],
                    stage="node",
                )
            )
        print(f"Traced node stream: {definition['label']} ({len(grid)} points)", flush=True)

    for stream_id, definition in edge_definitions.items():
        state = config["edge_streams"][stream_id]
        grid = multiplier_grid(definition, args.grid_points)
        for multiplier in grid:
            varied_state = dict(state)
            varied_state["tq_multiplier"] = float(multiplier)
            bundle = edge_upper_bundle(
                project,
                definition,
                varied_state,
                symbols,
                upper_i,
                upper_j,
                continuous=False,
            )
            rows.append(
                response_row(
                    definition,
                    state,
                    multiplier,
                    bundle["bayes_factors"],
                    bundle["eligible"],
                    stage="edge",
                )
            )
        print(f"Traced edge stream: {definition['label']} ({len(grid)} points)", flush=True)

    curves = pd.DataFrame(rows)
    maxima = maxima_table(curves)
    node_curves = curves.loc[curves["stage"].eq("node")].reset_index(drop=True)
    edge_curves = curves.loc[curves["stage"].eq("edge")].reset_index(drop=True)
    node_curves.to_csv(output / "node_tq_response_curves.tsv", sep="\t", index=False)
    edge_curves.to_csv(
        output / "edge_tq_response_curves_seed891.tsv", sep="\t", index=False
    )
    maxima.to_csv(output / "tq_response_curve_maxima.tsv", sep="\t", index=False)

    summary = {
        "mode": MODE,
        "weight_source": "default_configuration(registry)",
        "x_axis": "stream-specific Tq or reference-score multiplier",
        "grid_points_requested": int(args.grid_points),
        "node_stream_count": int(node_curves["stream_id"].nunique()),
        "edge_stream_count": int(edge_curves["stream_id"].nunique()),
        "bf_metric": "geometric mean of BF**weight over eligible hypotheses",
        "information_metric": (
            "sum over eligible hypotheses of D_KL(Bernoulli(BF**weight/(1+BF**weight)) "
            "|| Bernoulli(0.5)), in bits"
        ),
        "bf_boundary_maximum_stream_count": int(
            maxima["bf_maximum_at_search_boundary"].sum()
        ),
        "information_boundary_maximum_stream_count": int(
            maxima["information_maximum_at_search_boundary"].sum()
        ),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    methods = """# Tq-versus-evidence response curves

The response curves vary one stream's scale multiplier at a time over the same
bounds permitted by calibration in the GUI. Values are log-spaced and multiplier
1 is included exactly. All other settings are irrelevant to this standalone
calculation; each stream is evaluated on its own, with its registry-default
weight. The current positive-or-neutral Bayes-factor mapping is used.

The Bayes-factor curve reports the geometric mean of the effective Bayes factor
`BF^weight` over hypotheses eligible for that stream. A geometric mean is used
because Bayes factors combine multiplicatively; a product would overflow and an
arithmetic mean would be dominated by a few very large factors.

The information curve reports total realized KL information from an equal prior:
for each eligible hypothesis, `p = BF^weight / (1 + BF^weight)` and information is
`D_KL(Bernoulli(p) || Bernoulli(0.5))`, summed in bits.

The x-axis is a multiplier, not one universal raw Tq. Several streams have a
different empirical Tq for each protein or site. STRING and STITCH instead scale
their reference probability, so their x-axis is labelled `Ref x`.

Important: this is a sensitivity curve, not a valid unsupervised Tq optimizer.
Under the positive-or-neutral Gaussian-complement mapping, reducing Tq cannot
weaken a Bayes factor. Maximizing either Bayes-factor strength or realized KL
information alone therefore tends to select the smallest allowed multiplier.
An interior Tq must be chosen from an external criterion, such as predictive
performance on held-out known nodes/edges plus regularization.
"""
    (output / "METHODS.md").write_text(methods, encoding="utf-8")
    print(f"Wrote Tq response curves to {output}", flush=True)


if __name__ == "__main__":
    main()
