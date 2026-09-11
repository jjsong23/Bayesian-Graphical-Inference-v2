#!/usr/bin/env python3
"""Small, deterministic optimizer for positive-control Bayesian calibration.

The calibration problem is intentionally narrow: supplied nodes or edges are
treated as known-present positive controls.  Powell searches bounded evidence
weights and Tq/reference multipliers for settings that make those controls
probable, while one quadratic penalty keeps the settings near scientifically
chosen anchors.  Unknown hypotheses are not labeled absent and do not enter the
fit loss.
"""

from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np
import pandas as pd


ProbabilityFunction = Callable[[dict[str, float]], np.ndarray]
ProgressFunction = Callable[[str, float], None]
CancellationFunction = Callable[[], None]


def _to_internal(spec: dict[str, Any], value: float) -> float:
    return math.log(value) if spec.get("transform") == "log" else value


def _from_internal(spec: dict[str, Any], value: float) -> float:
    return math.exp(value) if spec.get("transform") == "log" else value


def _decode(
    specs: list[dict[str, Any]], values: np.ndarray
) -> dict[str, float]:
    return {
        spec["key"]: _from_internal(spec, float(value))
        for spec, value in zip(specs, values, strict=True)
    }


def _start_vectors(
    specs: list[dict[str, Any]], count: int
) -> list[np.ndarray]:
    current = np.asarray(
        [_to_internal(spec, float(spec["current"])) for spec in specs],
        dtype=float,
    )
    preferred = np.asarray(
        [_to_internal(spec, float(spec["preferred"])) for spec in specs],
        dtype=float,
    )
    lower = np.asarray(
        [_to_internal(spec, float(spec["bounds"][0])) for spec in specs],
        dtype=float,
    )
    upper = np.asarray(
        [_to_internal(spec, float(spec["bounds"][1])) for spec in specs],
        dtype=float,
    )
    current = np.clip(current, lower, upper)
    preferred = np.clip(preferred, lower, upper)
    candidates = [current, preferred, (current + preferred) / 2.0]
    starts: list[np.ndarray] = []
    for candidate in candidates:
        if not any(np.allclose(candidate, prior, rtol=0.0, atol=1e-12) for prior in starts):
            starts.append(candidate)
        if len(starts) >= count:
            break
    return starts


def bounded_powell_positive_calibration(
    specs: list[dict[str, Any]],
    probability_function: ProbabilityFunction,
    *,
    target_count: int,
    regularization_strength: float,
    multistart_count: int = 2,
    progress: ProgressFunction | None = None,
    check_cancel: CancellationFunction | None = None,
    stage_label: str = "parameters",
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Fit positive controls with bounded, deterministic multi-start Powell.

    ``specs`` use physical values for ``current``, ``preferred``, and ``bounds``.
    A spec with ``transform='log'`` is optimized on the log scale.  Its ``scale``
    is consequently also interpreted on the log scale (normally ``log(2)``).
    """
    if target_count < 1:
        raise ValueError("calibration requires at least one positive-control target")
    if not specs:
        raise ValueError("calibration requires at least one enabled primary evidence stream")
    try:
        from scipy.optimize import minimize
    except ImportError as exc:  # pragma: no cover - exercised only in incomplete installs
        raise RuntimeError(
            "Calibration requires SciPy. Install the project requirements and restart "
            "the workbench."
        ) from exc

    regularization_strength = float(regularization_strength)
    preferred = np.asarray(
        [_to_internal(spec, float(spec["preferred"])) for spec in specs],
        dtype=float,
    )
    scales = np.asarray([float(spec["scale"]) for spec in specs], dtype=float)
    if (~np.isfinite(scales) | (scales <= 0)).any():
        raise ValueError("calibration parameter scales must be finite and positive")
    bounds = [
        (
            _to_internal(spec, float(spec["bounds"][0])),
            _to_internal(spec, float(spec["bounds"][1])),
        )
        for spec in specs
    ]
    starts = _start_vectors(specs, multistart_count)
    trace_rows: list[dict[str, Any]] = []
    evaluation = 0
    best_payload: dict[str, Any] | None = None

    def evaluate(values: np.ndarray, start_number: int) -> tuple[float, np.ndarray, float, float]:
        nonlocal evaluation, best_payload
        if check_cancel:
            check_cancel()
        decoded = _decode(specs, values)
        probabilities = np.asarray(probability_function(decoded), dtype=float)
        if probabilities.shape != (target_count,):
            raise ValueError("calibration probability function returned the wrong shape")
        if (~np.isfinite(probabilities)).any():
            raise ValueError("calibration produced non-finite target probabilities")
        fit_loss = -float(np.mean(np.log(np.clip(probabilities, 1e-12, 1.0))))
        standardized = (np.asarray(values, dtype=float) - preferred) / scales
        regularization_loss = regularization_strength * float(
            np.mean(np.square(standardized))
        )
        total_loss = fit_loss + regularization_loss
        evaluation += 1
        row = {
            "start_number": start_number,
            "evaluation": evaluation,
            "fit_loss": fit_loss,
            "regularization_loss": regularization_loss,
            "total_loss": total_loss,
            "mean_target_probability": float(probabilities.mean()),
            "minimum_target_probability": float(probabilities.min()),
        }
        row.update(decoded)
        trace_rows.append(row)
        if best_payload is None or total_loss < best_payload["total_loss"]:
            best_payload = {
                "internal": np.asarray(values, dtype=float).copy(),
                "settings": decoded,
                "probabilities": probabilities.copy(),
                "fit_loss": fit_loss,
                "regularization_loss": regularization_loss,
                "total_loss": total_loss,
            }
        if progress and (evaluation == 1 or evaluation % 20 == 0):
            # Powell has no exact up-front evaluation count.  This asymptotic
            # display remains responsive without falsely promising a deadline.
            fraction = min(0.94, evaluation / (evaluation + 180.0))
            progress(
                f"Calibrating {stage_label} ({evaluation:,} objective evaluations)",
                fraction,
            )
        return total_loss, probabilities, fit_loss, regularization_loss

    current_internal = np.asarray(
        [_to_internal(spec, float(spec["current"])) for spec in specs],
        dtype=float,
    )
    initial_settings = _decode(specs, current_internal)
    initial_probabilities = np.asarray(
        probability_function(initial_settings), dtype=float
    )
    initial_fit_loss = -float(
        np.mean(np.log(np.clip(initial_probabilities, 1e-12, 1.0)))
    )
    optimizer_runs: list[dict[str, Any]] = []
    for start_number, start in enumerate(starts, start=1):
        if check_cancel:
            check_cancel()

        def objective(values: np.ndarray) -> float:
            return evaluate(values, start_number)[0]

        result = minimize(
            objective,
            start,
            method="Powell",
            bounds=bounds,
            options={
                "xtol": 1e-4,
                "ftol": 1e-7,
                "maxiter": 160,
                "maxfev": 5000,
                "disp": False,
            },
        )
        optimizer_runs.append(
            {
                "start_number": start_number,
                "success": bool(result.success),
                "message": str(result.message),
                "iterations": int(result.nit),
                "function_evaluations": int(result.nfev),
                "reported_total_loss": float(result.fun),
            }
        )

    assert best_payload is not None
    if progress:
        progress(f"Calibrated {stage_label}", 1.0)
    parameters = []
    for spec in specs:
        parameters.append(
            {
                "key": spec["key"],
                "stream_id": spec["stream_id"],
                "parameter": spec["parameter"],
                "current": float(spec["current"]),
                "preferred": float(spec["preferred"]),
                "fitted": float(best_payload["settings"][spec["key"]]),
                "lower_bound": float(spec["bounds"][0]),
                "upper_bound": float(spec["bounds"][1]),
                "regularization_scale": float(spec["scale"]),
                "optimization_scale": spec.get("transform", "linear"),
            }
        )
    summary = {
        "engine": "scipy.optimize.minimize",
        "method": "Powell",
        "deterministic": True,
        "target_label_model": "positive controls only; unknown hypotheses are unlabeled",
        "target_count": int(target_count),
        "regularization_strength": regularization_strength,
        "regularization_definition": (
            "lambda times the mean squared standardized displacement from preferred "
            "weights and log Tq/reference multipliers"
        ),
        "multistart_count_requested": int(multistart_count),
        "multistart_count_executed": int(len(starts)),
        "initial_fit_loss": initial_fit_loss,
        "initial_mean_target_probability": float(initial_probabilities.mean()),
        "initial_minimum_target_probability": float(initial_probabilities.min()),
        "final_fit_loss": float(best_payload["fit_loss"]),
        "final_regularization_loss": float(best_payload["regularization_loss"]),
        "final_total_loss": float(best_payload["total_loss"]),
        "final_mean_target_probability": float(best_payload["probabilities"].mean()),
        "final_minimum_target_probability": float(best_payload["probabilities"].min()),
        "parameters": parameters,
        "optimizer_runs": optimizer_runs,
        "objective_evaluations": int(evaluation),
        "initial_probabilities": initial_probabilities.tolist(),
        "final_probabilities": best_payload["probabilities"].tolist(),
    }
    return summary, pd.DataFrame(trace_rows)
