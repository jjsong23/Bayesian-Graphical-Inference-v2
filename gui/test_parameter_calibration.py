"""Regression tests for bounded positive-control calibration."""

from __future__ import annotations

import math
import importlib.util
import unittest

import numpy as np

from parameter_calibration import bounded_powell_positive_calibration


class ParameterCalibrationTests(unittest.TestCase):
    @unittest.skipUnless(
        importlib.util.find_spec("scipy") is not None,
        "SciPy is an optional runtime dependency required for calibration",
    )
    def test_regularized_powell_fit_is_deterministic_and_improves_controls(self) -> None:
        specs = [
            {
                "key": "example:weight",
                "stream_id": "example",
                "parameter": "weight",
                "current": 0.2,
                "preferred": 1.0,
                "bounds": (0.0, 3.0),
                "scale": 0.5,
                "transform": "linear",
            },
            {
                "key": "example:tq_multiplier",
                "stream_id": "example",
                "parameter": "tq_multiplier",
                "current": 1.0,
                "preferred": 1.0,
                "bounds": (0.25, 4.0),
                "scale": math.log(2.0),
                "transform": "log",
            },
        ]

        def probabilities(settings: dict[str, float]) -> np.ndarray:
            signal = (
                settings["example:weight"]
                * math.log(8.0)
                / settings["example:tq_multiplier"]
            )
            probability = 1.0 / (1.0 + math.exp(-signal))
            return np.asarray([probability, probability], dtype=float)

        first, first_trace = bounded_powell_positive_calibration(
            specs,
            probabilities,
            target_count=2,
            regularization_strength=0.1,
            multistart_count=2,
        )
        second, _ = bounded_powell_positive_calibration(
            specs,
            probabilities,
            target_count=2,
            regularization_strength=0.1,
            multistart_count=2,
        )
        self.assertGreater(
            first["final_mean_target_probability"],
            first["initial_mean_target_probability"],
        )
        self.assertEqual(first["engine"], "scipy.optimize.minimize")
        self.assertEqual(first["method"], "Powell")
        self.assertTrue(first["deterministic"])
        self.assertFalse(first_trace.empty)
        self.assertAlmostEqual(
            first["final_total_loss"], second["final_total_loss"], places=12
        )
        fitted = {row["parameter"]: row["fitted"] for row in first["parameters"]}
        self.assertGreaterEqual(fitted["weight"], 0.0)
        self.assertLessEqual(fitted["weight"], 3.0)
        self.assertGreaterEqual(fitted["tq_multiplier"], 0.25)
        self.assertLessEqual(fitted["tq_multiplier"], 4.0)


if __name__ == "__main__":
    unittest.main()
