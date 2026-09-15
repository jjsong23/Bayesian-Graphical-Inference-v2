from __future__ import annotations

import math
import shutil
import unittest
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

from .pipeline import _export_applied_edge_streams, normalize_backward_options


class CompletePipelineTests(unittest.TestCase):
    def test_exported_edge_factor_is_unweighted_and_weight_is_retained(self) -> None:
        output = Path.cwd() / ".test_work" / uuid.uuid4().hex
        output.mkdir(parents=True)
        try:
            contribution = np.zeros((3, 3), dtype=float)
            contribution[0, 1] = contribution[1, 0] = 2.0 * math.log(3.0)
            registry = {
                "edge_streams": [
                    {"id": "example", "label": "Example evidence"},
                    {"id": "disabled", "label": "Disabled evidence"},
                ]
            }
            configuration = {
                "edge_streams": {
                    "example": {
                        "enabled": True,
                        "weight": 2.0,
                        "tq_multiplier": 0.75,
                        "continuous_negative_evidence": True,
                    },
                    "disabled": {"enabled": False, "weight": 1.0},
                }
            }
            streams, manifest = _export_applied_edge_streams(
                output,
                ["A", "B", "C"],
                registry,
                configuration,
                {"example": contribution},
            )
            table = pd.read_csv(streams[0]["file"], sep="\t")
            self.assertEqual(streams[0]["weight"], 2.0)
            self.assertEqual(streams[0]["missing_bayes_factor"], 1.0)
            self.assertEqual(table[["node_a", "node_b"]].values.tolist(), [["A", "B"]])
            self.assertAlmostEqual(float(table.loc[0, "bayes_factor"]), 3.0)
            self.assertEqual(manifest[0]["non_neutral_pair_count"], 1)
        finally:
            shutil.rmtree(output, ignore_errors=True)

    def test_backward_options_validate_and_merge_structural_defaults(self) -> None:
        result = normalize_backward_options({
            "target": "Aqp2",
            "receptor": "Avpr2",
            "beam_width": 7,
            "structural": {"reference_score": 0.55, "maximum_concurrent_jobs": 4},
        })
        self.assertEqual(result["beam_width"], 7)
        self.assertEqual(result["structural"]["reference_score"], 0.55)
        self.assertEqual(result["structural"]["maximum_concurrent_jobs"], 4)
        self.assertEqual(result["structural"]["bayes_factor_floor"], 0.1)
        with self.assertRaisesRegex(ValueError, "must be different"):
            normalize_backward_options({"target": "Aqp2", "receptor": "Aqp2"})


if __name__ == "__main__":
    unittest.main()
