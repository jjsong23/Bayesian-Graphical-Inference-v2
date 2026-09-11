"""Fast structural tests for the all-stream end-to-end ablation design."""

from __future__ import annotations

import unittest

import analyze_end_to_end_ablation as ablation
import numpy as np
import pandas as pd

from workflow_engine import (
    UNIVERSE_RELATIVE,
    combine_edge_factors,
    default_configuration,
    load_registry,
    normalize_configuration,
)


class AllStreamAblationConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = load_registry(ablation.PROJECT_ROOT)
        cls.supplied = normalize_configuration(
            default_configuration(cls.registry),
            cls.registry,
        )

    def test_optional_node_gets_a_matched_add_one_context(self) -> None:
        context_config, context = ablation.focal_context_configuration(
            self.supplied,
            self.registry,
            stage="node",
            stream_id="collecting_duct_rna_imcd",
        )
        self.assertTrue(
            context_config["node_streams"]["collecting_duct_rna_imcd"]["enabled"]
        )
        self.assertFalse(
            context_config["node_streams"]["collecting_duct_rna_ccd"]["enabled"]
        )
        self.assertEqual(context, "supplied_plus[node:collecting_duct_rna_imcd]")

    def test_exclusive_alternative_gets_a_matched_full_model(self) -> None:
        high, context = ablation.focal_context_configuration(
            self.supplied,
            self.registry,
            stage="edge",
            stream_id="hpa_high_confidence",
        )
        self.assertTrue(high["edge_streams"]["hpa_high_confidence"]["enabled"])
        self.assertFalse(high["edge_streams"]["hpa_primary"]["enabled"])
        self.assertEqual(context, "supplied_plus[edge:hpa_high_confidence]")

    def test_primary_context_is_shared_by_nonexclusive_streams(self) -> None:
        _, node_context = ablation.focal_context_configuration(
            self.supplied,
            self.registry,
            stage="node",
            stream_id="protein_abundance",
        )
        _, primary_context = ablation.focal_context_configuration(
            self.supplied,
            self.registry,
            stage="edge",
            stream_id="hpa_primary",
        )
        self.assertEqual(node_context, "supplied_configuration")
        self.assertEqual(primary_context, "supplied_configuration")

    def test_cached_edge_integrator_matches_workflow_with_scaffold_closure(self) -> None:
        config, _ = ablation.focal_context_configuration(
            self.supplied,
            self.registry,
            stage="edge",
            stream_id="scaffold_triadic_closure",
        )
        universe = pd.read_csv(
            ablation.PROJECT_ROOT / UNIVERSE_RELATIVE,
            sep="\t",
            dtype=str,
        ).fillna("")
        requested = ["Prkaca", "Prkacb", "Creb1", "Hsp90aa1", "Calm1", "Pde4d"]
        metadata = universe.loc[universe["symbol"].isin(requested)].copy()
        symbols = metadata["symbol"].astype(str).tolist()
        self.assertEqual(set(symbols), set(requested))
        expected, _ = combine_edge_factors(
            ablation.PROJECT_ROOT,
            self.registry,
            config,
            symbols,
            graph_metadata=metadata,
        )
        observed = ablation.Evaluator(
            ablation.PROJECT_ROOT,
            self.registry,
        )._edge_matrix(config, symbols, metadata)
        self.assertTrue(
            np.allclose(expected.to_numpy(float), observed.to_numpy(float), atol=1e-12)
        )

    def test_optional_stream_grid_spans_low_reference_and_high_values(self) -> None:
        definition = next(
            item
            for item in self.registry["node_streams"]
            if item["id"] == "collecting_duct_rna_ccd"
        )
        self.assertEqual(
            ablation.tq_grid(definition, 1.0, compact=True),
            [0.25, 1.0, 4.0],
        )


if __name__ == "__main__":
    unittest.main()
