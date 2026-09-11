"""Regression tests for per-hypothesis GUI evidence inspection."""

from __future__ import annotations

import math
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import evidence_inspector
from evidence_inspector import inspect_edge_evidence, inspect_node_evidence
from workflow_engine import PROJECT_ROOT, UNIVERSE_RELATIVE, default_configuration, load_registry


class EvidenceInspectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = load_registry(PROJECT_ROOT)

    def test_node_ledger_reconciles_exact_recorded_contributions(self) -> None:
        config = default_configuration(self.registry)
        enabled = [
            stream_id
            for stream_id, state in config["node_streams"].items()
            if state["enabled"]
        ]
        assigned = dict(zip(enabled, [2.0, 0.5, 1.0, 1.0]))
        row: dict[str, object] = {
            "gene_symbol": "TestNode",
            "node_name": "Test node",
            "node_classes": "kinase",
            "gui_posterior": 0.5,
            "gui_rank": 1,
            "gui_selected_in_graph": False,
            "kinase_absolute_lfc": 0.2,
            "maximum_absolute_site_lfc": 0.3,
        }
        for stream_id in enabled:
            factor = assigned[stream_id]
            prefix = f"gui_{stream_id}_"
            row.update(
                {
                    prefix + "factor": factor,
                    prefix + "source_bayes_factor": factor,
                    prefix + "bayes_factor": factor,
                    prefix + "weighted_log_bayes_factor": math.log(factor),
                    prefix + "negative_evidence_eligible": True,
                    prefix + "observed": True,
                    prefix + "unobserved_penalty_applied": False,
                    prefix + "continuous_negative_evidence_applied": False,
                }
            )
        with (
            patch.object(Path, "is_file", return_value=True),
            patch.object(evidence_inspector, "_read_configuration", return_value=config),
            patch.object(evidence_inspector, "load_registry", return_value=self.registry),
            patch.object(evidence_inspector.pd, "read_csv", return_value=pd.DataFrame([row])),
        ):
            payload = inspect_node_evidence(Path("mock-run"), "testnode")
        self.assertEqual(payload["symbol"], "TestNode")
        self.assertAlmostEqual(payload["reconstructed_posterior_probability"], 0.5)
        self.assertLess(payload["reconciliation_absolute_difference"], 1e-12)
        statuses = {item["stream_id"]: item["status"] for item in payload["streams"]}
        self.assertEqual(statuses[enabled[0]], "supports")
        self.assertEqual(statuses[enabled[1]], "refutes")
        enabled_rows = [item for item in payload["streams"] if item["enabled"]]
        self.assertTrue(all(item["factor_distribution"] for item in enabled_rows))
        self.assertTrue(all(item["distribution_position"]["exact"] for item in enabled_rows))
        self.assertTrue(
            all(item["factor_distribution"]["hypothesis_count"] == 1 for item in enabled_rows)
        )
        disabled = next(
            stream_id
            for stream_id, state in config["node_streams"].items()
            if not state["enabled"]
        )
        self.assertEqual(statuses[disabled], "disabled")

    def test_edge_ledger_treats_pair_as_undirected_and_shows_disabled_streams(self) -> None:
        config = default_configuration(self.registry)
        for state in config["edge_streams"].values():
            state["enabled"] = False
        universe = pd.read_csv(
            PROJECT_ROOT / UNIVERSE_RELATIVE, sep="\t", dtype=str
        ).fillna("")
        selected = universe.loc[universe["symbol"].isin(["Prkaca", "Prkacb"])]
        self.assertEqual(len(selected), 2)
        with (
            patch.object(Path, "is_file", return_value=True),
            patch.object(evidence_inspector, "_read_configuration", return_value=config),
            patch.object(evidence_inspector, "load_registry", return_value=self.registry),
            patch.object(
                evidence_inspector, "_stream_distribution_catalog", return_value={}
            ),
            patch.object(
                evidence_inspector,
                "_matrix_symbols_and_probability",
                return_value=("Prkacb", "Prkaca", ["Prkaca", "Prkacb"], 0.5),
            ),
            patch.object(evidence_inspector, "_graph_metadata", return_value=selected),
            patch.object(
                evidence_inspector.pd,
                "read_csv",
                return_value=pd.DataFrame({"symbol": ["Prkaca", "Prkacb"]}),
            ),
        ):
            payload = inspect_edge_evidence(Path("mock-run"), "prkacb", "PRKACA")
        self.assertEqual(
            payload["canonical_undirected_pair"], ["Prkaca", "Prkacb"]
        )
        self.assertAlmostEqual(payload["stored_posterior_probability"], 0.5)
        self.assertAlmostEqual(payload["reconstructed_posterior_probability"], 0.5)
        self.assertTrue(all(item["status"] == "disabled" for item in payload["streams"]))

    def test_frontend_exposes_node_and_edge_inspector(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        for control in (
            "evidence-inspector",
            "node-evidence-form",
            "edge-evidence-form",
            "evidence-ledger-body",
            "evidence-stream-distribution-chart",
            "evidence-stream-distribution-summary",
        ):
            self.assertIn(f'id="{control}"', html)
        self.assertIn("async function inspectEvidence(kind)", javascript)
        self.assertIn("weighted_log2_odds_contribution", javascript)
        self.assertIn("reconciliation_absolute_difference", javascript)
        self.assertIn("factor_distribution", javascript)
        self.assertIn("distribution_position", javascript)


if __name__ == "__main__":
    unittest.main()
