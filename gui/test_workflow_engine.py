"""Scientific regression tests for the configurable GUI workflow engine."""

from __future__ import annotations

import re
import unittest

import numpy as np
import pandas as pd

from workflow_engine import (
    DEFAULT_SIGNAL_RELAY_CLASSES,
    PROJECT_ROOT,
    _calibration_parameter_specs,
    append_external_target,
    build_path_network_payload,
    combine_edge_factors,
    default_configuration,
    edge_stream_factor_table,
    evidence_factor_distribution_summary,
    load_registry,
    node_stream_values,
    normalize_configuration,
    scaffold_triadic_closure_factors,
    select_nodes,
)
from incremental_edge_cache import (
    ensure_incremental_pairs,
    incremental_pair_count,
    iter_incremental_pairs,
)
from ontology_directionality import (
    apply_ontology_directionality,
    build_omnipath_direction_evidence,
    build_complete_class_pair_catalog,
    load_direction_rule_catalog,
)


class WorkflowEngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = load_registry(PROJECT_ROOT)

    def test_default_configuration_is_current_workflow(self) -> None:
        config = normalize_configuration(None, self.registry)
        self.assertTrue(config["path"]["omnipath_directionality_enabled"])
        enabled_nodes = {
            key for key, state in config["node_streams"].items() if state["enabled"]
        }
        enabled_edges = {
            key for key, state in config["edge_streams"].items() if state["enabled"]
        }
        self.assertEqual(
            enabled_nodes,
            {
                "protein_abundance",
                "pc_transcript",
                "kinase_activity",
                "phosphoprotein_response",
                "pka_ca_ko_phosphoprotein_response",
                "pka_cb_ko_phosphoprotein_response",
                "collecting_duct_proteome_ccd",
                "collecting_duct_proteome_omcd",
                "collecting_duct_proteome_imcd",
                "collecting_duct_rna_ccd",
                "collecting_duct_rna_omcd",
                "collecting_duct_rna_imcd",
            },
        )
        self.assertEqual(
            {
                stream_id: state["tq_multiplier"]
                for stream_id, state in config["node_streams"].items()
            },
            {
                "protein_abundance": 0.1,
                "pc_transcript": 0.6,
                "kinase_activity": 0.1,
                "phosphoprotein_response": 0.1,
                "pka_ca_ko_phosphoprotein_response": 0.1,
                "pka_cb_ko_phosphoprotein_response": 0.05,
                "collecting_duct_proteome_ccd": 0.1,
                "collecting_duct_proteome_omcd": 0.1,
                "collecting_duct_proteome_imcd": 0.1,
                "collecting_duct_rna_ccd": 0.5,
                "collecting_duct_rna_omcd": 0.5,
                "collecting_duct_rna_imcd": 0.5,
            },
        )
        for stream_id, state in config["edge_streams"].items():
            definition = next(
                item for item in self.registry["edge_streams"]
                if item["id"] == stream_id
            )
            if definition["normalization"].get("user_control", True):
                expected = 0.5 if stream_id in {"hpa_primary", "hpa_high_confidence"} else 1.0
                self.assertEqual(state["tq_multiplier"], expected)
                if stream_id in {"hpa_primary", "hpa_high_confidence"}:
                    self.assertEqual(state["preferred_tq_multiplier"], 0.5)
            else:
                self.assertNotIn("tq_multiplier", state)
        self.assertEqual(
            enabled_edges,
            {
                "mpkccd_localization",
                "kinase_predictor",
                "string_v12",
                "hpa_primary",
                "omnipath_core",
                "stitch_secondary_messenger",
                "scaffold_triadic_closure",
            },
        )
        self.assertEqual(
            tuple(config["path"]["allowed_intermediate_classes"]),
            DEFAULT_SIGNAL_RELAY_CLASSES,
        )
        self.assertFalse(config["edge_integration"]["penalize_unsupported"])
        self.assertEqual(
            config["edge_integration"]["unsupported_bayes_factor"], 0.5
        )
        self.assertFalse(
            config["edge_integration"]["continuous_negative_evidence"]
        )
        self.assertEqual(
            config["edge_integration"]["continuous_bayes_factor_floor"], 1e-6
        )
        self.assertTrue(config["path"]["ontology_directionality_enabled"])
        self.assertNotIn(
            "adaptor_scaffold", config["path"]["allowed_intermediate_classes"]
        )
        self.assertEqual(
            config["temporal_validation"],
            {
                "enabled": False,
                "prior_df": 4.0,
                "alpha": 0.05,
                "monte_carlo_draws": 2000,
                "random_seed": 0,
                "p_adjust_method": "within_gene_bonferroni",
                "minimum_scored_nodes": 3,
            },
        )
        self.assertEqual(
            config["edge_streams"]["scaffold_triadic_closure"]["parameters"],
            {"anchor_probability_cutoff": 0.9, "closure_likelihood": 0.9},
        )
        self.assertEqual(
            config["calibration"],
            {
                "enabled": False,
                "known_nodes": [],
                "known_edges": [],
                "regularization_strength": 0.1,
                "multistart_count": 2,
            },
        )
        self.assertTrue(config["node_integration"]["penalize_unobserved"])
        self.assertEqual(
            config["node_integration"]["unobserved_bayes_factor"], 0.5
        )
        self.assertFalse(
            config["node_integration"]["continuous_negative_evidence"]
        )
        self.assertEqual(
            config["node_integration"]["continuous_bayes_factor_floor"], 1e-6
        )

    def test_exact_half_priors_are_valid_html_values(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        for field_id in ("node-prior", "edge-prior"):
            match = re.search(rf'<input id="{field_id}"[^>]*>', html)
            self.assertIsNotNone(match, field_id)
            tag = match.group(0)
            self.assertIn('min="0.000001"', tag)
            self.assertIn('max="0.999999"', tag)
            self.assertIn('step="any"', tag)
        nondetection = re.search(r'<input id="unobserved-bf"[^>]*>', html)
        self.assertIsNotNone(nondetection)
        self.assertIn('min="0.000001"', nondetection.group(0))
        self.assertIn('max="1"', nondetection.group(0))

    def test_gui_exposes_cooperative_cancellation(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        self.assertIn('id="cancel-job-button"', html)
        self.assertIn("async function cancelRun()", javascript)
        self.assertIn("/cancel`, { method: \"POST\" }", javascript)
        self.assertIn('job.status === "cancelled"', javascript)
        self.assertIn('/api/jobs/active', javascript)

    def test_gui_exposes_merged_interactive_path_network(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        for control in (
            "path-network-section",
            "path-network-limit",
            "path-network-relayout",
            "path-network-svg",
            "path-network-detail",
        ):
            self.assertIn(f'id="{control}"', html)
        self.assertIn("function renderPathNetwork(network)", javascript)
        self.assertIn("function layoutPathNetwork", javascript)
        self.assertIn("directionality === \"uniquely_directed\"", javascript)
        self.assertIn("Math.log1p(absoluteLogOdds)", javascript)
        self.assertIn('markerUnits: "userSpaceOnUse"', javascript)
        self.assertIn("function networkEvidenceColor", javascript)
        self.assertIn("const strokeWidth = 2.6", javascript)
        self.assertIn('inspectEvidence("node")', javascript)
        self.assertIn('inspectEvidence("edge")', javascript)

    def test_gui_renders_calibrated_parameters_and_factor_positions(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        for control in (
            "calibration-parameters",
            "node-calibration-parameters-body",
            "edge-calibration-parameters-body",
            "evidence-stream-distribution-chart",
            "evidence-stream-distribution-summary",
        ):
            self.assertIn(f'id="{control}"', html)
        self.assertIn("function renderCalibrationParameters", javascript)
        self.assertIn("function renderEvidenceFactorDistribution", javascript)
        self.assertIn("distribution_position", javascript)

    def test_evidence_factor_distribution_is_neutral_centered_and_compact(self) -> None:
        summary = evidence_factor_distribution_summary(
            np.asarray([0.5, 1.0, 1.0, 1.5, 2.0]),
            distribution_scope="test hypotheses",
            bin_count=8,
        )
        self.assertEqual(summary["hypothesis_count"], 5)
        self.assertEqual(summary["refuting_count"], 1)
        self.assertEqual(summary["neutral_count"], 2)
        self.assertEqual(summary["supporting_count"], 2)
        self.assertEqual(sum(summary["bin_counts"]), 5)
        self.assertEqual(summary["distribution_scope"], "test hypotheses")
        self.assertEqual(len(summary["quantile_percentages"]), 101)
        self.assertEqual(sum(item["count"] for item in summary["value_counts"]), 5)

    def test_gui_distinguishes_server_disconnect_from_bad_configuration(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        self.assertIn('id="error-title"', html)
        self.assertIn('id="reconnect-server"', html)
        self.assertIn("function showWorkflowError", javascript)
        self.assertIn("Local analysis server is not running", javascript)
        self.assertIn("failed to fetch", javascript.casefold())
        self.assertIn("window.location.reload()", javascript)

    def test_gui_exposes_optional_temporal_path_validation(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        for control in (
            "temporal-enabled",
            "temporal-prior-df",
            "temporal-alpha",
            "temporal-draws",
            "temporal-seed",
            "temporal-p-adjust",
            "temporal-min-scored",
        ):
            self.assertIn(f'id="{control}"', html)
        self.assertIn("temporal_validation:", javascript)
        self.assertIn("preview.temporal_validation", javascript)
        self.assertIn("Primary Bayesian ranks were preserved", javascript)
        for control in ("temporal-prior-df", "temporal-alpha"):
            match = re.search(rf'<input id="{control}"[^>]*>', html)
            self.assertIsNotNone(match, control)
            self.assertIn('step="any"', match.group(0))

    def test_gui_exposes_positive_control_parameter_calibration(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        for control in (
            "calibration-enabled",
            "calibration-known-nodes",
            "calibration-known-edges",
            "calibration-lambda",
            "calibration-result",
        ):
            self.assertIn(f'id="{control}"', html)
        self.assertNotIn('id="calibration-phospho-tq"', html)
        self.assertIn("calibration:", javascript)
        self.assertIn("stream-preferred-tq", javascript)
        self.assertIn("preferred_tq_multiplier", javascript)
        self.assertIn("setCalibrationControls", javascript)
        self.assertIn("preview.calibration", javascript)
        self.assertIn("Unknown hypotheses were not treated as negatives", javascript)

    def test_calibration_targets_are_normalized_without_negative_labels(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["calibration"].update(
            {
                "enabled": True,
                "known_nodes": ["Prkaca", "prkaca", " Pde4d "],
                "known_edges": ["Prkaca,Aqp2", ["Aqp2", "Prkaca"]],
            }
        )
        config = normalize_configuration(supplied, self.registry)
        self.assertEqual(config["calibration"]["known_nodes"], ["Prkaca", "Pde4d"])
        self.assertEqual(config["calibration"]["known_edges"], [["Prkaca", "Aqp2"]])

    def test_every_primary_stream_has_an_independent_preferred_tq(self) -> None:
        supplied = default_configuration(self.registry)
        for group in ("node_streams", "edge_streams"):
            definitions = {item["id"]: item for item in self.registry[group]}
            for stream_id, state in supplied[group].items():
                definition = definitions[stream_id]
                if definition.get("normalization", {}).get("user_control", True):
                    self.assertIn("preferred_tq_multiplier", state, stream_id)
                else:
                    self.assertNotIn("preferred_tq_multiplier", state, stream_id)

        supplied["node_streams"]["protein_abundance"]["preferred_tq_multiplier"] = 2.0
        supplied["node_streams"]["kinase_activity"]["preferred_tq_multiplier"] = 0.2
        supplied["edge_streams"]["string_v12"]["preferred_tq_multiplier"] = 0.5
        config = normalize_configuration(supplied, self.registry)
        self.assertEqual(
            config["node_streams"]["protein_abundance"]["preferred_tq_multiplier"],
            2.0,
        )
        self.assertEqual(
            config["node_streams"]["kinase_activity"]["preferred_tq_multiplier"],
            0.2,
        )
        self.assertEqual(
            config["edge_streams"]["string_v12"]["preferred_tq_multiplier"],
            0.5,
        )

        node_specs = _calibration_parameter_specs(
            self.registry, config, "node_streams"
        )
        edge_specs = _calibration_parameter_specs(
            self.registry, config, "edge_streams"
        )
        preferred = {
            row["key"]: row["preferred"] for row in node_specs + edge_specs
        }
        self.assertEqual(preferred["protein_abundance:tq_multiplier"], 2.0)
        self.assertEqual(preferred["kinase_activity:tq_multiplier"], 0.2)
        self.assertEqual(preferred["string_v12:tq_multiplier"], 0.5)

    def test_legacy_global_phosphoproteomic_preference_is_migrated(self) -> None:
        supplied = {"calibration": {"phosphoproteomic_preferred_tq_multiplier": 0.3}}
        config = normalize_configuration(supplied, self.registry)
        self.assertNotIn(
            "phosphoproteomic_preferred_tq_multiplier", config["calibration"]
        )
        for group in ("node_streams", "edge_streams"):
            for stream_id, state in config[group].items():
                if stream_id in {
                    "kinase_activity",
                    "phosphoprotein_response",
                    "pka_ca_ko_phosphoprotein_response",
                    "pka_cb_ko_phosphoprotein_response",
                    "kinase_predictor",
                }:
                    self.assertEqual(state["preferred_tq_multiplier"], 0.3)

    def test_run_button_waits_for_configuration_initialization(self) -> None:
        html = (PROJECT_ROOT / "gui/web/index.html").read_text(encoding="utf-8")
        javascript = (PROJECT_ROOT / "gui/web/app.js").read_text(encoding="utf-8")
        match = re.search(r'<button id="run-button"[^>]*>', html)
        self.assertIsNotNone(match)
        self.assertIn("disabled", match.group(0))
        self.assertIn('$("run-button").disabled = false;', javascript)

    def test_ontology_directionality_orients_only_unambiguous_role_pairs(self) -> None:
        symbols = ["Ligand", "Receptor", "Kinase", "Binder", "Broad", "MultiA", "MultiB"]
        values = np.full((len(symbols), len(symbols)), 0.9, dtype=float)
        np.fill_diagonal(values, 0.0)
        matrix = pd.DataFrame(values, index=symbols, columns=symbols)
        metadata = pd.DataFrame(
            {
                "symbol": symbols,
                "classes": [
                    "ligand",
                    "receptor",
                    "kinase",
                    "kinase_phosphatase_binding",
                    "signaling_process",
                    "kinase;kinase_phosphatase_binding",
                    "kinase;kinase_phosphatase_binding",
                ],
            }
        )
        catalog = load_direction_rule_catalog(
            known_classes=[item["id"] for item in self.registry["path_ontology_classes"]]
        )
        directed, audit, summary = apply_ontology_directionality(
            matrix,
            metadata,
            catalog,
            audit_probability_cutoff=0.5,
            edge_output_cutoff=0.5,
            path_probability_cutoff=0.5,
        )
        self.assertEqual(directed.loc["Ligand", "Receptor"], 0.9)
        self.assertEqual(directed.loc["Receptor", "Ligand"], 0.0)
        self.assertEqual(directed.loc["Kinase", "Binder"], 0.9)
        self.assertEqual(directed.loc["Binder", "Kinase"], 0.0)
        self.assertEqual(directed.loc["Kinase", "Broad"], 0.9)
        self.assertEqual(directed.loc["Broad", "Kinase"], 0.9)
        self.assertEqual(directed.loc["MultiA", "MultiB"], 0.9)
        self.assertEqual(directed.loc["MultiB", "MultiA"], 0.9)
        conflict = audit.loc[
            audit["node_a"].eq("MultiA") & audit["node_b"].eq("MultiB")
        ].iloc[0]
        self.assertEqual(
            conflict["directionality_status"], "unresolved_conflicting_directions"
        )
        self.assertGreater(summary["edge_output_graph"]["uniquely_oriented_edge_count"], 0)
        class_pairs = build_complete_class_pair_catalog(
            catalog,
            [item["id"] for item in self.registry["path_ontology_classes"]],
        )
        self.assertEqual(len(class_pairs), 153)

    def test_omnipath_directions_orient_only_unidirectional_pairs(self) -> None:
        symbols = ["Ligand", "Receptor", "A", "B", "C", "D"]
        raw = pd.DataFrame(
            [
                {
                    "source_genesymbol": "Receptor",
                    "target_genesymbol": "Ligand",
                    "is_directed": "True",
                    "consensus_direction": "False",
                    "sources": "KEA",
                    "references": "KEA:2",
                },
                {
                    "source_genesymbol": "A",
                    "target_genesymbol": "B",
                    "is_directed": "True",
                    "consensus_direction": "False",
                    "sources": "HPRD",
                    "references": "HPRD:3",
                },
                {
                    "source_genesymbol": "C",
                    "target_genesymbol": "D",
                    "is_directed": "True",
                    "consensus_direction": "False",
                    "sources": "X",
                    "references": "X:4",
                },
                {
                    "source_genesymbol": "D",
                    "target_genesymbol": "C",
                    "is_directed": "True",
                    "consensus_direction": "False",
                    "sources": "Y",
                    "references": "Y:5",
                },
            ]
        )
        omnipath = build_omnipath_direction_evidence(raw, symbols)
        self.assertEqual(len(omnipath), 3)
        ab = omnipath.loc[
            omnipath["node_a"].eq("A") & omnipath["node_b"].eq("B")
        ].iloc[0]
        self.assertTrue(ab["omnipath_supports_a_to_b"])
        self.assertFalse(ab["omnipath_supports_b_to_a"])
        cd = omnipath.loc[
            omnipath["node_a"].eq("C") & omnipath["node_b"].eq("D")
        ].iloc[0]
        self.assertTrue(cd["omnipath_bidirectional"])

        values = np.full((len(symbols), len(symbols)), 0.9, dtype=float)
        np.fill_diagonal(values, 0.0)
        matrix = pd.DataFrame(values, index=symbols, columns=symbols)
        metadata = pd.DataFrame(
            {
                "symbol": symbols,
                "classes": ["ligand", "receptor", "", "", "", ""],
            }
        )
        catalog = load_direction_rule_catalog(
            known_classes=[item["id"] for item in self.registry["path_ontology_classes"]]
        )
        directed, audit, summary = apply_ontology_directionality(
            matrix,
            metadata,
            catalog,
            audit_probability_cutoff=0.5,
            edge_output_cutoff=0.5,
            path_probability_cutoff=0.5,
            omnipath_directions=omnipath,
        )
        self.assertEqual(directed.loc["A", "B"], 0.9)
        self.assertEqual(directed.loc["B", "A"], 0.0)
        self.assertEqual(directed.loc["C", "D"], 0.9)
        self.assertEqual(directed.loc["D", "C"], 0.9)
        self.assertEqual(directed.loc["Ligand", "Receptor"], 0.9)
        self.assertEqual(directed.loc["Receptor", "Ligand"], 0.0)
        ligand_receptor = audit.loc[
            audit["node_a"].eq("Ligand") & audit["node_b"].eq("Receptor")
        ].iloc[0]
        self.assertEqual(
            ligand_receptor["directionality_status"],
            "oriented_a_to_b",
        )
        self.assertEqual(
            ligand_receptor["direction_evidence_sources"],
            "ontology_precedence_over_omnipath_conflict",
        )
        edge_summary = summary["edge_output_graph"]
        self.assertEqual(edge_summary["uniquely_oriented_by_omnipath_only_count"], 1)
        self.assertEqual(edge_summary["unresolved_conflicting_direction_count"], 1)
        self.assertEqual(
            edge_summary["ontology_precedence_over_opposing_omnipath_count"], 1
        )

    def test_path_network_merges_paths_and_marks_only_constrained_arrows(self) -> None:
        path_rows = pd.DataFrame(
            [
                {"rank": 1, "path_symbols": "A -> B -> D"},
                {"rank": 2, "path_symbols": "A -> C -> B -> D"},
            ]
        )
        path_edges = pd.DataFrame(
            [
                {"path_rank": 1, "step": 1, "source_symbol": "A", "target_symbol": "B", "edge_probability": 0.9},
                {"path_rank": 1, "step": 2, "source_symbol": "B", "target_symbol": "D", "edge_probability": 0.8},
                {"path_rank": 2, "step": 1, "source_symbol": "A", "target_symbol": "C", "edge_probability": 0.7},
                {"path_rank": 2, "step": 2, "source_symbol": "C", "target_symbol": "B", "edge_probability": 0.75},
                {"path_rank": 2, "step": 3, "source_symbol": "B", "target_symbol": "D", "edge_probability": 0.8},
            ]
        )
        metadata = pd.DataFrame(
            {
                "symbol": ["A", "B", "C", "D"],
                "display_symbol": ["A", "B", "C", "D"],
                "name": ["start", "hub", "branch", "target"],
                "classes": ["kinase", "signaling_process", "receptor", "external_target"],
                "node_type": ["protein"] * 4,
                "gui_posterior": [0.9, 0.7, 0.6, np.nan],
            }
        )
        propagation = pd.DataFrame(
            np.zeros((4, 4)), index=["A", "B", "C", "D"], columns=["A", "B", "C", "D"]
        )
        for left, right, probability in (
            ("A", "B", 0.9),
            ("A", "C", 0.7),
            ("B", "C", 0.75),
            ("B", "D", 0.8),
        ):
            propagation.loc[left, right] = probability
            propagation.loc[right, left] = probability
        propagation.loc["B", "A"] = 0.0

        payload = build_path_network_payload(
            path_rows,
            path_edges,
            metadata,
            propagation,
            directed_graph=True,
        )
        self.assertEqual(payload["visualized_path_count"], 2)
        self.assertEqual(len(payload["nodes"]), 4)
        self.assertEqual(len(payload["edges"]), 4)
        ab = next(edge for edge in payload["edges"] if edge["id"] == "A--B")
        self.assertEqual(ab["directionality"], "uniquely_directed")
        self.assertEqual((ab["source"], ab["target"]), ("A", "B"))
        bd = next(edge for edge in payload["edges"] if edge["id"] == "B--D")
        self.assertEqual(bd["directionality"], "unresolved_bidirectional")
        self.assertEqual(bd["path_ranks"], [1, 2])
        node_b = next(node for node in payload["nodes"] if node["id"] == "B")
        self.assertEqual(node_b["path_count"], 2)
        self.assertAlmostEqual(node_b["mean_path_position"], (0.5 + 2 / 3) / 2)
        node_d = next(node for node in payload["nodes"] if node["id"] == "D")
        self.assertIsNone(node_d["posterior_probability"])
        self.assertTrue(node_d["is_target"])

    def test_incremental_pairs_stream_without_materializing_full_graph(self) -> None:
        symbols = ["SeedA", "AddedA", "SeedB", "AddedB"]
        base = {"SeedA", "SeedB"}
        pairs = list(iter_incremental_pairs(symbols, {"AddedA", "AddedB"}))
        self.assertEqual(incremental_pair_count(symbols, base), 5)
        self.assertEqual(len(pairs), 5)
        self.assertEqual(len(set(pairs)), 5)
        self.assertNotIn(("SeedA", "SeedB"), set(pairs))

    def test_all_collecting_duct_streams_have_a_large_but_valid_graph(self) -> None:
        supplied = default_configuration(self.registry)
        collecting_duct_streams = {
            stream_id
            for stream_id in supplied["node_streams"]
            if stream_id.startswith("collecting_duct_")
        }
        self.assertEqual(len(collecting_duct_streams), 6)
        for stream_id in collecting_duct_streams:
            supplied["node_streams"][stream_id]["enabled"] = True
        supplied["edge_streams"]["scaffold_triadic_closure"]["enabled"] = True
        config = normalize_configuration(supplied, self.registry)
        _, selected, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        self.assertEqual(len(selected), 1506)
        self.assertEqual(len(selected) * (len(selected) - 1) // 2, 1_133_265)

    def test_path_ontology_class_selection_is_normalized_and_validated(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["path"]["allowed_intermediate_classes"] = [
            "second_messenger",
            "kinase",
            "ligand",
            "kinase",
        ]
        config = normalize_configuration(supplied, self.registry)
        self.assertEqual(
            config["path"]["allowed_intermediate_classes"],
            ["ligand", "kinase", "second_messenger"],
        )

        supplied["path"]["allowed_intermediate_classes"] = ["not_a_real_class"]
        with self.assertRaisesRegex(ValueError, "unknown intermediate ontology"):
            normalize_configuration(supplied, self.registry)

    def test_default_node_selection_uses_current_site_level_evidence(self) -> None:
        config = normalize_configuration(None, self.registry)
        factors, selected, summary = select_nodes(PROJECT_ROOT, self.registry, config)
        self.assertEqual(len(factors), 9170)
        self.assertEqual(len(selected), 1506)
        self.assertEqual(summary["selected_protein_count"], 1486)
        self.assertEqual(summary["incrementally_added_protein_count"], 881)
        self.assertEqual(summary["curated_second_messenger_count"], 20)
        self.assertTrue(summary["posterior_probabilities_are_independent"])
        self.assertEqual(summary["node_prior_probability"], 0.5)
        self.assertEqual(summary["node_output_probability_cutoff"], 0.5)
        self.assertTrue((factors["initial_prior"] == 0.5).all())
        self.assertTrue((factors["gui_initial_prior_probability"] == 0.5).all())
        self.assertTrue((factors["gui_posterior"] > 0.0).all())
        self.assertTrue((factors["gui_posterior"] < 1.0).all())
        self.assertEqual(
            int((factors["gui_posterior"] > 0.5).sum()),
            summary["selected_protein_count"],
        )
        self.assertGreater(float(factors["gui_posterior"].sum()), 1.0)
        distribution = summary["probability_distribution"]
        self.assertEqual(distribution["hypothesis_count"], len(factors))
        self.assertEqual(sum(distribution["bin_counts"]), len(factors))
        self.assertEqual(
            distribution["at_exact_prior_count"],
            summary["neutral_posterior_count"],
        )
        self.assertEqual(
            distribution["above_output_cutoff_count"],
            summary["selected_protein_count"],
        )
        self.assertEqual(distribution["below_prior_count"], 7681)
        for stream in summary["active_streams"]:
            factor_distribution = stream["factor_distribution"]
            self.assertEqual(factor_distribution["hypothesis_count"], len(factors))
            self.assertEqual(sum(factor_distribution["bin_counts"]), len(factors))

    def test_neutral_node_evidence_preserves_independent_half_priors(self) -> None:
        config = normalize_configuration(None, self.registry)
        for state in config["node_streams"].values():
            state["enabled"] = False

        factors, selected, summary = select_nodes(PROJECT_ROOT, self.registry, config)

        self.assertTrue(np.array_equal(
            factors["gui_posterior"].to_numpy(float),
            np.full(len(factors), 0.5),
        ))
        self.assertEqual(summary["selected_protein_count"], 0)
        self.assertEqual(summary["neutral_posterior_count"], len(factors))
        self.assertEqual(len(selected), 20)

    def test_optional_nondetection_evidence_can_lower_node_posteriors(self) -> None:
        neutral_supplied = default_configuration(self.registry)
        neutral_supplied["node_integration"]["penalize_unobserved"] = False
        neutral_config = normalize_configuration(neutral_supplied, self.registry)
        _, neutral_selected, neutral_summary = select_nodes(
            PROJECT_ROOT, self.registry, neutral_config
        )

        supplied = default_configuration(self.registry)
        supplied["node_integration"]["penalize_unobserved"] = True
        supplied["node_integration"]["unobserved_bayes_factor"] = 0.5
        config = normalize_configuration(supplied, self.registry)

        factors, selected, summary = select_nodes(
            PROJECT_ROOT, self.registry, config
        )

        self.assertTrue(summary["penalize_unobserved"])
        self.assertEqual(summary["unobserved_bayes_factor"], 0.5)
        self.assertGreater(summary["total_unobserved_penalties_applied"], 0)
        self.assertGreater(summary["candidates_below_prior"], 0)
        self.assertEqual(
            summary["probability_distribution"]["below_prior_count"],
            summary["candidates_below_prior"],
        )
        self.assertLess(summary["posterior_minimum"], 0.5)
        self.assertLess(
            summary["selected_protein_count"],
            neutral_summary["selected_protein_count"],
        )
        self.assertLess(len(selected), len(neutral_selected))
        self.assertEqual(len(selected), summary["selected_protein_count"] + 20)

        active_ids = [stream["id"] for stream in summary["active_streams"]]
        effective_bfs = factors[
            [f"gui_{stream_id}_bayes_factor" for stream_id in active_ids]
        ].to_numpy(float)
        expected = 1.0 / (1.0 + np.exp(-np.log(effective_bfs).sum(axis=1)))
        self.assertTrue(
            np.allclose(expected, factors["gui_posterior"].to_numpy(float))
        )

    def test_kinase_nondetection_penalty_is_scope_aware(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["node_streams"].values():
            state["enabled"] = False
        supplied["node_streams"]["kinase_activity"]["enabled"] = True
        supplied["node_integration"]["penalize_unobserved"] = True
        config = normalize_configuration(supplied, self.registry)

        factors, _, summary = select_nodes(PROJECT_ROOT, self.registry, config)
        eligible = factors["gui_kinase_activity_negative_evidence_eligible"]
        observed = factors["gui_kinase_activity_observed"]
        applied = factors["gui_kinase_activity_unobserved_penalty_applied"]
        effective = factors["gui_kinase_activity_bayes_factor"]

        self.assertTrue((~applied[~eligible]).all())
        self.assertTrue((effective[~eligible] == 1.0).all())
        self.assertTrue((effective[eligible & ~observed] == 0.5).all())
        self.assertEqual(
            summary["active_streams"][0]["unobserved_penalties_applied"],
            int((eligible & ~observed).sum()),
        )

    def test_continuous_node_evidence_scores_weak_values_and_zero_nondetections(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["node_streams"].values():
            state["enabled"] = False
        supplied["node_streams"]["pc_transcript"]["enabled"] = True
        supplied["node_streams"]["pc_transcript"][
            "continuous_negative_evidence"
        ] = True
        supplied["node_integration"]["continuous_bayes_factor_floor"] = 1e-5
        supplied["node_integration"]["penalize_unobserved"] = True
        config = normalize_configuration(supplied, self.registry)
        self.assertTrue(config["node_integration"]["penalize_unobserved"])

        factors, _, summary = select_nodes(PROJECT_ROOT, self.registry, config)
        observed = factors["gui_pc_transcript_observed"]
        bayes_factors = factors["gui_pc_transcript_bayes_factor"]
        eligible = factors["gui_pc_transcript_negative_evidence_eligible"]

        self.assertTrue((bayes_factors[eligible & ~observed] == 1e-5).all())
        self.assertTrue((bayes_factors[observed] < 1.0).any())
        self.assertTrue((bayes_factors[observed] > 1.0).any())
        self.assertGreater(summary["total_continuous_negative_applications"], 0)
        self.assertGreater(summary["candidates_below_prior"], 0)

    def test_continuous_kinase_evidence_leaves_nonkinases_neutral(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["node_streams"].values():
            state["enabled"] = False
        supplied["node_streams"]["kinase_activity"]["enabled"] = True
        supplied["node_streams"]["kinase_activity"][
            "continuous_negative_evidence"
        ] = True
        config = normalize_configuration(supplied, self.registry)

        factors, _, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        eligible = factors["gui_kinase_activity_negative_evidence_eligible"]
        bayes_factors = factors["gui_kinase_activity_bayes_factor"]
        self.assertTrue((bayes_factors[~eligible] == 1.0).all())
        self.assertTrue((bayes_factors[eligible] < 1.0).any())

    def test_nondetection_bayes_factor_must_be_positive_and_at_most_one(self) -> None:
        for invalid in (0.0, 1.000001):
            supplied = default_configuration(self.registry)
            supplied["node_integration"]["unobserved_bayes_factor"] = invalid
            with self.assertRaisesRegex(ValueError, "unobserved node Bayes factor"):
                normalize_configuration(supplied, self.registry)

    def test_protein_and_pc_only_configuration_is_valid_subset(self) -> None:
        supplied = default_configuration(self.registry)
        for stream_id in supplied["node_streams"]:
            supplied["node_streams"][stream_id]["enabled"] = stream_id in {
                "protein_abundance", "pc_transcript"
            }
        supplied["node_integration"]["penalize_unobserved"] = False
        config = normalize_configuration(supplied, self.registry)
        _, selected, summary = select_nodes(PROJECT_ROOT, self.registry, config)
        self.assertEqual(len(selected), 1254)
        self.assertEqual(summary["selected_protein_count"], 1234)
        self.assertEqual(
            {item["id"] for item in summary["active_streams"]},
            {"protein_abundance", "pc_transcript"},
        )

    def test_default_edge_combination_reconstructs_current_graph(self) -> None:
        config = normalize_configuration(None, self.registry)
        _, selected, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        matrix, summary = combine_edge_factors(
            PROJECT_ROOT,
            self.registry,
            config,
            selected["symbol"].tolist(),
        )
        values = matrix.to_numpy(float)
        self.assertEqual(matrix.shape, (1506, 1506))
        self.assertTrue(np.array_equal(values, values.T))
        self.assertTrue(np.array_equal(np.diag(values), np.zeros(1506)))
        self.assertGreater(summary["pairs_above_output_cutoff"], 169414)
        distribution = summary["probability_distribution"]
        self.assertEqual(distribution["hypothesis_count"], 1_133_265)
        self.assertEqual(sum(distribution["bin_counts"]), 1_133_265)
        self.assertEqual(
            distribution["at_exact_prior_count"], summary["pairs_at_exact_prior"]
        )
        self.assertEqual(
            distribution["above_output_cutoff_count"],
            summary["pairs_above_output_cutoff"],
        )
        for stream in summary["active_streams"]:
            factor_distribution = stream["factor_distribution"]
            self.assertEqual(factor_distribution["hypothesis_count"], 1_133_265)
            self.assertEqual(sum(factor_distribution["bin_counts"]), 1_133_265)

    def test_optional_edge_absence_penalty_lowers_only_eligible_pairs(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["edge_streams"].values():
            state["enabled"] = False
        supplied["edge_streams"]["stitch_secondary_messenger"]["enabled"] = True
        supplied["edge_integration"]["penalize_unsupported"] = True
        supplied["edge_integration"]["unsupported_bayes_factor"] = 0.5
        config = normalize_configuration(supplied, self.registry)
        symbols = ["Prkar2a", "Actn1", "SM_CAMP"]
        matrix, summary = combine_edge_factors(
            PROJECT_ROOT,
            self.registry,
            config,
            symbols,
        )
        values = matrix.to_numpy(float)
        self.assertTrue(np.array_equal(values, values.T))
        self.assertEqual(matrix.loc["Prkar2a", "Actn1"], 0.5)
        self.assertGreater(matrix.loc["Prkar2a", "SM_CAMP"], 0.5)
        self.assertLess(matrix.loc["Actn1", "SM_CAMP"], 0.5)
        self.assertEqual(summary["pairs_below_prior"], 1)
        self.assertEqual(summary["negative_penalty_applications"], 1)

    def test_continuous_edge_evidence_scores_low_and_zero_support(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["edge_streams"].values():
            state["enabled"] = False
        supplied["edge_streams"]["stitch_secondary_messenger"]["enabled"] = True
        supplied["edge_streams"]["stitch_secondary_messenger"][
            "continuous_negative_evidence"
        ] = True
        supplied["edge_integration"]["continuous_bayes_factor_floor"] = 1e-5
        supplied["edge_integration"]["penalize_unsupported"] = True
        config = normalize_configuration(supplied, self.registry)
        self.assertTrue(config["edge_integration"]["penalize_unsupported"])

        symbols = ["Prkar2a", "Actn1", "SM_CAMP"]
        matrix, summary = combine_edge_factors(
            PROJECT_ROOT,
            self.registry,
            config,
            symbols,
        )
        self.assertEqual(matrix.loc["Prkar2a", "Actn1"], 0.5)
        self.assertGreater(matrix.loc["Prkar2a", "SM_CAMP"], 0.5)
        self.assertAlmostEqual(
            matrix.loc["Actn1", "SM_CAMP"],
            1e-5 / (1.0 + 1e-5),
        )
        self.assertEqual(
            summary["continuous_negative_streams"],
            ["stitch_secondary_messenger"],
        )
        self.assertEqual(summary["negative_penalty_applications"], 1)

    def test_unsupported_edge_bayes_factor_is_validated(self) -> None:
        for invalid in (0.0, 1.000001):
            supplied = default_configuration(self.registry)
            supplied["edge_integration"]["unsupported_bayes_factor"] = invalid
            with self.assertRaisesRegex(ValueError, "unsupported edge Bayes factor"):
                normalize_configuration(supplied, self.registry)

    def test_external_target_edge_penalty_uses_source_scope(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["edge_streams"].values():
            state["enabled"] = False
        supplied["edge_streams"]["stitch_secondary_messenger"]["enabled"] = True
        supplied["edge_integration"]["penalize_unsupported"] = True
        config = normalize_configuration(supplied, self.registry)
        symbols = ["SM_CAMP", "SM_AA", "SM_PIP3"]
        matrix = pd.DataFrame(
            np.zeros((len(symbols), len(symbols))),
            index=symbols,
            columns=symbols,
        )
        vector = pd.read_csv(
            PROJECT_ROOT
            / "results/path_finding/target_extensions/Aqp2/target_adjacency_vector.tsv",
            sep="\t",
        )
        extended, warnings = append_external_target(
            PROJECT_ROOT,
            matrix,
            "Aqp2",
            vector,
            self.registry,
            config,
        )
        self.assertGreater(extended.loc["SM_CAMP", "Aqp2"], 0.5)
        self.assertEqual(extended.loc["SM_AA", "Aqp2"], 1.0 / 3.0)
        self.assertEqual(extended.loc["SM_PIP3", "Aqp2"], 0.5)
        self.assertTrue(any("applied unsupported-pair BF" in item for item in warnings))

    def test_hpa_alternatives_cannot_be_enabled_together(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["edge_streams"]["hpa_high_confidence"]["enabled"] = True
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            normalize_configuration(supplied, self.registry)

    def test_zero_weight_stream_does_not_count_as_effective(self) -> None:
        supplied = default_configuration(self.registry)
        for state in supplied["node_streams"].values():
            state["enabled"] = True
            state["weight"] = 0
        with self.assertRaisesRegex(ValueError, "positive weight"):
            normalize_configuration(supplied, self.registry)

    def test_path_limit_allows_500_but_rejects_more(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["path"]["top_k"] = 500
        config = normalize_configuration(supplied, self.registry)
        self.assertEqual(config["path"]["top_k"], 500)
        supplied["path"]["top_k"] = 501
        with self.assertRaisesRegex(ValueError, "top paths"):
            normalize_configuration(supplied, self.registry)

    def test_temporal_validation_configuration_is_normalized_and_requires_paths(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["temporal_validation"].update(
            {
                "enabled": True,
                "prior_df": 8,
                "alpha": 0.1,
                "monte_carlo_draws": 500,
                "random_seed": 22,
                "p_adjust_method": "none",
                "minimum_scored_nodes": 4,
            }
        )
        config = normalize_configuration(supplied, self.registry)
        self.assertEqual(
            config["temporal_validation"], supplied["temporal_validation"]
        )

        supplied["path"]["enabled"] = False
        with self.assertRaisesRegex(ValueError, "requires path finding"):
            normalize_configuration(supplied, self.registry)

        supplied["path"]["enabled"] = True
        supplied["temporal_validation"]["p_adjust_method"] = "unknown"
        with self.assertRaisesRegex(ValueError, "p-adjustment"):
            normalize_configuration(supplied, self.registry)

    def test_each_node_stream_has_an_independent_tq_control(self) -> None:
        config = normalize_configuration(None, self.registry)
        factors, _, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        for definition in self.registry["node_streams"]:
            stream_id = definition["id"]
            state = dict(config["node_streams"][stream_id])
            state["tq_multiplier"] = 0.8
            rescored = node_stream_values(
                PROJECT_ROOT, factors, definition, state
            )
            original = factors[definition["column"]].fillna(
                definition["neutral_value"]
            ).to_numpy(float)
            self.assertFalse(
                np.allclose(rescored, original, rtol=0, atol=1e-14),
                stream_id,
            )
            self.assertTrue(np.isfinite(rescored).all(), stream_id)
            self.assertTrue((rescored > 0).all(), stream_id)

    def test_collecting_duct_streams_are_neutral_for_nondetection(self) -> None:
        config = normalize_configuration(None, self.registry)
        factors, _, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        definitions = {
            item["id"]: item
            for item in self.registry["node_streams"]
            if item["id"].startswith("collecting_duct_")
        }
        self.assertEqual(len(definitions), 6)
        self.assertTrue(
            all(config["node_streams"][stream_id]["enabled"] for stream_id in definitions)
        )
        self.assertEqual(
            {
                item["dependence_group"]
                for item in definitions.values()
            },
            {"ktea_collecting_duct_proteome", "mouse_renal_tubule_rna_seq"},
        )
        for stream_id, definition in definitions.items():
            original = pd.to_numeric(
                factors[definition["column"]], errors="raise"
            ).to_numpy(float)
            observed = factors[definition["observed_column"]].astype(bool).to_numpy()
            self.assertTrue(np.all(original[~observed] == 0.5), stream_id)
            self.assertTrue(np.all((original >= 0.5) & (original <= 1.0)), stream_id)
            sensitive_state = dict(config["node_streams"][stream_id])
            sensitive_state["tq_multiplier"] = 0.8
            sensitive = node_stream_values(
                PROJECT_ROOT, factors, definition, sensitive_state
            )
            self.assertTrue(np.all(sensitive >= original - 1e-14), stream_id)
            self.assertTrue(np.any(sensitive > original + 1e-14), stream_id)

    def test_selective_pka_subunit_ko_streams_are_separate(self) -> None:
        stream_ids = {
            "pka_ca_ko_phosphoprotein_response",
            "pka_cb_ko_phosphoprotein_response",
        }
        definitions = {
            item["id"]: item
            for item in self.registry["node_streams"]
            if item["id"] in stream_ids
        }
        self.assertEqual(set(definitions), stream_ids)
        self.assertEqual(
            {item["dependence_group"] for item in definitions.values()},
            {"pka_subunit_ko_phosphoproteomics"},
        )

        supplied = default_configuration(self.registry)
        self.assertTrue(
            all(supplied["node_streams"][stream_id]["enabled"] for stream_id in stream_ids)
        )
        for stream_id in supplied["node_streams"]:
            supplied["node_streams"][stream_id]["enabled"] = stream_id in stream_ids
        config = normalize_configuration(supplied, self.registry)
        factors, selected, summary = select_nodes(PROJECT_ROOT, self.registry, config)

        self.assertEqual(len(factors), 9170)
        self.assertEqual(len(selected), summary["selected_protein_count"] + 20)
        expected_counts = {
            "pka_ca_ko_phosphoprotein_response": (776, 708),
            "pka_cb_ko_phosphoprotein_response": (776, 748),
        }
        active = {item["id"]: item for item in summary["active_streams"]}
        self.assertEqual(set(active), stream_ids)
        for stream_id, (observed_count, positive_count) in expected_counts.items():
            definition = definitions[stream_id]
            self.assertEqual(
                int(factors[definition["observed_column"]].sum()), observed_count
            )
            self.assertEqual(active[stream_id]["observed_eligible_candidates"], observed_count)
            self.assertEqual(active[stream_id]["positive_support_candidates"], positive_count)
            self.assertTrue(
                (factors[f"gui_{stream_id}_source_bayes_factor"] >= 1.0).all()
            )
            self.assertTrue(
                (factors[f"gui_{stream_id}_source_bayes_factor"] <= 2.0).all()
            )

    def test_phosphoprotein_streams_use_one_site_tq_and_maximum_site_bf(self) -> None:
        double_nodes = pd.read_csv(
            PROJECT_ROOT
            / "results/phosphoprotein_evidence/"
            "node_selection_protein_pc_phosphosite_posterior.tsv",
            sep="\t",
        )
        double_sites = pd.read_csv(
            PROJECT_ROOT / "results/phosphoprotein_evidence/phosphosite_audit.tsv",
            sep="\t",
        )
        observed_double = double_nodes.loc[
            double_nodes["phosphosite_evidence_observed"]
        ].set_index("gene_symbol")
        expected_double = (
            double_sites.loc[double_sites["in_signaling_universe"]]
            .groupby("gene_symbol")["site_bayes_factor"]
            .max()
            .reindex(observed_double.index)
        )
        self.assertEqual(observed_double["site_T_q"].nunique(), 1)
        self.assertNotIn("matched_T_q", observed_double.columns)
        self.assertTrue(
            np.allclose(
                observed_double["relative_multiplier"],
                expected_double,
                rtol=0,
                atol=1e-14,
            )
        )

        subunit_nodes = pd.read_csv(
            PROJECT_ROOT
            / "data/node_selection/pka_subunit_ko/processed/"
            "pka_subunit_ko_node_factors.tsv.gz",
            sep="\t",
        )
        subunit_sites = pd.read_csv(
            PROJECT_ROOT
            / "data/node_selection/pka_subunit_ko/processed/"
            "pka_subunit_ko_site_audit.tsv.gz",
            sep="\t",
        )
        for prefix in ("pka_ca_ko", "pka_cb_ko"):
            observed = subunit_nodes.loc[
                subunit_nodes[f"{prefix}_observed"]
            ].set_index("gene_symbol")
            expected = (
                subunit_sites.loc[subunit_sites["mapped_gene_symbol"].notna()]
                .groupby("mapped_gene_symbol")[f"{prefix}_site_bayes_factor"]
                .max()
                .reindex(observed.index)
            )
            self.assertEqual(observed[f"{prefix}_site_tq"].nunique(), 1)
            self.assertNotIn(f"{prefix}_matched_tq", observed.columns)
            self.assertTrue(
                np.allclose(
                    observed[f"{prefix}_phosphoprotein_factor"],
                    expected,
                    rtol=0,
                    atol=1e-14,
                ),
                prefix,
            )

    def test_each_edge_stream_can_be_rescored_from_raw_evidence(self) -> None:
        config = normalize_configuration(None, self.registry)
        _, selected, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        symbols = selected["symbol"].tolist()
        for definition in self.registry["edge_streams"]:
            if definition.get("derived"):
                continue
            state = dict(config["edge_streams"][definition["id"]])
            state["tq_multiplier"] = 0.8
            table = edge_stream_factor_table(
                PROJECT_ROOT, definition, state, symbols
            )
            self.assertGreater(len(table), 0, definition["id"])
            self.assertTrue(np.isfinite(table["bayes_factor"]).all(), definition["id"])
            self.assertTrue((table["bayes_factor"] > 0).all(), definition["id"])

    def test_scaffold_closure_is_binary_one_pass_protein_evidence(self) -> None:
        symbols = ["A", "B", "C", "S", "cyclic AMP"]
        physical_scores = np.zeros((5, 5), dtype=float)
        physical_scores[0, 3] = physical_scores[3, 0] = 0.98
        physical_scores[1, 3] = physical_scores[3, 1] = 0.97
        physical_scores[2, 3] = physical_scores[3, 2] = 0.91
        metadata = pd.DataFrame(
            {
                "symbol": symbols,
                "classes": ["kinase", "phosphatase", "gtpase", "adaptor_scaffold", "secondary_messenger"],
                "node_type": ["protein", "protein", "protein", "protein", "molecule"],
            }
        )
        state = {
            "tq_multiplier": 1.0,
            "parameters": {
                "anchor_probability_cutoff": 0.9,
                "closure_likelihood": 0.9,
            },
        }
        audit, summary = scaffold_triadic_closure_factors(
            PROJECT_ROOT,
            symbols,
            physical_scores,
            state,
            graph_metadata=metadata,
        )
        pairs = set(zip(audit["node_a"], audit["node_b"], strict=True))
        self.assertEqual(pairs, {("A", "B"), ("A", "C"), ("B", "C")})
        self.assertTrue(np.allclose(audit["bayes_factor"], 1.8))
        self.assertTrue(np.allclose(audit["closure_support_likelihood"], 0.9))
        self.assertEqual(summary["scaffold_node_count"], 1)
        self.assertEqual(summary["anchored_protein_scaffold_associations"], 3)
        self.assertEqual(summary["degree_adjustment"], "none")
        self.assertEqual(summary["closure_bayes_factor"], 1.8)
        self.assertNotIn("cyclic AMP", set(audit["node_a"]) | set(audit["node_b"]))

    def test_scaffold_closure_can_be_enabled_without_changing_defaults(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["edge_streams"]["scaffold_triadic_closure"]["enabled"] = True
        config = normalize_configuration(supplied, self.registry)
        _, selected, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        audits: dict[str, pd.DataFrame] = {}
        matrix, summary = combine_edge_factors(
            PROJECT_ROOT,
            self.registry,
            config,
            selected["symbol"].tolist(),
            graph_metadata=selected,
            audit_collector=audits,
        )
        values = matrix.to_numpy(float)
        self.assertTrue(np.array_equal(values, values.T))
        self.assertIn("scaffold_triadic_closure", audits)
        closure = next(
            stream for stream in summary["active_streams"]
            if stream["id"] == "scaffold_triadic_closure"
        )
        self.assertTrue(closure["derived"])
        self.assertEqual(len(audits["scaffold_triadic_closure"]), 0)
        self.assertEqual(
            closure["derivation_summary"]["status"],
            "deferred_no_physical_anchor_matrix",
        )

    def test_lower_hpa_tq_increases_supported_pairs(self) -> None:
        config = normalize_configuration(None, self.registry)
        symbols = pd.read_csv(
            PROJECT_ROOT / "data/node_selection/node_universe_combined_nonzero.tsv",
            sep="\t",
            dtype=str,
        )["symbol"].drop_duplicates().tolist()
        definition = next(
            stream
            for stream in self.registry["edge_streams"]
            if stream["id"] == "hpa_primary"
        )
        reference_state = dict(config["edge_streams"]["hpa_primary"])
        reference_state["tq_multiplier"] = 1.000001
        reference_table = edge_stream_factor_table(
            PROJECT_ROOT,
            definition,
            reference_state,
            symbols,
        )
        sensitive_state = dict(config["edge_streams"]["hpa_primary"])
        self.assertEqual(sensitive_state["tq_multiplier"], 0.5)
        sensitive_table = edge_stream_factor_table(
            PROJECT_ROOT,
            definition,
            sensitive_state,
            symbols,
        )
        self.assertGreater(len(sensitive_table), len(reference_table))

    def test_more_sensitive_node_setting_adds_nodes_beyond_seed(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["node_streams"]["protein_abundance"]["tq_multiplier"] = 0.05
        config = normalize_configuration(supplied, self.registry)
        factors, selected, summary = select_nodes(PROJECT_ROOT, self.registry, config)
        # Lowering the Version 1 protein-abundance multiplier from 0.10 to the
        # validated minimum 0.05 makes that stream more sensitive.
        self.assertEqual(len(selected), 1518)
        self.assertEqual(summary["incrementally_added_protein_count"], 892)
        self.assertTrue(
            factors.loc[
                factors["gene_symbol"].isin(selected["symbol"]),
                "gui_selected_in_graph",
            ].all()
        )

    def test_incremental_pair_cache_is_complete_and_reused(self) -> None:
        supplied = default_configuration(self.registry)
        supplied["node_streams"]["protein_abundance"]["tq_multiplier"] = 0.95
        config = normalize_configuration(supplied, self.registry)
        _, selected, _ = select_nodes(PROJECT_ROOT, self.registry, config)
        seed = pd.read_csv(
            PROJECT_ROOT / "data/node_selection/node_universe_combined_nonzero.tsv",
            sep="\t",
            dtype=str,
        )["symbol"].tolist()
        ensure_incremental_pairs(PROJECT_ROOT, selected, seed)
        update = ensure_incremental_pairs(PROJECT_ROOT, selected, seed)
        expected_incremental_pairs = incremental_pair_count(
            selected["symbol"].astype(str).tolist(), set(seed)
        )
        self.assertEqual(update.requested_incremental_pairs, expected_incremental_pairs)
        self.assertEqual(update.newly_characterized_pairs, 0)
        self.assertEqual(update.cached_pairs_reused, expected_incremental_pairs)
        matrix, summary = combine_edge_factors(
            PROJECT_ROOT,
            self.registry,
            config,
            selected["symbol"].tolist(),
        )
        values = matrix.to_numpy(float)
        self.assertEqual(matrix.shape, (len(selected), len(selected)))
        self.assertTrue(np.array_equal(values, values.T))
        self.assertEqual(
            summary["unique_pair_count"],
            len(selected) * (len(selected) - 1) // 2,
        )


if __name__ == "__main__":
    unittest.main()
