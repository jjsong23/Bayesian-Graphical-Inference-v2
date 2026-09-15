from __future__ import annotations

import json
import shutil
import unittest
import uuid
from pathlib import Path

import pandas as pd

from .engine import (
    PairCache,
    build_node_sequence_fasta,
    cached_structural_result,
    cleanup_run_features,
    collect_alphapulldown_scores,
    initialize_run,
    import_structural_scores,
    load_run_state,
    step_run,
    structural_bayes_factor,
)
from .development import record_structural_import


class BackwardSearchTests(unittest.TestCase):
    def setUp(self) -> None:
        self.project = Path.cwd() / ".test_work" / uuid.uuid4().hex
        self.project.mkdir(parents=True)
        nodes = pd.DataFrame(
            {
                "symbol": ["END", "REC", "A", "B", "C"],
                "node_type": ["protein"] * 5,
                "classes": ["target", "receptor", "kinase", "kinase", "kinase"],
            }
        )
        nodes.to_csv(self.project / "nodes.tsv", sep="\t", index=False)
        with (self.project / "sequences.fasta").open("w", encoding="utf-8") as handle:
            for symbol in nodes["symbol"]:
                handle.write(f">{symbol}\nMSTNPKPQRKTKRNTNRRPQ\n")
        pd.DataFrame(
            [
                ("END", "A", 9.0),
                ("END", "B", 4.0),
                ("END", "C", 0.5),
                ("END", "REC", 1.0),
            ],
            columns=["node_a", "node_b", "bayes_factor"],
        ).to_csv(self.project / "cheap.tsv", sep="\t", index=False)
        configuration = {
            "project_root": ".",
            "target": "END",
            "receptor": "REC",
            "node_table": "nodes.tsv",
            "sequence_fasta": "sequences.fasta",
            "allowed_intermediate_classes": ["kinase", "receptor"],
            "cheap_top_n_per_frontier": 2,
            "beam_width": 2,
            "maximum_depth": 4,
            "cheap_streams": [
                {
                    "id": "cheap",
                    "file": "cheap.tsv",
                    "weight": 1.0,
                    "missing_bayes_factor": 1.0,
                }
            ],
            "structural": {
                "enabled": True,
                "metric": "iptm",
                "reference_score": 0.6,
                "weight": 1.0,
                "bayes_factor_floor": 0.1,
                "bayes_factor_ceiling": 10.0,
                "protocol_id": "test_protocol",
            },
        }
        self.config_path = self.project / "config.json"
        self.config_path.write_text(json.dumps(configuration), encoding="utf-8")
        self.run_dir = self.project / "run"

    def tearDown(self) -> None:
        shutil.rmtree(self.project, ignore_errors=True)

    def _import_scores(self, rows: list[tuple[str, str, float]], name: str) -> None:
        path = self.project / name
        pd.DataFrame(rows, columns=["node_a", "node_b", "score"]).to_csv(
            path, sep="\t", index=False
        )
        import_structural_scores(self.run_dir, path)

    def test_resumable_two_round_search_and_pair_cache(self) -> None:
        state = initialize_run(self.config_path, self.run_dir)
        self.assertEqual(state["status"], "ready")

        state = step_run(self.run_dir)
        self.assertEqual(state["status"], "waiting_for_structural")
        self.assertEqual(state["pending_unique_structural_pairs"], 3)
        pairs = pd.read_csv(self.run_dir / "round_000/pairs.tsv", sep="\t")
        self.assertEqual(set(pairs["node_b"]) | set(pairs["node_a"]), {"END", "A", "B", "REC"})

        self._import_scores(
            [("END", "A", 0.55), ("END", "B", 0.90), ("END", "REC", 0.10)],
            "round0_scores.tsv",
        )
        state = step_run(self.run_dir)
        self.assertEqual(state["status"], "ready")
        self.assertEqual(state["round"], 1)
        self.assertEqual([path["nodes_backward"][-1] for path in state["beam"]], ["B", "A"])

        state = step_run(self.run_dir)
        self.assertEqual(state["status"], "waiting_for_structural")
        pending = pd.read_csv(self.run_dir / "round_001/pairs.tsv", sep="\t")
        score_rows = []
        for row in pending.to_dict("records"):
            pair = {row["node_a"], row["node_b"]}
            score = 0.95 if pair == {"B", "REC"} else 0.10
            score_rows.append((row["node_a"], row["node_b"], score))
        self._import_scores(score_rows, "round1_scores.tsv")
        state = step_run(self.run_dir)
        self.assertEqual(state["status"], "complete")
        self.assertEqual(state["completion_reason"], "receptor_reached")
        self.assertEqual(state["solutions"][0]["nodes_forward"], ["REC", "B", "END"])

        second_run = self.project / "second_run"
        initialize_run(self.config_path, second_run)
        second_state = step_run(second_run)
        self.assertNotEqual(second_state["status"], "waiting_for_structural")

    def test_development_mode_pauses_after_each_auditable_decision(self) -> None:
        state = initialize_run(self.config_path, self.run_dir, development_mode=True)
        self.assertEqual(state["status"], "development_paused")
        self.assertEqual(state["development_stage"], "enumerate_candidates")
        self.assertTrue((self.run_dir / "development_trace.html").is_file())

        expected_stages = [
            "score_cheap_evidence",
            "select_cheap_shortlist",
            "prepare_structural",
            "integrate_edge_evidence",
        ]
        for expected in expected_stages:
            state = step_run(self.run_dir)
            self.assertEqual(state["development_stage"], expected)

        self.assertEqual(state["status"], "waiting_for_structural")
        round_dir = self.run_dir / "round_000"
        for name in (
            "candidate_enumeration.tsv",
            "cheap_evidence_ledger.tsv",
            "cheap_scored_all.tsv",
            "cheap_selection.tsv",
            "structural_cache_lookup.tsv",
            "pairs.tsv",
        ):
            self.assertTrue((round_dir / name).is_file(), name)
        enumeration = pd.read_csv(round_dir / "candidate_enumeration.tsv", sep="\t")
        self.assertIn("ontology_or_eligibility_filter", set(enumeration["exclusion_reason"].dropna()))
        selection = pd.read_csv(round_dir / "cheap_selection.tsv", sep="\t")
        self.assertEqual(int(selection["selected_for_structural_stage"].sum()), 3)

        self._import_scores(
            [("END", "A", 0.55), ("END", "B", 0.90), ("END", "REC", 0.10)],
            "development_round0_scores.tsv",
        )
        state = record_structural_import(self.run_dir, load_run_state(self.run_dir))
        self.assertEqual(state["pending_unique_structural_pairs"], 0)
        self.assertEqual(state["status"], "development_paused")
        state = step_run(self.run_dir)
        self.assertEqual(state["development_stage"], "select_frontier")
        self.assertTrue((round_dir / "edge_evidence_ledger.tsv").is_file())
        state = step_run(self.run_dir)
        self.assertEqual(state["round"], 1)
        self.assertEqual(state["development_stage"], "enumerate_candidates")
        self.assertEqual([path["nodes_backward"][-1] for path in state["beam"]], ["B", "A"])
        events = (self.run_dir / "development_events.jsonl").read_text(encoding="utf-8")
        self.assertIn('"action": "select_frontier"', events)
        self.assertIn("edge_evidence_ledger.tsv", (self.run_dir / "development_trace.html").read_text(encoding="utf-8"))

    def test_explicit_disallowed_backward_traversal_is_removed(self) -> None:
        direction = pd.DataFrame(
            [("A", "END", False)], columns=["source", "target", "allowed"]
        )
        direction.to_csv(self.project / "direction.tsv", sep="\t", index=False)
        configuration = json.loads(self.config_path.read_text(encoding="utf-8"))
        configuration["directionality_file"] = "direction.tsv"
        self.config_path.write_text(json.dumps(configuration), encoding="utf-8")
        initialize_run(self.config_path, self.run_dir)
        step_run(self.run_dir)
        candidates = pd.read_csv(self.run_dir / "round_000/cheap_candidates.tsv", sep="\t")
        self.assertNotIn("A", set(candidates["candidate_parent"]))

    def test_structural_reference_is_neutral_and_bounds_are_respected(self) -> None:
        self.assertAlmostEqual(structural_bayes_factor(0.6, 0.6, 0.1, 10.0), 1.0)
        self.assertEqual(structural_bayes_factor(0.0, 0.6, 0.1, 10.0), 0.1)
        self.assertEqual(structural_bayes_factor(1.0, 0.6, 0.1, 10.0), 10.0)

    def test_state_is_readable_after_atomic_writes(self) -> None:
        initialize_run(self.config_path, self.run_dir)
        self.assertEqual(load_run_state(self.run_dir)["schema_version"], 1)

    def test_round_package_requests_exactly_one_model_prediction(self) -> None:
        initialize_run(self.config_path, self.run_dir)
        step_run(self.run_dir)
        config = json.loads((self.run_dir / "round_000/round_config.json").read_text())
        self.assertEqual(config["num_predictions_per_model"], 1)
        self.assertEqual(config["model_names"], "model_1_multimer_v3")

    def test_legacy_protocol_is_reused_only_when_explicitly_compatible(self) -> None:
        cache = PairCache(self.project / "compatibility.sqlite3")
        try:
            frame = pd.DataFrame([{"node_a": "A", "node_b": "B", "score": 0.72}])
            cache.import_structural(frame, "legacy", "iptm", self.project / "legacy.tsv")
            structural = {"protocol_id": "new", "compatible_protocol_ids": ["legacy"], "metric": "iptm"}
            result = cached_structural_result(cache, ("B", "A"), structural)
            self.assertIsNotNone(result)
            self.assertEqual(result["protocol_id"], "legacy")
            self.assertAlmostEqual(result["score"], 0.72)
            structural["compatible_protocol_ids"] = []
            self.assertIsNone(cached_structural_result(cache, ("A", "B"), structural))
        finally:
            cache.close()

    def test_completed_frontier_cleanup_keeps_only_next_frontier_features(self) -> None:
        round_dir = self.run_dir / "round_000"
        feature_dir = self.run_dir / "features"
        round_dir.mkdir(parents=True)
        feature_dir.mkdir(parents=True)
        (round_dir / "sequences.fasta").write_text(">END\nMAAA\n>A\nMBBB\n>B\nMCCC\n", encoding="utf-8")
        for symbol in ("END", "A", "B"):
            (feature_dir / f"{symbol}.pkl").write_bytes(b"feature-object")
            msa = feature_dir / symbol
            msa.mkdir()
            (msa / "uniref90_hits.sto").write_bytes(b"msa")
        summary = cleanup_run_features(
            self.run_dir, {"A"}, 0, reason="frontier_advanced"
        )
        self.assertTrue((feature_dir / "A.pkl").is_file())
        self.assertTrue((feature_dir / "A").is_dir())
        self.assertFalse((feature_dir / "END.pkl").exists())
        self.assertFalse((feature_dir / "B").exists())
        self.assertEqual(summary["deleted_artifact_count"], 4)
        self.assertTrue((self.run_dir / "feature_cleanup.tsv").is_file())

    def test_collects_one_best_iptm_score_for_every_manifest_pair(self) -> None:
        initialize_run(self.config_path, self.run_dir)
        step_run(self.run_dir)
        round_dir = self.run_dir / "round_000"
        pairs = pd.read_csv(round_dir / "pairs.tsv", sep="\t")
        for number, pair_id in enumerate(pairs["pair_id"], start=1):
            output = round_dir / "models" / pair_id
            output.mkdir(parents=True)
            pd.DataFrame({"iptm": [0.1, 0.2 + number / 100.0]}).to_csv(
                output / "predictions_with_good_interpae.csv", index=False
            )
        score_path = collect_alphapulldown_scores(self.run_dir)
        scores = pd.read_csv(score_path, sep="\t")
        self.assertEqual(len(scores), len(pairs))
        self.assertTrue((scores["metric"] == "iptm").all())

    def test_partial_import_does_not_rewrite_pair_manifest(self) -> None:
        initialize_run(self.config_path, self.run_dir)
        step_run(self.run_dir)
        manifest = self.run_dir / "round_000/pairs.tsv"
        original = manifest.read_text(encoding="utf-8")
        pairs = pd.read_csv(manifest, sep="\t")
        first = pairs.iloc[0]
        self._import_scores(
            [(first["node_a"], first["node_b"], 0.8)],
            "partial_scores.tsv",
        )
        state = step_run(self.run_dir)
        self.assertEqual(state["status"], "waiting_for_structural")
        self.assertEqual(state["pending_unique_structural_pairs"], len(pairs) - 1)
        self.assertEqual(manifest.read_text(encoding="utf-8"), original)

    def test_build_fasta_prefers_explicit_accession_and_audits_mapping(self) -> None:
        nodes = pd.DataFrame(
            {"symbol": ["A", "B"], "selected_uniprot": ["P00001", ""]}
        )
        nodes.to_csv(self.project / "fasta_nodes.tsv", sep="\t", index=False)
        uniprot = pd.DataFrame(
            {
                "Entry": ["P00001", "P00002"],
                "Sequence": ["MAAA", "MBBB"],
                "Gene Names (primary)": ["WrongSymbol", "B"],
                "Gene Names (synonym)": ["A", ""],
                "Reviewed": ["reviewed", "reviewed"],
            }
        )
        uniprot.to_csv(self.project / "uniprot.tsv.gz", sep="\t", index=False, compression="gzip")
        summary = build_node_sequence_fasta(
            self.project / "fasta_nodes.tsv",
            self.project / "uniprot.tsv.gz",
            self.project / "built.fasta",
        )
        self.assertEqual(summary["sequence_count"], 2)
        fasta = (self.project / "built.fasta").read_text(encoding="utf-8")
        self.assertIn(">A\nMAAA", fasta)
        self.assertIn(">B\nMBBB", fasta)
        mapping = pd.read_csv(self.project / "built.mapping.tsv", sep="\t")
        self.assertEqual(mapping.loc[mapping["symbol"] == "A", "mapping_method"].iloc[0], "selected_uniprot_accession")


if __name__ == "__main__":
    unittest.main()
