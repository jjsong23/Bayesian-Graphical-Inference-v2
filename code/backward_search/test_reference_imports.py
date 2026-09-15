from __future__ import annotations

import shutil
import unittest
import uuid
from pathlib import Path

import pandas as pd

from .engine import PairCache
from .huri import build_huri_search_evidence
from .legacy_screen import import_legacy_ppi_screen


class ReferenceImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.root = Path.cwd() / ".test_work" / uuid.uuid4().hex
        self.root.mkdir(parents=True)

    def tearDown(self) -> None:
        shutil.rmtree(self.root, ignore_errors=True)

    def test_legacy_screen_requires_valid_score_and_structure_and_keeps_low_scores(self) -> None:
        screen = self.root / "ppi_screen"
        (screen / "scores").mkdir(parents=True)
        (screen / "hits/A_and_B").mkdir(parents=True)
        (screen / "hits/A_and_B/ranked_0.pdb").write_text("ATOM\n", encoding="utf-8")
        pd.DataFrame([
            {"pair": "A_and_B", "iptm": "0.12", "ranking_confidence": "0.2"},
            {"pair": "A_and_C", "iptm": "0.8", "ranking_confidence": "0.9"},
            {"pair": "A_and_D", "iptm": "bad", "ranking_confidence": "0.9"},
        ]).to_csv(screen / "scores/scores_1.csv", index=False)
        cache_path = self.root / "cache.sqlite3"
        summary = import_legacy_ppi_screen(screen, cache_path)
        self.assertEqual(summary["valid_unique_pair_count"], 1)
        self.assertTrue(summary["low_scores_retained"])
        cache = PairCache(cache_path)
        try:
            self.assertAlmostEqual(
                cache.structural_score(("A", "B"), summary["protocol_id"], "iptm"),
                0.12,
            )
            self.assertIsNone(cache.structural_score(("A", "C"), summary["protocol_id"], "iptm"))
        finally:
            cache.close()

    def test_huri_positive_and_scoped_nonreport_are_separate(self) -> None:
        interactions = self.root / "HuRI.tsv"
        interactions.write_text(
            "ENSG000001\tENSG000002\nENSG000002\tENSG000003\n",
            encoding="utf-8",
        )
        mapping = pd.DataFrame({
            "symbol": ["A", "B", "C", "D"],
            "hpa_ensembl_gene_id": ["ENSG000001", "ENSG000002", "ENSG000003", "ENSG999999"],
        })
        mapping.to_csv(self.root / "mapping.tsv", sep="\t", index=False)
        nodes = pd.DataFrame({"symbol": ["A", "B", "C", "D"], "node_type": ["protein"] * 4})
        stream, summary = build_huri_search_evidence(
            interactions,
            self.root / "mapping.tsv",
            nodes,
            self.root / "huri_output",
            positive_bayes_factor=5.0,
            nonreported_bayes_factor=0.9,
            nonreported_scope="huri_interactor_set",
        )
        factors = pd.read_csv(stream["file"], sep="\t")
        self.assertEqual(len(factors), 2)
        self.assertEqual(summary["nonreported_scope_node_count"], 3)
        scope = pd.read_csv(stream["scope_nodes_file"], sep="\t")
        self.assertEqual(set(scope["symbol"]), {"A", "B", "C"})
        self.assertEqual(stream["scoped_missing_bayes_factor"], 0.9)


if __name__ == "__main__":
    unittest.main()
