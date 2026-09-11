"""Regression tests for temporal path validation."""

from __future__ import annotations

import math
import unittest

import numpy as np
import pandas as pd

from temporal_path_ranking import (
    ReplicateGeneResponse,
    SiteTrajectory,
    build_replicate_responses,
    parse_path,
    rank_paths_replicate,
    score_path_replicate,
    student_t_two_sided_p,
)


class TemporalPathRankingTests(unittest.TestCase):
    def test_parse_path_is_whitespace_tolerant(self) -> None:
        self.assertEqual(parse_path("Prkaca -> Creb1-> Aqp2"), ["Prkaca", "Creb1", "Aqp2"])

    def test_student_t_p_values_match_reference_values(self) -> None:
        self.assertTrue(math.isclose(student_t_two_sided_p(0, 6), 1.0, abs_tol=1e-14))
        self.assertTrue(
            math.isclose(
                student_t_two_sided_p(2, 6),
                0.0924263115,
                rel_tol=1e-8,
            )
        )
        self.assertTrue(
            math.isclose(
                student_t_two_sided_p(3, 10),
                0.0133436550,
                rel_tol=1e-8,
            )
        )

    def test_within_gene_bonferroni_accounts_for_site_timepoint_search(self) -> None:
        sites = [
            SiteTrajectory(
                gene="GeneA",
                uniprot="A1",
                site="10",
                times=np.array([1.0, 5.0]),
                log2_ratios=np.array([[0.8, 1.0, 1.2], [1.8, 2.0, 2.2]]),
                intensity=np.array([1e5, 1e5]),
            ),
            SiteTrajectory(
                gene="GeneA",
                uniprot="A1",
                site="20",
                times=np.array([1.0, 5.0]),
                log2_ratios=np.array([[0.1, 0.2, 0.3], [0.2, 0.3, 0.4]]),
                intensity=np.array([1e5, 1e5]),
            ),
        ]
        unadjusted, _ = build_replicate_responses(
            sites,
            prior_df=4,
            alpha=1,
            n_boot=100,
            seed=0,
            p_adjust_method="none",
        )
        adjusted, _ = build_replicate_responses(
            sites,
            prior_df=4,
            alpha=1,
            n_boot=100,
            seed=0,
            p_adjust_method="within_gene_bonferroni",
        )
        self.assertEqual(adjusted["GENEA"].gene_test_count, 4)
        self.assertTrue(
            math.isclose(
                adjusted["GENEA"].peak_p,
                min(1.0, 4 * unadjusted["GENEA"].peak_p_raw),
                rel_tol=1e-12,
            )
        )

    def test_soft_precedence_and_temporal_rank_do_not_change_primary_rank(self) -> None:
        responses = {
            "A": ReplicateGeneResponse(
                gene="A",
                site="1",
                uniprot="A1",
                passed=True,
                peak_t=5,
                peak_p_raw=0.001,
                peak_p=0.002,
                p_adjust_method="within_gene_bonferroni",
                gene_test_count=2,
                center_of_mass_point=1,
                response_time_draws=np.ones(100),
                observed_timepoint_count=2,
            ),
            "B": ReplicateGeneResponse(
                gene="B",
                site="2",
                uniprot="B1",
                passed=True,
                peak_t=5,
                peak_p_raw=0.001,
                peak_p=0.002,
                p_adjust_method="within_gene_bonferroni",
                gene_test_count=2,
                center_of_mass_point=5,
                response_time_draws=np.full(100, 5.0),
                observed_timepoint_count=2,
            ),
        }
        score = score_path_replicate(["A", "B"], responses)
        self.assertEqual(score["temporal_soft_precedence"], 1.0)
        self.assertEqual(score["temporal_kendall_tau_mean"], 1.0)
        paths = pd.DataFrame(
            {
                "rank": [1, 2],
                "path_symbols": ["B -> A", "A -> B"],
                "path_probability_product": [0.99, 0.98],
            }
        )
        annotated = rank_paths_replicate(
            paths,
            responses,
            minimum_scored_nodes=2,
        )
        self.assertEqual(annotated["rank"].tolist(), [1, 2])
        self.assertEqual(annotated["temporal_evidence_rank"].tolist(), [2, 1])


if __name__ == "__main__":
    unittest.main()
