import math
import unittest

import pandas as pd

from bayes_factors import (
    bayes_update,
    binary_bayes_factor_update,
    binary_edge_bayes_update,
    complement_minimum_factors,
    signaling_bayes_factors,
)


class SignalingBayesFactorTests(unittest.TestCase):
    def test_threshold_aware_kernel_aligns_row_specific_thresholds(self):
        values = pd.Series([1.0, 2.0, None], index=["A", "B", "C"])
        thresholds = pd.Series([2.0, 1.0, 5.0], index=["A", "B", "C"])

        result = complement_minimum_factors(values, thresholds, minimum_factor=0.0)

        self.assertAlmostEqual(result["A"], 1.0 - math.exp(-0.5 * (1.0 / 2.0) ** 2))
        self.assertAlmostEqual(result["B"], 1.0 - math.exp(-0.5 * (2.0 / 1.0) ** 2))
        self.assertEqual(result["C"], 0.0)

    def test_threshold_aware_kernel_rejects_invalid_matched_threshold(self):
        values = pd.Series([1.0], index=["A"])
        thresholds = pd.Series([0.0], index=["A"])

        with self.assertRaisesRegex(ValueError, "positive finite thresholds"):
            complement_minimum_factors(values, thresholds)

    def test_filters_output_but_uses_full_series_for_background(self):
        values = pd.Series([1.0, 2.0, 3.0, 100.0], index=["A", "B", "C", "D"])

        result = signaling_bayes_factors(
            values,
            {"A", "C"},
            q=0.5,
            minimum_factor=0.0,
        )

        self.assertEqual(result.index.tolist(), ["A", "C"])
        self.assertEqual(result.attrs["background_size"], 4)
        self.assertEqual(result.attrs["candidate_size"], 2)
        self.assertEqual(result.attrs["T_q"], 2.5)
        self.assertAlmostEqual(result["A"], 1.0 - math.exp(-0.5 * (1.0 / 2.5) ** 2))
        self.assertAlmostEqual(result["C"], 1.0 - math.exp(-0.5 * (3.0 / 2.5) ** 2))

    def test_minimum_factor_applies_to_low_and_missing_candidate_values(self):
        values = pd.Series([1.0, 10.0, None], index=["A", "B", "C"])

        result = signaling_bayes_factors(values, {"A", "C"}, q=0.5)

        self.assertEqual(result["A"], 0.5)
        self.assertEqual(result["C"], 0.5)
        self.assertEqual(result.attrs["background_size"], 2)

    def test_duplicate_gene_symbols_are_rejected(self):
        values = pd.Series([1.0, 2.0], index=["A", "A"])

        with self.assertRaisesRegex(ValueError, "duplicate gene symbols"):
            signaling_bayes_factors(values, {"A"})

    def test_invalid_quantile_is_rejected(self):
        values = pd.Series([1.0], index=["A"])

        with self.assertRaisesRegex(ValueError, "q must be between 0 and 1"):
            signaling_bayes_factors(values, {"A"}, q=75)

    def test_bayes_update_uses_floor_for_unobserved_nodes(self):
        prior = pd.Series([1 / 3, 1 / 3, 1 / 3], index=["A", "B", "C"])
        factors = pd.Series([1.0], index=["B"])

        posterior = bayes_update(prior, factors, missing_factor=0.5)

        self.assertAlmostEqual(posterior["A"], 0.25)
        self.assertAlmostEqual(posterior["B"], 0.50)
        self.assertAlmostEqual(posterior["C"], 0.25)
        self.assertAlmostEqual(posterior.sum(), 1.0)

    def test_sequential_update_order_is_commutative(self):
        prior = pd.Series([0.5, 0.5], index=["A", "B"])
        first = pd.Series([0.6, 0.9], index=["A", "B"])
        second = pd.Series([0.8, 0.5], index=["A", "B"])

        forward = bayes_update(bayes_update(prior, first), second)
        reverse = bayes_update(bayes_update(prior, second), first)

        pd.testing.assert_series_equal(forward, reverse)

    def test_binary_bayes_factor_update_preserves_neutral_prior(self):
        prior = pd.Series([0.5, 0.5], index=["A", "B"])
        factors = pd.Series([1.0, 2.0], index=prior.index)

        posterior = binary_bayes_factor_update(prior, factors)

        self.assertEqual(posterior["A"], 0.5)
        self.assertAlmostEqual(posterior["B"], 2.0 / 3.0)

    def test_binary_bayes_factor_update_does_not_normalize_across_nodes(self):
        prior = pd.Series([0.5, 0.5], index=["A", "B"])
        factors = pd.Series([2.0, 2.0], index=prior.index)

        posterior = binary_bayes_factor_update(prior, factors)

        self.assertAlmostEqual(posterior["A"], 2.0 / 3.0)
        self.assertAlmostEqual(posterior["B"], 2.0 / 3.0)
        self.assertGreater(posterior.sum(), 1.0)

    def test_binary_bayes_factor_sequential_updates_are_commutative(self):
        prior = pd.Series([0.5, 0.5], index=["A", "B"])
        first = pd.Series([1.5, 2.0], index=prior.index)
        second = pd.Series([3.0, 1.0], index=prior.index)

        forward = binary_bayes_factor_update(
            binary_bayes_factor_update(prior, first), second
        )
        reverse = binary_bayes_factor_update(
            binary_bayes_factor_update(prior, second), first
        )

        pd.testing.assert_series_equal(forward, reverse)

    def test_binary_edge_update_preserves_neutral_prior(self):
        index = ["A--B", "A--C"]
        prior = pd.Series([0.5, 0.5], index=index)
        likelihood = pd.Series([0.5, 1.0], index=index)

        posterior = binary_edge_bayes_update(prior, likelihood)

        self.assertEqual(posterior["A--B"], 0.5)
        self.assertAlmostEqual(posterior["A--C"], 2.0 / 3.0)

    def test_binary_edge_update_does_not_normalize_across_edges(self):
        prior = pd.Series([0.2, 0.2], index=["A--B", "A--C"])
        likelihood = pd.Series([1.0, 1.0], index=prior.index)

        posterior = binary_edge_bayes_update(prior, likelihood)

        self.assertAlmostEqual(posterior["A--B"], 1.0 / 3.0)
        self.assertAlmostEqual(posterior["A--C"], 1.0 / 3.0)


if __name__ == "__main__":
    unittest.main()
