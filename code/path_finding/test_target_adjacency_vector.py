"""Fast unit tests for target-vector probability helpers."""

from __future__ import annotations

import unittest

import numpy as np

from build_target_adjacency_vector import (
    apply_bayes_factor,
    complement_minimum,
)


class TargetAdjacencyProbabilityTests(unittest.TestCase):
    def test_neutral_factor_is_exact_no_op(self) -> None:
        prior = np.array([0.5, 0.7, 0.125], dtype=float)
        updated = apply_bayes_factor(prior, np.ones(3, dtype=float))
        self.assertTrue(np.array_equal(updated, prior))

    def test_factor_multiplies_odds(self) -> None:
        updated = apply_bayes_factor(
            np.array([0.5, 0.25], dtype=float),
            np.array([2.0, 3.0], dtype=float),
        )
        np.testing.assert_allclose(updated, np.array([2.0 / 3.0, 0.5]))

    def test_sequential_factors_equal_product_of_factors(self) -> None:
        prior = np.array([0.5, 0.8], dtype=float)
        first = np.array([1.5, 2.0], dtype=float)
        second = np.array([3.0, 1.25], dtype=float)
        sequential = apply_bayes_factor(
            apply_bayes_factor(prior, first),
            second,
        )
        combined = apply_bayes_factor(prior, first * second)
        np.testing.assert_allclose(sequential, combined)

    def test_complement_minimum_floor_and_threshold(self) -> None:
        values = np.array([0.0, 1.0, np.nan], dtype=float)
        thresholds = np.array([1.0, 1.0, 1.0], dtype=float)
        observed = complement_minimum(values, thresholds)
        expected_at_threshold = max(0.5, 1.0 - np.exp(-0.5))
        np.testing.assert_allclose(
            observed,
            np.array([0.5, expected_at_threshold, 0.5]),
        )

    def test_invalid_shape_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            apply_bayes_factor(np.array([0.5]), np.array([1.0, 1.0]))


if __name__ == "__main__":
    unittest.main()
