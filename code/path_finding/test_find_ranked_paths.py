"""Unit and exhaustive-order tests for ranked path finding."""

from __future__ import annotations

import math
import unittest

import numpy as np
import pandas as pd

from find_ranked_paths import (
    DEFAULT_SIGNAL_RELAY_CLASSES,
    build_adjacency,
    build_intermediate_eligibility,
    k_shortest_simple_paths,
    resolve_existing_node,
)


def brute_force_paths(
    adjacency: list[list[tuple[int, float, float]]],
    start: int,
    target: int,
    max_hops: int,
) -> list[tuple[float, tuple[int, ...]]]:
    found: list[tuple[float, tuple[int, ...]]] = []

    def visit(path: tuple[int, ...], product: float) -> None:
        node = path[-1]
        if node == target:
            found.append((product, path))
            return
        if len(path) - 1 >= max_hops:
            return
        for neighbor, probability, _ in adjacency[node]:
            if neighbor not in path:
                visit(path + (neighbor,), product * probability)

    visit((start,), 1.0)
    return sorted(found, key=lambda item: (-item[0], item[1]))


class RankedPathTests(unittest.TestCase):
    def setUp(self) -> None:
        symbols = ["A", "B", "C", "D", "T"]
        values = np.zeros((5, 5), dtype=float)
        for left, right, probability in (
            (0, 4, 0.60),
            (0, 1, 0.90),
            (1, 4, 0.90),
            (0, 2, 0.80),
            (2, 4, 0.95),
            (0, 3, 0.99),
            (3, 1, 0.99),
            (1, 2, 0.70),
            (2, 3, 0.65),
        ):
            values[left, right] = values[right, left] = probability
        self.matrix = pd.DataFrame(values, index=symbols, columns=symbols)

    def test_rank_order_matches_exhaustive_enumeration(self) -> None:
        adjacency, _ = build_adjacency(self.matrix, 0.5)
        expected = brute_force_paths(adjacency, 0, 4, max_hops=4)
        observed = k_shortest_simple_paths(
            adjacency,
            self.matrix.to_numpy(float),
            0,
            4,
            top_k=len(expected),
            max_hops=4,
        )
        self.assertEqual([item.nodes for item in observed], [item[1] for item in expected])
        np.testing.assert_allclose(
            [math.exp(-item.cost) for item in observed],
            [item[0] for item in expected],
        )

    def test_cutoff_is_strict(self) -> None:
        adjacency, edge_count = build_adjacency(self.matrix, 0.9)
        self.assertEqual(edge_count, 3)
        retained = {
            tuple(sorted((source, target)))
            for source, neighbors in enumerate(adjacency)
            for target, _, _ in neighbors
        }
        self.assertEqual(retained, {(0, 3), (2, 4), (1, 3)})

    def test_curated_pka_alias_is_explicit(self) -> None:
        metadata = pd.DataFrame(
            {
                "symbol": ["Prkaca", "Aqp2"],
                "display_symbol": ["Prkaca", "Aqp2"],
                "name": ["cAMP-dependent protein kinase catalytic subunit alpha", "Aquaporin-2"],
                "stable_id": ["MGI_SYMBOL:Prkaca", "MGI_SYMBOL:Aqp2"],
            }
        )
        resolved = resolve_existing_node(
            "PKA subunit a",
            metadata,
            ["Prkaca", "Aqp2"],
            allow_curated_alias=True,
        )
        self.assertEqual(resolved, ("Prkaca", "curated_common_name_alias"))

    def test_sensitive_signal_policy_excludes_only_scaffold_only_nodes(self) -> None:
        metadata = pd.DataFrame(
            {
                "symbol": ["Start", "Scaffold", "Mapk", "Ligand", "Binder", "Generic", "SM_CAMP", "Target"],
                "display_symbol": ["Start", "Scaffold", "Mapk", "Ligand", "Binder", "Generic", "cAMP", "Target"],
                "name": ["start", "scaffold", "kinase scaffold", "ligand", "messenger binder", "signaling participant", "cAMP", "target"],
                "classes": [
                    "receptor",
                    "adaptor_scaffold",
                    "adaptor_scaffold;kinase;signaling_process",
                    "ligand;signaling_process",
                    "second_messenger_binding",
                    "signaling_process",
                    "second_messenger",
                    "external_target",
                ],
                "node_type": ["protein", "protein", "protein", "protein", "protein", "protein", "molecule", "protein"],
            }
        )
        audit = build_intermediate_eligibility(
            metadata,
            metadata["symbol"].tolist(),
            "Start",
            "Target",
            signaling_intermediates_only=True,
            relay_classes=DEFAULT_SIGNAL_RELAY_CLASSES,
            exclude_any_scaffold=False,
        ).set_index("symbol")
        self.assertFalse(bool(audit.loc["Scaffold", "allowed_as_intermediate"]))
        self.assertTrue(bool(audit.loc["Mapk", "allowed_as_intermediate"]))
        self.assertTrue(bool(audit.loc["Ligand", "allowed_as_intermediate"]))
        self.assertTrue(bool(audit.loc["Binder", "allowed_as_intermediate"]))
        self.assertTrue(bool(audit.loc["Generic", "allowed_as_intermediate"]))
        self.assertTrue(bool(audit.loc["SM_CAMP", "allowed_as_intermediate"]))
        self.assertFalse(bool(audit.loc["Target", "allowed_as_intermediate"]))
        self.assertTrue(bool(audit.loc["Target", "permitted_in_search"]))

    def test_optional_literal_scaffold_exclusion_overrides_kinase_role(self) -> None:
        metadata = pd.DataFrame(
            {
                "symbol": ["Start", "Mapk", "Target"],
                "display_symbol": ["Start", "Mapk", "Target"],
                "name": ["start", "kinase scaffold", "target"],
                "classes": ["receptor", "adaptor_scaffold;kinase", "external_target"],
                "node_type": ["protein", "protein", "protein"],
            }
        )
        audit = build_intermediate_eligibility(
            metadata,
            metadata["symbol"].tolist(),
            "Start",
            "Target",
            signaling_intermediates_only=True,
            relay_classes=DEFAULT_SIGNAL_RELAY_CLASSES,
            exclude_any_scaffold=True,
        ).set_index("symbol")
        self.assertFalse(bool(audit.loc["Mapk", "allowed_as_intermediate"]))

    def test_ineligible_high_probability_intermediate_is_not_searched(self) -> None:
        symbols = ["Start", "Scaffold", "Kinase", "Target"]
        values = np.zeros((4, 4), dtype=float)
        for left, right, probability in (
            (0, 1, 0.99),
            (1, 3, 0.99),
            (0, 2, 0.80),
            (2, 3, 0.80),
        ):
            values[left, right] = values[right, left] = probability
        matrix = pd.DataFrame(values, index=symbols, columns=symbols)
        adjacency, _ = build_adjacency(matrix, 0.5)
        observed = k_shortest_simple_paths(
            adjacency,
            values,
            0,
            3,
            top_k=5,
            max_hops=3,
            permitted_nodes=frozenset({0, 2, 3}),
        )
        self.assertEqual([item.nodes for item in observed], [(0, 2, 3)])


if __name__ == "__main__":
    unittest.main()
