#!/usr/bin/env python3
"""Sensitivity sweep for replicate-aware temporal path validation.

The sweep varies the variance-moderation prior degrees of freedom and the
gene-level significance gate.  It reports which gene calls and path-order
conclusions are robust to those analysis choices.  The largest default prior
(``prior_df=50``) is deliberately an over-shrinkage stress test, not a preferred
setting.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import pandas as pd

import temporal_path_ranking as temporal


def _floats(value: str) -> list[float]:
    return [float(item) for item in str(value).split(",") if item.strip()]


def path_genes(paths: pd.DataFrame, path_col: str, delimiter: str) -> list[str]:
    genes: set[str] = set()
    for path in paths[path_col]:
        genes.update(symbol.upper() for symbol in temporal.parse_path(path, delimiter))
    return sorted(genes)


def build_by_prior_df(
    sites,
    prior_dfs: list[float],
    n_boot: int,
    seed: int,
    p_adjust_method: str,
):
    """Build once per prior; alpha is applied later without redrawing times."""

    return {
        prior_df: temporal.build_replicate_responses(
            sites,
            prior_df=prior_df,
            alpha=1.0,
            n_boot=n_boot,
            seed=seed,
            p_adjust_method=p_adjust_method,
        )[0]
        for prior_df in prior_dfs
    }


def regate(responses: dict, alpha: float) -> dict:
    return {
        gene: replace(response, passed=response.peak_p < alpha)
        for gene, response in responses.items()
    }


def _cell(
    paths: pd.DataFrame,
    responses: dict,
    rank: int,
    path_col: str,
    delimiter: str,
) -> dict[str, object]:
    scored = temporal.rank_paths_replicate(
        paths,
        responses,
        path_col=path_col,
        delimiter=delimiter,
        minimum_scored_nodes=2,
    )
    row = scored.loc[scored["rank"].eq(rank)].iloc[0]
    return {
        "n_scored": int(row.temporal_n_scored),
        "kendall_tau_mean": row.temporal_kendall_tau_mean,
        "kendall_tau_low": row.temporal_kendall_tau_low,
        "kendall_tau_high": row.temporal_kendall_tau_high,
        "soft_precedence": row.temporal_soft_precedence,
    }


def _format_cell(cell: dict[str, object]) -> str:
    if pd.isna(cell["kendall_tau_mean"]):
        return f"{cell['n_scored']}|NA"
    return (
        f"{cell['n_scored']}|{cell['kendall_tau_mean']:+.2f}|"
        f"[{cell['kendall_tau_low']:+.2f},{cell['kendall_tau_high']:+.2f}]"
    )


def gene_p_table(builds: dict, prior_dfs: list[float], genes: list[str]) -> pd.DataFrame:
    present = [gene for gene in genes if gene in builds[prior_dfs[0]]]
    return pd.DataFrame(
        [
            [gene]
            + [round(builds[prior_df][gene].peak_p, 6) for prior_df in prior_dfs]
            for gene in present
        ],
        columns=["gene"] + [f"prior_df={prior_df:g}" for prior_df in prior_dfs],
    )


def select_target_ranks(
    paths: pd.DataFrame,
    builds: dict,
    prior_dfs: list[float],
    alphas: list[float],
    minimum_scored: int,
    path_col: str,
    delimiter: str,
) -> list[int]:
    maximum_scored = {int(rank): 0 for rank in paths["rank"]}
    for prior_df in prior_dfs:
        for alpha in alphas:
            scored = temporal.rank_paths_replicate(
                paths,
                regate(builds[prior_df], alpha),
                path_col=path_col,
                delimiter=delimiter,
                minimum_scored_nodes=minimum_scored,
            )
            for row in scored.itertuples():
                maximum_scored[int(row.rank)] = max(
                    maximum_scored[int(row.rank)], int(row.temporal_n_scored)
                )
    return sorted(
        rank for rank, scored in maximum_scored.items() if scored >= minimum_scored
    )


def path_table_vary_prior(
    paths: pd.DataFrame,
    builds: dict,
    prior_dfs: list[float],
    fixed_alpha: float,
    ranks: list[int],
    path_col: str,
    delimiter: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            [rank]
            + [
                _format_cell(
                    _cell(
                        paths,
                        regate(builds[prior_df], fixed_alpha),
                        rank,
                        path_col,
                        delimiter,
                    )
                )
                for prior_df in prior_dfs
            ]
            for rank in ranks
        ],
        columns=["rank"] + [f"prior_df={value:g}" for value in prior_dfs],
    )


def path_table_vary_alpha(
    paths: pd.DataFrame,
    builds: dict,
    fixed_prior_df: float,
    alphas: list[float],
    ranks: list[int],
    path_col: str,
    delimiter: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            [rank]
            + [
                _format_cell(
                    _cell(
                        paths,
                        regate(builds[fixed_prior_df], alpha),
                        rank,
                        path_col,
                        delimiter,
                    )
                )
                for alpha in alphas
            ]
            for rank in ranks
        ],
        columns=["rank"] + [f"alpha={value:g}" for value in alphas],
    )


def gate_stability(
    builds: dict,
    prior_dfs: list[float],
    alphas: list[float],
    genes: list[str],
) -> tuple[list[str], list[str], list[str]]:
    present = [gene for gene in genes if gene in builds[prior_dfs[0]]]
    always_pass: list[str] = []
    always_fail: list[str] = []
    threshold_dependent: list[str] = []
    for gene in present:
        calls = {
            builds[prior_df][gene].peak_p < alpha
            for prior_df in prior_dfs
            for alpha in alphas
        }
        if calls == {True}:
            always_pass.append(gene)
        elif calls == {False}:
            always_fail.append(gene)
        else:
            threshold_dependent.append(gene)
    return always_pass, always_fail, threshold_dependent


def cells_long(
    paths: pd.DataFrame,
    builds: dict,
    prior_dfs: list[float],
    alphas: list[float],
    ranks: list[int],
    path_col: str,
    delimiter: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for prior_df in prior_dfs:
        for alpha in alphas:
            responses = regate(builds[prior_df], alpha)
            for rank in ranks:
                rows.append(
                    {
                        "rank": rank,
                        "prior_df": prior_df,
                        "alpha": alpha,
                        **_cell(paths, responses, rank, path_col, delimiter),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-phospho", required=True)
    parser.add_argument("--paths", required=True)
    parser.add_argument("--path-col", default=temporal.DEFAULT_PATH_COL)
    parser.add_argument("--delim", default=temporal.DEFAULT_PATH_DELIM)
    parser.add_argument("--prior-dfs", type=_floats, default=[1, 2, 4, 8, 50])
    parser.add_argument("--alphas", type=_floats, default=[0.01, 0.05, 0.1, 0.2])
    parser.add_argument("--fixed-alpha", type=float, default=0.05)
    parser.add_argument("--fixed-prior-df", type=float, default=4.0)
    parser.add_argument("--min-scored", type=int, default=3)
    parser.add_argument(
        "--target-ranks",
        type=lambda value: [int(item) for item in value.split(",") if item.strip()],
    )
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--p-adjust-method",
        choices=temporal.P_ADJUST_METHODS,
        default="within_gene_bonferroni",
    )
    parser.add_argument("--out")
    parser.add_argument("--cells-tsv")
    arguments = parser.parse_args()

    prior_dfs = sorted(set(arguments.prior_dfs) | {arguments.fixed_prior_df})
    alphas = sorted(set(arguments.alphas) | {arguments.fixed_alpha})
    sites = temporal.load_site_trajectories(arguments.raw_phospho)
    paths = pd.read_csv(arguments.paths, sep="\t")
    genes = path_genes(paths, arguments.path_col, arguments.delim)
    builds = build_by_prior_df(
        sites,
        prior_dfs,
        arguments.n_boot,
        arguments.seed,
        arguments.p_adjust_method,
    )
    ranks = arguments.target_ranks or select_target_ranks(
        paths,
        builds,
        prior_dfs,
        alphas,
        arguments.min_scored,
        arguments.path_col,
        arguments.delim,
    )
    table_a = gene_p_table(builds, prior_dfs, genes)
    table_b1 = path_table_vary_prior(
        paths,
        builds,
        prior_dfs,
        arguments.fixed_alpha,
        ranks,
        arguments.path_col,
        arguments.delim,
    )
    table_b2 = path_table_vary_alpha(
        paths,
        builds,
        arguments.fixed_prior_df,
        alphas,
        ranks,
        arguments.path_col,
        arguments.delim,
    )
    always_pass, always_fail, flips = gate_stability(builds, prior_dfs, alphas, genes)
    header = (
        f"replicate-mode sensitivity: {len(sites)} sites, prior_dfs={prior_dfs}, "
        f"alphas={alphas}, draws={arguments.n_boot}, seed={arguments.seed}, "
        f"p_adjust={arguments.p_adjust_method}\n"
        f"target paths (maximum n_scored >= {arguments.min_scored}): {ranks}\n"
        "cells: n_scored|kendall_tau_mean|[95% Monte-Carlo interval]"
    )
    print(header)
    print("\nA. adjusted peak p by prior_df")
    print(table_a.to_string(index=False))
    print(f"\nB1. alpha={arguments.fixed_alpha:g}; vary prior_df")
    print(table_b1.to_string(index=False))
    print(f"\nB2. prior_df={arguments.fixed_prior_df:g}; vary alpha")
    print(table_b2.to_string(index=False))
    print("\nC. gate stability")
    print(f"always pass: {always_pass}")
    print(f"always fail: {always_fail}")
    print(f"threshold-dependent: {flips}")

    if arguments.cells_tsv:
        cells_long(
            paths,
            builds,
            prior_dfs,
            alphas,
            ranks,
            arguments.path_col,
            arguments.delim,
        ).to_csv(arguments.cells_tsv, sep="\t", index=False)
    if arguments.out:
        report = [
            "# Temporal-validation sensitivity\n",
            "```\n" + header + "\n```\n",
            "## A. Adjusted peak p by prior df\n",
            "```\n" + table_a.to_string(index=False) + "\n```\n",
            f"## B1. Alpha={arguments.fixed_alpha:g}; vary prior df\n",
            "```\n" + table_b1.to_string(index=False) + "\n```\n",
            f"## B2. Prior df={arguments.fixed_prior_df:g}; vary alpha\n",
            "```\n" + table_b2.to_string(index=False) + "\n```\n",
            "## C. Gate stability\n",
            f"- Always pass: {always_pass}\n",
            f"- Always fail: {always_fail}\n",
            f"- Threshold-dependent: {flips}\n",
        ]
        Path(arguments.out).write_text("\n".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
