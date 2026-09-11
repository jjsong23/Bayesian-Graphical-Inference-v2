"""Build the first undirected edge-probability matrix from localization data.

This adapts the Deshpande et al. kinase/substrate localization protocol to all
pairs in the expanded signaling-node universe:

1. Represent each detected protein by its abundance in five mpkCCD
   differential-centrifugation fractions.
2. Calculate all pairwise dot products.
3. For each target node, calculate T_q from its dot products with all other
   detected nodes (q=0.75).
4. Convert dot products to complement-of-minimum factors.
5. Average the two reciprocal directed factors to obtain an undirected
   localization likelihood.
6. Update an independent Bernoulli edge prior of 0.5, using 0.5 as the
   no-edge reference likelihood.

Pairs with missing localization data retain the 0.5 prior. Diagonal values are
set to zero by adjacency-matrix convention.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
CODE_DIR = PROJECT / "code"
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from bayes_factors import binary_edge_bayes_update  # noqa: E402


NODE_UNIVERSE = PROJECT / "data" / "node_selection" / "node_universe_combined_nonzero.tsv"
SOURCE_WORKBOOK = (
    PROJECT
    / "data"
    / "edge_characterization"
    / "localization"
    / "raw"
    / "Proteomics_of_subcellular_fractions.xlsx"
)
RESULT_DIR = PROJECT / "results" / "edge_characterization" / "localization"
DATA_DIR = PROJECT / "data" / "edge_characterization" / "localization" / "processed"

FRACTIONS = ["1K", "4K", "17K", "200Kp", "200Ks"]
Q = 0.75
MINIMUM_FACTOR = 0.5
EDGE_PRIOR = 0.5
NO_EDGE_LIKELIHOOD = 0.5


def collapse_source_profiles(source: pd.DataFrame) -> pd.DataFrame:
    """Select one nonredundant source profile per gene symbol.

    The source workbook contains repeated GI/isoform entries, frequently with
    identical five-fraction profiles. Exact duplicate profiles are collapsed;
    when distinct profiles remain, the entry with the largest total detected
    abundance is selected. This avoids artificially summing repeated isoforms.
    """
    frame = source[["GI Number", "Gene Symbol", *FRACTIONS, "Annotation"]].copy()
    frame["Gene Symbol"] = frame["Gene Symbol"].astype("string").str.strip()
    frame = frame[frame["Gene Symbol"].notna() & frame["Gene Symbol"].ne("")]
    for column in FRACTIONS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    frame["symbol_key"] = frame["Gene Symbol"].str.upper()
    frame["total_fraction_abundance"] = frame[FRACTIONS].sum(axis=1)

    rows: list[dict[str, object]] = []
    for symbol_key, group in frame.groupby("symbol_key", sort=False):
        unique_profiles = group.drop_duplicates(subset=FRACTIONS)
        selected = unique_profiles.sort_values(
            ["total_fraction_abundance", "GI Number"],
            ascending=[False, True],
            kind="stable",
        ).iloc[0]
        rows.append(
            {
                "symbol_key": symbol_key,
                "source_gene_symbol": selected["Gene Symbol"],
                "selected_gi_number": str(selected["GI Number"]),
                "annotation": selected["Annotation"],
                "source_entry_count": int(len(group)),
                "distinct_profile_count": int(len(unique_profiles)),
                **{fraction: float(selected[fraction]) for fraction in FRACTIONS},
            }
        )
    return pd.DataFrame(rows).set_index("symbol_key")


def write_matrix(path: Path, symbols: list[str], matrix: np.ndarray) -> None:
    frame = pd.DataFrame(matrix, index=symbols, columns=symbols)
    frame.index.name = "symbol"
    frame.to_csv(path, sep="\t", float_format="%.9g")


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    universe = pd.read_csv(NODE_UNIVERSE, sep="\t", dtype={"symbol": str})
    if universe["symbol"].duplicated().any():
        raise ValueError("node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    symbol_keys = pd.Index([symbol.upper() for symbol in symbols])

    source = pd.read_excel(
        SOURCE_WORKBOOK,
        sheet_name="Web Page w Links",
        header=3,
        usecols="A,C,E:I,Q",
    )
    collapsed = collapse_source_profiles(source)

    profiles = np.full((len(symbols), len(FRACTIONS)), np.nan, dtype=float)
    profile_rows: list[dict[str, object]] = []
    for index, (symbol, symbol_key) in enumerate(zip(symbols, symbol_keys, strict=True)):
        matched = collapsed.loc[symbol_key] if symbol_key in collapsed.index else None
        is_protein = universe.iloc[index]["node_type"] == "protein"
        observed = bool(
            is_protein
            and matched is not None
            and np.asarray([matched[f] for f in FRACTIONS], dtype=float).sum() > 0
        )
        if observed:
            profiles[index] = np.asarray([matched[f] for f in FRACTIONS], dtype=float)
        profile_rows.append(
            {
                "symbol": symbol,
                "display_symbol": universe.iloc[index].get("display_symbol", symbol),
                "node_type": universe.iloc[index]["node_type"],
                "localization_observed": observed,
                "selected_gi_number": matched["selected_gi_number"] if matched is not None else "",
                "source_gene_symbol": matched["source_gene_symbol"] if matched is not None else "",
                "source_entry_count": int(matched["source_entry_count"]) if matched is not None else 0,
                "distinct_profile_count": int(matched["distinct_profile_count"]) if matched is not None else 0,
                **{
                    fraction: float(matched[fraction]) if matched is not None else np.nan
                    for fraction in FRACTIONS
                },
                "annotation": matched["annotation"] if matched is not None else "",
            }
        )

    observed_mask = np.isfinite(profiles).all(axis=1)
    observed_indices = np.flatnonzero(observed_mask)
    observed_profiles = profiles[observed_mask]
    dot_observed = observed_profiles @ observed_profiles.T
    np.fill_diagonal(dot_observed, np.nan)

    thresholds = np.nanquantile(dot_observed, Q, axis=0)
    if (~np.isfinite(thresholds) | (thresholds <= 0)).any():
        bad = np.flatnonzero(~np.isfinite(thresholds) | (thresholds <= 0))
        raise ValueError(f"nonpositive or missing T_q for observed nodes: {bad[:5].tolist()}")

    z = dot_observed / thresholds[np.newaxis, :]
    directed = np.maximum(MINIMUM_FACTOR, 1.0 - np.exp(-0.5 * np.square(z)))
    np.fill_diagonal(directed, MINIMUM_FACTOR)
    symmetric_factor_observed = (directed + directed.T) / 2.0

    n = len(symbols)
    localization_factor = np.full((n, n), MINIMUM_FACTOR, dtype=float)
    localization_factor[np.ix_(observed_indices, observed_indices)] = symmetric_factor_observed
    np.fill_diagonal(localization_factor, 0.0)

    upper_i, upper_j = np.triu_indices(n, k=1)
    edge_ids = pd.Index(
        [f"{symbols[i]}--{symbols[j]}" for i, j in zip(upper_i, upper_j, strict=True)]
    )
    prior = pd.Series(EDGE_PRIOR, index=edge_ids, dtype=float)
    likelihood = pd.Series(localization_factor[upper_i, upper_j], index=edge_ids, dtype=float)
    posterior_upper = binary_edge_bayes_update(
        prior,
        likelihood,
        no_link_likelihood=NO_EDGE_LIKELIHOOD,
    ).to_numpy()
    posterior = np.zeros((n, n), dtype=float)
    posterior[upper_i, upper_j] = posterior_upper
    posterior[upper_j, upper_i] = posterior_upper

    profile_frame = pd.DataFrame(profile_rows)
    threshold_lookup = pd.Series(np.nan, index=symbols, dtype=float)
    threshold_lookup.iloc[observed_indices] = thresholds
    profile_frame["localization_Tq"] = profile_frame["symbol"].map(threshold_lookup)
    profile_frame.to_csv(DATA_DIR / "node_localization_profiles.tsv", sep="\t", index=False)

    write_matrix(RESULT_DIR / "localization_evidence_matrix.tsv", symbols, localization_factor)
    write_matrix(RESULT_DIR / "localization_adjacency_matrix.tsv", symbols, posterior)

    dot_full = np.full((n, n), np.nan, dtype=float)
    dot_full[np.ix_(observed_indices, observed_indices)] = observed_profiles @ observed_profiles.T
    observed_pair = observed_mask[upper_i] & observed_mask[upper_j]
    edge_frame = pd.DataFrame(
        {
            "node_a": np.asarray(symbols)[upper_i],
            "node_b": np.asarray(symbols)[upper_j],
            "localization_data_observed": observed_pair,
            "dot_product": dot_full[upper_i, upper_j],
            "Tq_node_a": threshold_lookup.to_numpy()[upper_i],
            "Tq_node_b": threshold_lookup.to_numpy()[upper_j],
            "localization_likelihood": localization_factor[upper_i, upper_j],
            "prior_edge_probability": EDGE_PRIOR,
            "posterior_edge_probability": posterior_upper,
        }
    )
    with gzip.open(RESULT_DIR / "localization_edges.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        edge_frame.to_csv(handle, sep="\t", index=False, float_format="%.9g")
    edge_frame.sort_values(
        ["posterior_edge_probability", "node_a", "node_b"],
        ascending=[False, True, True],
    ).head(1000).to_csv(RESULT_DIR / "top_localization_edges.tsv", sep="\t", index=False)

    observed_pairs = int(observed_pair.sum())
    possible_pairs = int(len(edge_frame))
    summary = {
        "node_count": n,
        "protein_gene_nodes": int((universe["node_type"] == "protein").sum()),
        "molecule_nodes": int((universe["node_type"] == "molecule").sum()),
        "nodes_with_localization_profiles": int(observed_mask.sum()),
        "nodes_without_localization_profiles": int((~observed_mask).sum()),
        "possible_undirected_pairs": possible_pairs,
        "pairs_with_localization_evidence": observed_pairs,
        "pairs_strengthened_by_localization": int((posterior_upper > EDGE_PRIOR + 1e-12).sum()),
        "pairs_at_neutral_prior": int((~observed_pair).sum()),
        "q": Q,
        "minimum_factor": MINIMUM_FACTOR,
        "edge_prior": EDGE_PRIOR,
        "no_edge_likelihood": NO_EDGE_LIKELIHOOD,
        "posterior_min_off_diagonal": float(posterior_upper.min()),
        "posterior_max_off_diagonal": float(posterior_upper.max()),
        "posterior_median_observed_pairs": float(np.median(posterior_upper[observed_pair])),
        "matrix_symmetric": bool(np.allclose(posterior, posterior.T)),
        "diagonal_all_zero": bool(np.allclose(np.diag(posterior), 0.0)),
        "source_workbook": str(SOURCE_WORKBOOK),
        "source_url": "https://esbl.nhlbi.nih.gov/Databases/mpkFractions/proteomic_fractions_log_files/Proteomics_of_subcellular_fractions.xlsx",
        "duplicate_resolution": "collapse exact five-fraction duplicates; select the distinct profile with greatest total abundance",
        "background_definition": (
            "For each observed target node, dot products with every other "
            f"observed node in the {n}-node universe; self excluded."
        ),
        "symmetrization": "Arithmetic mean of reciprocal directed complement-of-minimum factors.",
        "posterior_definition": "Independent Bernoulli update with P(edge)=0.5 and P(localization evidence|no edge)=0.5.",
    }
    (RESULT_DIR / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
