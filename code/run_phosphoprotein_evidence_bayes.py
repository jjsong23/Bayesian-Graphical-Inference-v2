"""Integrate differential-phosphosite evidence as a separate posterior branch."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from bayes_factors import binary_bayes_factor_update, complement_minimum_factors


PROJECT = Path(__file__).resolve().parents[1]
SITE_DATA = PROJECT / "results" / "kinase_predictor" / "top10_hit_analysis" / "site_data.json"
BASE_POSTERIOR = PROJECT / "results" / "node_selection_protein_pc_transcript_posterior.tsv"
UNIVERSE_TABLE = PROJECT / "data" / "node_selection" / "mouse_signaling_nodes_liberal.tsv"
OUTPUT_DIR = PROJECT / "results" / "phosphoprotein_evidence"
Q = 0.75
MINIMUM_FACTOR = 0.5


def main() -> None:
    sites = pd.DataFrame(json.loads(SITE_DATA.read_text(encoding="utf-8")))
    sites["absolute_log2_change"] = sites["raw_log2_change"].abs()
    sites["site_p_value_raw"] = pd.to_numeric(sites["site_p_value"], errors="coerce")
    sites["site_p_value_valid"] = sites["site_p_value_raw"].where(
        sites["site_p_value_raw"].between(0.0, 1.0)
    )
    sites["site_p_value_valid_flag"] = sites["site_p_value_valid"].notna()
    base = pd.read_csv(BASE_POSTERIOR, sep="\t")
    universe_info = pd.read_csv(UNIVERSE_TABLE, sep="\t", usecols=["symbol", "name", "classes"])
    universe_info = universe_info.drop_duplicates("symbol").set_index("symbol")
    universe = set(base["gene_symbol"])
    sites["in_signaling_universe"] = sites["gene_symbol"].isin(universe)
    background = sites["absolute_log2_change"].astype(float)
    site_t_q = float(background.quantile(Q))
    sites["site_T_q"] = site_t_q
    sites["site_evidence_likelihood"] = complement_minimum_factors(
        sites["absolute_log2_change"],
        site_t_q,
        minimum_factor=MINIMUM_FACTOR,
    )
    sites["site_bayes_factor"] = (
        sites["site_evidence_likelihood"] / MINIMUM_FACTOR
    )
    sites["site_positive_hit"] = sites["site_bayes_factor"] > 1.0 + 1e-12
    observed_sites = sites.loc[sites["in_signaling_universe"]].copy()

    grouped = observed_sites.groupby("gene_symbol", sort=False)
    protein_rows = []
    for gene, group in grouped:
        group = group.sort_values(
            ["absolute_log2_change", "site_p_value_valid", "site_key"],
            ascending=[False, True, True],
            na_position="last",
            kind="stable",
        )
        selected = group.iloc[0]
        n_sites = len(group)
        protein_rows.append({
            "gene_symbol": gene,
            "node_name": universe_info["name"].get(gene),
            "node_classes": universe_info["classes"].get(gene),
            "detected_phosphosites": n_sites,
            "max_absolute_lfc": float(selected["absolute_log2_change"]),
            "selected_signed_lfc": float(selected["raw_log2_change"]),
            "selected_site_key": selected["site_key"],
            "selected_uniprot": selected["uniprot"],
            "selected_site": selected["site"],
            "selected_site_p_value_raw": selected["site_p_value_raw"],
            "selected_site_p_value_valid": selected["site_p_value_valid"],
            "mean_absolute_lfc": float(group["absolute_log2_change"].mean()),
            "median_absolute_lfc": float(group["absolute_log2_change"].median()),
            "valid_site_p_values": int(group["site_p_value_valid"].notna().sum()),
            "sites_with_p_lt_0_05": int((group["site_p_value_valid"] < 0.05).sum()),
            "sites_with_positive_evidence": int(group["site_positive_hit"].sum()),
            "percent_valid_sites_with_p_lt_0_05": float(
                (group["site_p_value_valid"] < 0.05).sum()
                / group["site_p_value_valid"].notna().sum()
            ) if group["site_p_value_valid"].notna().any() else np.nan,
            "site_T_q": site_t_q,
            "selected_site_evidence_likelihood": float(
                selected["site_evidence_likelihood"]
            ),
            "selected_site_bayes_factor": float(selected["site_bayes_factor"]),
        })

    protein = pd.DataFrame(protein_rows).set_index("gene_symbol")
    factor_values = protein["max_absolute_lfc"]
    factor_thresholds = protein["site_T_q"]
    factors = complement_minimum_factors(
        factor_values,
        factor_thresholds,
        minimum_factor=MINIMUM_FACTOR,
    )
    protein["raw_evidence_factor"] = factors
    protein["relative_multiplier"] = protein["raw_evidence_factor"] / MINIMUM_FACTOR
    protein["at_neutral_floor"] = protein["relative_multiplier"] == 1.0

    base_indexed = base.set_index("gene_symbol")
    prior = base_indexed["posterior_probability"]
    relative_factors = protein["relative_multiplier"]
    updated = binary_bayes_factor_update(prior, relative_factors, missing_factor=1.0)

    integrated = base_indexed.copy()
    integrated["node_name"] = universe_info["name"].reindex(integrated.index)
    integrated["node_classes"] = universe_info["classes"].reindex(integrated.index)
    integrated["phosphosite_evidence_observed"] = integrated.index.isin(protein.index)
    for column in [
        "detected_phosphosites", "max_absolute_lfc", "selected_signed_lfc",
        "selected_site_key", "selected_uniprot", "selected_site",
        "selected_site_p_value_raw", "selected_site_p_value_valid",
        "mean_absolute_lfc", "median_absolute_lfc", "valid_site_p_values",
        "sites_with_p_lt_0_05", "sites_with_positive_evidence",
        "percent_valid_sites_with_p_lt_0_05", "site_T_q",
        "selected_site_evidence_likelihood", "selected_site_bayes_factor",
        "raw_evidence_factor", "relative_multiplier", "at_neutral_floor",
    ]:
        integrated[column] = protein[column].reindex(integrated.index)
    integrated["relative_multiplier"] = integrated["relative_multiplier"].fillna(1.0)
    integrated["posterior_before_phosphosite"] = prior
    integrated["posterior_after_phosphosite"] = updated
    integrated["posterior_change"] = integrated["posterior_after_phosphosite"] - prior
    integrated["posterior_fold_change"] = integrated["posterior_after_phosphosite"] / prior
    integrated["rank_before_phosphosite"] = prior.rank(method="min", ascending=False).astype(int)
    integrated["rank_after_phosphosite"] = updated.rank(method="min", ascending=False).astype(int)
    integrated["rank_change"] = integrated["rank_before_phosphosite"] - integrated["rank_after_phosphosite"]
    integrated = integrated.sort_values(
        ["posterior_after_phosphosite", integrated.index.name],
        ascending=[False, True],
        kind="stable",
    ).reset_index()

    protein_output = integrated.loc[integrated["phosphosite_evidence_observed"]].copy()
    selected_site_keys = set(protein["selected_site_key"])
    sites["selected_for_protein_statistic"] = sites["site_key"].isin(selected_site_keys)
    threshold = pd.DataFrame(
        [{
            "background_phosphosites": int(len(background)),
            "q": Q,
            "site_T_q": site_t_q,
            "site_count_adjustment": "none",
            "protein_rule": "maximum site Bayes factor; at least one positive site is a hit",
        }]
    )

    summary = {
        "signaling_nodes": len(integrated),
        "usable_single_phosphosites": len(sites),
        "unique_genes_with_usable_sites": int(sites["gene_symbol"].nunique()),
        "signaling_nodes_with_detected_sites": len(protein),
        "signaling_nodes_without_detected_sites": len(integrated) - len(protein),
        "observed_proteins_above_neutral": int((protein["relative_multiplier"] > 1).sum()),
        "observed_proteins_at_neutral": int((protein["relative_multiplier"] == 1).sum()),
        "background_size": len(background),
        "invalid_site_p_values": int((~sites["site_p_value_valid_flag"]).sum()),
        "sites_at_maximum_absolute_lfc": int((sites["absolute_log2_change"] == background.max()).sum()),
        "selected_proteins_at_maximum_absolute_lfc": int((protein["max_absolute_lfc"] == background.max()).sum()),
        "q": Q,
        "minimum_factor": MINIMUM_FACTOR,
        "non_detected_multiplier": 1.0,
        "protein_statistic": "maximum absolute raw phosphosite LFC per mouse gene",
        "site_threshold_method": "one empirical absolute-site-LFC q75 shared by every site; no adjustment for the number of sites on a protein",
        "protein_evidence_rule": "maximum site-level Bayes factor; at least one site with BF > 1 makes the protein a positive hit",
        "posterior_probabilities_are_independent": True,
        "posterior_minimum": float(updated.min()),
        "posterior_mean": float(updated.mean()),
        "posterior_maximum": float(updated.max()),
        "base_posterior_file": str(BASE_POSTERIOR),
        "site_data_file": str(SITE_DATA),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    integrated.to_csv(OUTPUT_DIR / "node_selection_protein_pc_phosphosite_posterior.tsv", sep="\t", index=False)
    protein_output.to_csv(OUTPUT_DIR / "phosphoprotein_evidence_nodes.tsv", sep="\t", index=False)
    threshold.to_csv(OUTPUT_DIR / "site_level_threshold.tsv", sep="\t", index=False)
    (OUTPUT_DIR / "site_count_matched_thresholds.tsv").unlink(missing_ok=True)
    sites.to_csv(OUTPUT_DIR / "phosphosite_audit.tsv", sep="\t", index=False)
    (OUTPUT_DIR / "integrated_posterior.json").write_text(integrated.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "protein_nodes.json").write_text(protein_output.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "thresholds.json").write_text(threshold.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "site_audit.json").write_text(sites.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("Top 15 integrated nodes:")
    print(integrated[["gene_symbol", "phosphosite_evidence_observed", "relative_multiplier", "posterior_after_phosphosite", "rank_change"]].head(15).to_string(index=False))
    print("Top 15 phosphoprotein evidence multipliers:")
    print(protein_output[["gene_symbol", "detected_phosphosites", "sites_with_positive_evidence", "max_absolute_lfc", "site_T_q", "relative_multiplier", "selected_site_key"]].sort_values("relative_multiplier", ascending=False).head(15).to_string(index=False))


if __name__ == "__main__":
    main()
