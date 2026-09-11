"""Combine kinase-activity and phosphoprotein evidence on the protein + PC posterior."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from bayes_factors import binary_bayes_factor_update


PROJECT = Path(__file__).resolve().parents[1]
BASE_FILE = PROJECT / "results" / "node_selection_protein_pc_transcript_posterior.tsv"
KINASE_FILE = PROJECT / "results" / "integrated_kinase_evidence" / "node_selection_protein_pc_kinase_posterior.tsv"
PHOSPHO_FILE = PROJECT / "results" / "phosphoprotein_evidence" / "node_selection_protein_pc_phosphosite_posterior.tsv"
OUTPUT_DIR = PROJECT / "results" / "combined_kinase_phosphoprotein_evidence"
USER_OUTPUT_DIR = PROJECT.parents[0] / "outputs" / "combined-kinase-phosphoprotein-evidence"
NEUTRAL_FACTOR_FLOOR = 0.5
EPSILON = 1e-12


def ranked(values: pd.Series) -> pd.Series:
    return values.rank(method="min", ascending=False).astype(int)


def main() -> None:
    base = pd.read_csv(BASE_FILE, sep="\t").set_index("gene_symbol")
    kinase = pd.read_csv(KINASE_FILE, sep="\t").set_index("gene_symbol")
    phospho = pd.read_csv(PHOSPHO_FILE, sep="\t").set_index("gene_symbol")

    if not base.index.equals(kinase.reindex(base.index).index) or not base.index.equals(phospho.reindex(base.index).index):
        raise ValueError("Evidence tables do not contain the same signaling-node universe")
    if not np.allclose(base["posterior_probability"], kinase.reindex(base.index)["posterior_before_kinase"], rtol=0, atol=1e-14):
        raise ValueError("Kinase evidence table does not start from the expected protein + PC posterior")
    if not np.allclose(base["posterior_probability"], phospho.reindex(base.index)["posterior_before_phosphosite"], rtol=0, atol=1e-14):
        raise ValueError("Phosphoprotein evidence table does not start from the expected protein + PC posterior")

    kinase = kinase.reindex(base.index)
    phospho = phospho.reindex(base.index)
    kinase_multiplier = kinase["kinase_relative_multiplier"].fillna(1.0)
    phosphoprotein_multiplier = phospho["relative_multiplier"].fillna(1.0)

    posterior_after_kinase = binary_bayes_factor_update(
        base["posterior_probability"], kinase_multiplier, missing_factor=1.0
    )
    combined_multiplier = kinase_multiplier * phosphoprotein_multiplier
    posterior_after_both = binary_bayes_factor_update(
        base["posterior_probability"], combined_multiplier, missing_factor=1.0
    )

    result = base.copy()
    result["node_name"] = kinase["node_name"].combine_first(phospho["node_name"])
    result["node_classes"] = kinase["node_classes"].combine_first(phospho["node_classes"])
    result["protein_non_neutral"] = result["protein_factor"] > NEUTRAL_FACTOR_FLOOR + EPSILON
    result["pc_transcript_non_neutral"] = result["pc_transcript_factor"] > NEUTRAL_FACTOR_FLOOR + EPSILON
    result["kinase_evidence_observed"] = kinase["kinase_evidence_observed"].fillna(False).astype(bool)
    result["kinase_non_neutral"] = kinase_multiplier > 1.0 + EPSILON
    result["kinase_relative_multiplier"] = kinase_multiplier
    result["kinase_predictor_labels"] = kinase["kinase_predictor_labels"]
    result["kinase_absolute_lfc"] = kinase["kinase_absolute_lfc"]
    result["phosphoprotein_evidence_observed"] = phospho["phosphosite_evidence_observed"].fillna(False).astype(bool)
    result["phosphoprotein_non_neutral"] = phosphoprotein_multiplier > 1.0 + EPSILON
    result["phosphoprotein_relative_multiplier"] = phosphoprotein_multiplier
    result["detected_phosphosites"] = phospho["detected_phosphosites"]
    result["maximum_absolute_site_lfc"] = phospho["max_absolute_lfc"]
    result["selected_site_key"] = phospho["selected_site_key"]
    result["combined_relative_multiplier"] = combined_multiplier
    result["posterior_after_kinase"] = posterior_after_kinase
    result["posterior_after_both_phosphoproteomic_streams"] = posterior_after_both
    result["posterior_fold_change_vs_protein_pc"] = posterior_after_both / result["posterior_probability"]
    result["rank_after_protein"] = ranked(result["posterior_after_protein"])
    result["rank_after_pc_transcript"] = ranked(result["posterior_probability"])
    result["rank_after_kinase"] = ranked(posterior_after_kinase)
    result["rank_after_both_phosphoproteomic_streams"] = ranked(posterior_after_both)
    result["rank_change_vs_protein_pc"] = result["rank_after_pc_transcript"] - result["rank_after_both_phosphoproteomic_streams"]
    result["any_non_neutral_evidence"] = result[[
        "protein_non_neutral", "pc_transcript_non_neutral", "kinase_non_neutral", "phosphoprotein_non_neutral"
    ]].any(axis=1)
    result["shared_phosphoproteomic_source_warning"] = True

    result = result.sort_values(
        ["posterior_after_both_phosphoproteomic_streams", result.index.name],
        ascending=[False, True], kind="stable"
    ).reset_index()
    non_neutral = result.loc[result["any_non_neutral_evidence"]].copy()

    kinase_supported = set(result.loc[result["kinase_non_neutral"], "gene_symbol"])
    phospho_supported = set(result.loc[result["phosphoprotein_non_neutral"], "gene_symbol"])
    summary = {
        "signaling_nodes": len(result),
        "non_neutral_components_plotted": len(non_neutral),
        "protein_non_neutral": int(result["protein_non_neutral"].sum()),
        "pc_transcript_non_neutral": int(result["pc_transcript_non_neutral"].sum()),
        "kinase_non_neutral": len(kinase_supported),
        "phosphoprotein_non_neutral": len(phospho_supported),
        "kinase_and_phosphoprotein_overlap": len(kinase_supported & phospho_supported),
        "node_prior_probability": 0.5,
        "posterior_probabilities_are_independent": True,
        "posterior_minimum": float(posterior_after_both.min()),
        "posterior_mean": float(posterior_after_both.mean()),
        "posterior_maximum": float(posterior_after_both.max()),
        "combination_rule": "kinase BF multiplied by phosphoprotein BF, then applied to each node's independent prior odds without cross-node normalization",
        "dependence_caveat": "Both multipliers derive from the same phosphoproteomic experiment; the combined result is provisional and may overstate evidence if treated as independent.",
        "non_zero_definition": "A node is plotted when at least one evidence factor is above its neutral floor (0.5 for protein/PC factors; 1 after rescaling for kinase/phosphoprotein multipliers).",
        "stage_order": ["protein abundance", "PC transcript", "kinase activity", "phosphoprotein evidence"],
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    USER_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT_DIR / "combined_all_nodes.tsv", sep="\t", index=False)
    non_neutral.to_csv(OUTPUT_DIR / "combined_non_neutral_components.tsv", sep="\t", index=False)
    result.to_csv(USER_OUTPUT_DIR / "combined_all_nodes.tsv", sep="\t", index=False)
    non_neutral.to_csv(USER_OUTPUT_DIR / "combined_non_neutral_components.tsv", sep="\t", index=False)
    (OUTPUT_DIR / "combined_all_nodes.json").write_text(result.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "combined_non_neutral_components.json").write_text(non_neutral.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(result[["gene_symbol", "posterior_after_both_phosphoproteomic_streams", "rank_after_both_phosphoproteomic_streams", "rank_change_vs_protein_pc"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
