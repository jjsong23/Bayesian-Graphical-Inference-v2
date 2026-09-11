"""Build the active node universe from independent binary node posteriors.

At the default 0.5 prior and with positive-or-neutral evidence only, a protein
has posterior greater than 0.5 exactly when at least one evidence stream is
above its neutral floor. Curated second-messenger nodes are appended unchanged.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
COMBINED_DIR = PROJECT / "results" / "combined_kinase_phosphoprotein_evidence"
COMBINED_ALL = COMBINED_DIR / "combined_all_nodes.tsv"
COMBINED_NON_NEUTRAL = COMBINED_DIR / "combined_non_neutral_components.tsv"
UNIVERSE_INFO = (
    PROJECT / "data" / "node_selection" / "mouse_signaling_nodes_liberal.tsv"
)
SECOND_MESSENGERS = (
    PROJECT / "data" / "node_selection" / "second_messenger_nodes.tsv"
)
OUTPUT = (
    PROJECT / "data" / "node_selection" / "node_universe_combined_nonzero.tsv"
)
RESULT_DIR = PROJECT / "results" / "node_selection"
CLASS_OUTPUT = RESULT_DIR / "node_universe_class_composition.tsv"
SUMMARY_OUTPUT = RESULT_DIR / "analysis_summary.json"


NODE_COLUMNS = [
    "symbol",
    "display_symbol",
    "name",
    "classes",
    "node_type",
    "stable_id",
    "scope_tier",
    "selection_basis",
    "bayesian_score_status",
    "included_by_curated_rule",
    "source_url",
    "scope_note",
    "final_posterior",
    "final_rank",
    "protein_factor",
    "protein_non_neutral",
    "pc_transcript_factor",
    "pc_transcript_non_neutral",
    "kinase_multiplier",
    "kinase_non_neutral",
    "phosphoprotein_multiplier",
    "phosphoprotein_non_neutral",
    "combined_multiplier",
    "rank_after_protein",
    "rank_after_pc_transcript",
    "rank_after_kinase",
    "detected_phosphosites",
    "maximum_absolute_site_lfc",
    "selected_site_key",
    "kinase_predictor_labels",
]


def main() -> None:
    all_nodes = pd.read_csv(COMBINED_ALL, sep="\t")
    selected = pd.read_csv(COMBINED_NON_NEUTRAL, sep="\t")
    info = pd.read_csv(UNIVERSE_INFO, sep="\t").drop_duplicates("symbol")
    info = info.set_index("symbol")

    required = {
        "gene_symbol",
        "any_non_neutral_evidence",
        "posterior_after_both_phosphoproteomic_streams",
        "rank_after_both_phosphoproteomic_streams",
        "protein_factor",
        "protein_non_neutral",
        "pc_transcript_factor",
        "pc_transcript_non_neutral",
        "kinase_relative_multiplier",
        "kinase_non_neutral",
        "phosphoprotein_relative_multiplier",
        "phosphoprotein_non_neutral",
        "combined_relative_multiplier",
        "rank_after_protein",
        "rank_after_pc_transcript",
        "rank_after_kinase",
        "detected_phosphosites",
        "maximum_absolute_site_lfc",
        "selected_site_key",
        "kinase_predictor_labels",
    }
    missing = required - set(selected.columns)
    if missing:
        raise ValueError(f"combined node table is missing columns: {sorted(missing)}")
    if not selected["any_non_neutral_evidence"].astype(bool).all():
        raise ValueError("combined non-neutral table contains neutral-only rows")
    if selected["gene_symbol"].duplicated().any():
        raise ValueError("combined non-neutral table contains duplicate gene symbols")

    genes = pd.DataFrame(index=selected.index)
    genes["symbol"] = selected["gene_symbol"]
    genes["display_symbol"] = selected["gene_symbol"]
    genes["name"] = selected["node_name"].combine_first(
        selected["gene_symbol"].map(info["name"])
    )
    genes["classes"] = selected["node_classes"].combine_first(
        selected["gene_symbol"].map(info["classes"])
    )
    genes["node_type"] = "protein"
    genes["stable_id"] = "MGI_SYMBOL:" + selected["gene_symbol"].astype(str)
    genes["scope_tier"] = "bayesian_selected"
    genes["selection_basis"] = "binary_node_posterior_above_0.5"
    genes["bayesian_score_status"] = "scored"
    genes["included_by_curated_rule"] = False
    genes["source_url"] = ""
    genes["scope_note"] = ""
    genes["final_posterior"] = selected[
        "posterior_after_both_phosphoproteomic_streams"
    ]
    genes["final_rank"] = selected[
        "rank_after_both_phosphoproteomic_streams"
    ]
    genes["protein_factor"] = selected["protein_factor"]
    genes["protein_non_neutral"] = selected["protein_non_neutral"]
    genes["pc_transcript_factor"] = selected["pc_transcript_factor"]
    genes["pc_transcript_non_neutral"] = selected[
        "pc_transcript_non_neutral"
    ]
    genes["kinase_multiplier"] = selected["kinase_relative_multiplier"]
    genes["kinase_non_neutral"] = selected["kinase_non_neutral"]
    genes["phosphoprotein_multiplier"] = selected[
        "phosphoprotein_relative_multiplier"
    ]
    genes["phosphoprotein_non_neutral"] = selected[
        "phosphoprotein_non_neutral"
    ]
    genes["combined_multiplier"] = selected["combined_relative_multiplier"]
    genes["rank_after_protein"] = selected["rank_after_protein"]
    genes["rank_after_pc_transcript"] = selected[
        "rank_after_pc_transcript"
    ]
    genes["rank_after_kinase"] = selected["rank_after_kinase"]
    genes["detected_phosphosites"] = selected["detected_phosphosites"]
    genes["maximum_absolute_site_lfc"] = selected[
        "maximum_absolute_site_lfc"
    ]
    genes["selected_site_key"] = selected["selected_site_key"]
    genes["kinase_predictor_labels"] = selected["kinase_predictor_labels"]
    genes = genes[NODE_COLUMNS].sort_values(
        ["final_rank", "symbol"],
        ascending=[True, True],
        kind="stable",
    )

    messengers = pd.read_csv(SECOND_MESSENGERS, sep="\t")
    if messengers["symbol"].duplicated().any():
        raise ValueError("second-messenger table contains duplicate symbols")
    overlap = set(genes["symbol"]) & set(messengers["symbol"])
    if overlap:
        raise ValueError(f"gene and second-messenger symbols overlap: {sorted(overlap)}")
    for column in NODE_COLUMNS:
        if column not in messengers.columns:
            messengers[column] = pd.NA
    messengers = messengers[NODE_COLUMNS]

    node_universe = pd.concat([genes, messengers], ignore_index=True)
    if node_universe["symbol"].duplicated().any():
        raise ValueError("final node universe contains duplicate symbols")

    class_rows = []
    for node_type, subset in node_universe.groupby("node_type", dropna=False):
        memberships: dict[str, int] = {}
        for value in subset["classes"].dropna():
            for node_class in str(value).split(";"):
                node_class = node_class.strip()
                if node_class:
                    memberships[node_class] = memberships.get(node_class, 0) + 1
        for node_class, count in sorted(
            memberships.items(), key=lambda item: (-item[1], item[0])
        ):
            class_rows.append(
                {
                    "node_type": node_type,
                    "class": node_class,
                    "node_count": count,
                    "percent_of_node_type": 100.0 * count / len(subset),
                    "membership_note": (
                        "Classes are non-exclusive; one node may contribute to "
                        "multiple rows."
                    ),
                }
            )
    class_composition = pd.DataFrame(class_rows)

    summary = {
        "candidate_signaling_nodes": len(all_nodes),
        "bayesian_selected_gene_protein_nodes": len(genes),
        "curated_second_messenger_nodes": len(messengers),
        "final_node_universe_nodes": len(node_universe),
        "protein_non_neutral": int(all_nodes["protein_non_neutral"].sum()),
        "pc_transcript_non_neutral": int(
            all_nodes["pc_transcript_non_neutral"].sum()
        ),
        "kinase_non_neutral": int(all_nodes["kinase_non_neutral"].sum()),
        "phosphoprotein_non_neutral": int(
            all_nodes["phosphoprotein_non_neutral"].sum()
        ),
        "selection_rule": (
            "Retain a gene/protein when its independent present-versus-absent "
            "posterior is strictly above 0.5, then append all curated second "
            "messengers. With the current positive-or-neutral evidence, this "
            "is equivalent to at least one evidence stream being above its "
            "neutral floor."
        ),
        "protein_preprocessing": (
            "10 ** supplied Log10 Abundance before empirical-threshold scoring"
        ),
        "pc_transcript_preprocessing": (
            "finite PC median TPM > 0 only; zero or missing values neutral in "
            "integration"
        ),
        "combined_result_file": str(COMBINED_ALL),
        "selected_result_file": str(COMBINED_NON_NEUTRAL),
        "second_messenger_file": str(SECOND_MESSENGERS),
        "node_universe_file": str(OUTPUT),
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    node_universe.to_csv(OUTPUT, sep="\t", index=False)
    class_composition.to_csv(CLASS_OUTPUT, sep="\t", index=False)
    SUMMARY_OUTPUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("Top 15 selected gene/protein nodes:")
    print(
        genes[
            [
                "symbol",
                "final_posterior",
                "final_rank",
                "protein_non_neutral",
                "pc_transcript_non_neutral",
                "kinase_non_neutral",
                "phosphoprotein_non_neutral",
            ]
        ]
        .head(15)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
