"""Integrate kinase absolute-LFC evidence with the protein + PC posterior.

Kinase evidence factors are converted to multipliers relative to their 0.5
minimum. Thus a factor of 0.5 is neutral (multiplier 1), and all non-kinase
nodes also receive multiplier 1. Predictor labels are mapped to mouse genes
using the official Kinase Logos metadata table.
"""

from __future__ import annotations

import json
import ssl
import sys
from io import StringIO
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "code"))
from bayes_factors import binary_bayes_factor_update  # noqa: E402
BASE_POSTERIOR = PROJECT / "results" / "node_selection_protein_pc_transcript_posterior.tsv"
UNIVERSE_TABLE = PROJECT / "data" / "node_selection" / "mouse_signaling_nodes_liberal.tsv"
KINASE_DIR = PROJECT / "results" / "kinase_predictor" / "absolute_lfc_bayes_factors"
KINASE_RESULTS = KINASE_DIR / "kinase_absolute_lfc_bayes_factors.json"
KINASE_SUMMARY = KINASE_DIR / "analysis_summary.json"
OUTPUT_DIR = PROJECT / "results" / "integrated_kinase_evidence"
FAMILY_URL = "https://esbl.nhlbi.nih.gov/Databases/Kinase_Logos/"
EXPLICIT_MOUSE_ALIASES = {
    "CDC2": "Cdk1",
    "EGFR[T790M-L858R]": "Egfr",
    "AMPKa1": "Prkaa1",
    "AMPKa2": "Prkaa2",
    "PKACg": "Prkacg",
    "GPRK7": "Grk7",
    "ACTR2": "Actr2",
}


def clean_symbol(value) -> str | None:
    text = str(value).strip()
    if not text or text.lower() == "nan" or text in {"-", "None"}:
        return None
    return text.split("/")[0].strip()


def load_metadata() -> pd.DataFrame:
    request = Request(FAMILY_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request, context=ssl.create_default_context(), timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")
    table = max(pd.read_html(StringIO(html)), key=len).copy()
    table = table.rename(columns={
        table.columns[0]: "family",
        table.columns[1]: "predictor_label",
        table.columns[3]: "human_gene_symbol",
        table.columns[7]: "mouse_gene_symbol",
    })
    return table[["family", "predictor_label", "human_gene_symbol", "mouse_gene_symbol"]]


def main() -> None:
    base = pd.read_csv(BASE_POSTERIOR, sep="\t")
    universe_info = pd.read_csv(UNIVERSE_TABLE, sep="\t", usecols=["symbol", "name", "classes"])
    universe_info = universe_info.drop_duplicates("symbol").set_index("symbol")
    universe = set(base["gene_symbol"])
    universe_ci = {symbol.upper(): symbol for symbol in universe}

    kinase_rows = json.loads(KINASE_RESULTS.read_text(encoding="utf-8"))
    kinase_summary = json.loads(KINASE_SUMMARY.read_text(encoding="utf-8"))
    minimum_factor = float(kinase_summary["minimum_factor"])
    metadata = load_metadata()
    metadata["predictor_label"] = metadata["predictor_label"].astype(str).str.strip()

    audit_rows = []
    for row in kinase_rows:
        label = row["kinase"]
        matches = metadata.loc[metadata["predictor_label"] == label]
        mapping_method = "official mouse symbol"
        mouse_symbol = None
        human_symbol = None
        metadata_family = None
        if not matches.empty:
            for _, match in matches.iterrows():
                metadata_family = clean_symbol(match["family"])
                human_symbol = clean_symbol(match["human_gene_symbol"])
                mouse_symbol = clean_symbol(match["mouse_gene_symbol"])
                if mouse_symbol:
                    break
        if mouse_symbol and mouse_symbol not in universe:
            ci_match = universe_ci.get(mouse_symbol.upper())
            if ci_match:
                mouse_symbol = ci_match
                mapping_method = "official mouse symbol; case-normalized to universe"
        if (not mouse_symbol or mouse_symbol not in universe) and human_symbol:
            ci_match = universe_ci.get(human_symbol.upper())
            if ci_match:
                mouse_symbol = ci_match
                mapping_method = "human symbol; case-normalized to universe"
        if (not mouse_symbol or mouse_symbol not in universe):
            ci_match = universe_ci.get(label.upper())
            if ci_match:
                mouse_symbol = ci_match
                mapping_method = "predictor label; case-normalized to universe"
        if (not mouse_symbol or mouse_symbol not in universe) and label in EXPLICIT_MOUSE_ALIASES:
            mouse_symbol = EXPLICIT_MOUSE_ALIASES[label]
            mapping_method = "explicit predictor alias"

        in_universe = bool(mouse_symbol and mouse_symbol in universe)
        audit_rows.append({
            "predictor_label": label,
            "predictor_family": row["kinase_family"],
            "metadata_family": metadata_family,
            "human_gene_symbol": human_symbol,
            "mouse_gene_symbol": mouse_symbol,
            "mapping_method": mapping_method if mouse_symbol else "unmapped",
            "in_signaling_universe": in_universe,
            "absolute_lfc": row["absolute_lfc"],
            "raw_evidence_factor": row["bayes_factor"],
            "minimum_factor": minimum_factor,
            "relative_multiplier": row["bayes_factor"] / minimum_factor,
            "signed_profile_classification": row["signed_profile_classification"],
            "signed_profile_bh_fdr": row["signed_profile_bh_fdr"],
            "unique_hit_sites": row["unique_hit_sites"],
        })

    audit = pd.DataFrame(audit_rows)
    mapped = audit.loc[audit["in_signaling_universe"]].copy()
    mapped = mapped.sort_values(
        ["mouse_gene_symbol", "raw_evidence_factor", "absolute_lfc", "unique_hit_sites", "predictor_label"],
        ascending=[True, False, False, False, True],
        kind="stable",
    )
    mapped["selected_for_gene_factor"] = ~mapped["mouse_gene_symbol"].duplicated(keep="first")
    selected = mapped.loc[mapped["selected_for_gene_factor"]].set_index("mouse_gene_symbol")

    labels_by_gene = mapped.groupby("mouse_gene_symbol")["predictor_label"].agg(lambda s: "; ".join(s)).rename("kinase_predictor_labels")
    counts_by_gene = mapped.groupby("mouse_gene_symbol").size().rename("predictor_labels_mapped")
    selected = selected.join(labels_by_gene).join(counts_by_gene)
    audit = audit.merge(
        mapped[["predictor_label", "selected_for_gene_factor"]],
        on="predictor_label",
        how="left",
    )
    audit["selected_for_gene_factor"] = audit["selected_for_gene_factor"].fillna(False).astype(bool)

    integrated = base.set_index("gene_symbol").copy()
    integrated["node_name"] = universe_info["name"].reindex(integrated.index)
    integrated["node_classes"] = universe_info["classes"].reindex(integrated.index)
    integrated["kinase_evidence_observed"] = integrated.index.isin(selected.index)
    integrated["kinase_predictor_labels"] = labels_by_gene.reindex(integrated.index)
    integrated["kinase_raw_evidence_factor"] = selected["raw_evidence_factor"].reindex(integrated.index)
    integrated["kinase_minimum_factor"] = minimum_factor
    integrated["kinase_relative_multiplier"] = selected["relative_multiplier"].reindex(integrated.index).fillna(1.0)
    integrated["kinase_absolute_lfc"] = selected["absolute_lfc"].reindex(integrated.index)
    integrated["kinase_signed_profile_classification"] = selected["signed_profile_classification"].reindex(integrated.index)
    integrated["kinase_signed_profile_bh_fdr"] = selected["signed_profile_bh_fdr"].reindex(integrated.index)
    integrated["posterior_before_kinase"] = integrated["posterior_probability"]
    integrated["posterior_after_kinase"] = binary_bayes_factor_update(
        integrated["posterior_before_kinase"],
        integrated["kinase_relative_multiplier"],
        missing_factor=1.0,
    )
    integrated["posterior_change"] = integrated["posterior_after_kinase"] - integrated["posterior_before_kinase"]
    integrated["posterior_fold_change"] = integrated["posterior_after_kinase"] / integrated["posterior_before_kinase"]
    integrated["rank_before_kinase"] = integrated["posterior_before_kinase"].rank(method="min", ascending=False).astype(int)
    integrated["rank_after_kinase"] = integrated["posterior_after_kinase"].rank(method="min", ascending=False).astype(int)
    integrated["rank_change"] = integrated["rank_before_kinase"] - integrated["rank_after_kinase"]
    integrated = integrated.sort_values(["posterior_after_kinase", integrated.index.name], ascending=[False, True], kind="stable").reset_index()

    kinase_nodes = integrated.loc[integrated["kinase_evidence_observed"]].copy()
    summary = {
        "signaling_nodes": len(integrated),
        "input_kinase_predictor_labels": len(audit),
        "predictor_labels_mapped_to_universe": int(audit["in_signaling_universe"].sum()),
        "predictor_labels_unmapped_or_outside_universe": int((~audit["in_signaling_universe"]).sum()),
        "unique_signaling_nodes_with_kinase_evidence": len(selected),
        "mouse_gene_nodes_with_multiple_predictor_labels": int((counts_by_gene > 1).sum()),
        "duplicate_mapping_rule": "maximum raw evidence factor per mouse gene; ties resolved by absolute LFC, hit count, then label",
        "minimum_factor": minimum_factor,
        "non_kinase_multiplier": 1.0,
        "kinase_multiplier_formula": "raw evidence factor / minimum factor",
        "posterior_probabilities_are_independent": True,
        "posterior_minimum": float(integrated["posterior_after_kinase"].min()),
        "posterior_mean": float(integrated["posterior_after_kinase"].mean()),
        "posterior_maximum": float(integrated["posterior_after_kinase"].max()),
        "kinase_metadata_source": FAMILY_URL,
        "base_posterior_file": str(BASE_POSTERIOR),
        "kinase_factor_file": str(KINASE_RESULTS),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    integrated.to_csv(OUTPUT_DIR / "node_selection_protein_pc_kinase_posterior.tsv", sep="\t", index=False)
    kinase_nodes.to_csv(OUTPUT_DIR / "kinase_nodes_integrated_posterior.tsv", sep="\t", index=False)
    audit.to_csv(OUTPUT_DIR / "kinase_predictor_to_mouse_mapping.tsv", sep="\t", index=False)
    (OUTPUT_DIR / "integrated_posterior.json").write_text(integrated.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "kinase_nodes.json").write_text(kinase_nodes.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "mapping_audit.json").write_text(audit.to_json(orient="records", double_precision=15), encoding="utf-8")
    (OUTPUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("Top 15 integrated nodes:")
    print(integrated[["gene_symbol", "kinase_evidence_observed", "kinase_relative_multiplier", "posterior_after_kinase", "rank_change"]].head(15).to_string(index=False))
    print("Unmapped/outside predictor labels:")
    print(audit.loc[~audit["in_signaling_universe"], ["predictor_label", "mouse_gene_symbol", "human_gene_symbol"]].to_string(index=False))


if __name__ == "__main__":
    main()
