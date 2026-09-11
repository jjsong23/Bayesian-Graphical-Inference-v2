"""Integrate principal-cell scRNA-seq with mpkCCD protein-abundance evidence."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from bayes_factors import binary_bayes_factor_update, signaling_bayes_factors


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = PROJECT_ROOT / "data" / "pc_scrna_seq_raw.xlsx"
UNIVERSE = (
    PROJECT_ROOT
    / "data"
    / "node_selection"
    / "mouse_signaling_symbols_liberal.txt"
)
PROTEIN_RESULTS = (
    PROJECT_ROOT / "results" / "mpkccd_protein_abundance_bayes_factors.tsv"
)
TRANSCRIPT_RESULTS = (
    PROJECT_ROOT / "results" / "pc_median_tpm_bayes_factors.tsv"
)
COMBINED_RESULTS = (
    PROJECT_ROOT / "results" / "node_selection_protein_pc_transcript_posterior.tsv"
)
SUMMARY_OUTPUT = (
    PROJECT_ROOT / "results" / "node_selection_protein_pc_transcript_summary.json"
)


def main() -> None:
    universe = UNIVERSE.read_text(encoding="utf-8").splitlines()
    universe_index = pd.Index(universe, name="gene_symbol")

    transcript_table = pd.read_excel(
        WORKBOOK,
        sheet_name="Median TPM",
        usecols=["Gene Symbol", "PC (Median TPM,n=74)"],
    )
    pc_tpm_all = pd.to_numeric(
        transcript_table.set_index("Gene Symbol")["PC (Median TPM,n=74)"],
        errors="coerce",
    )
    finite_mask = np.isfinite(pc_tpm_all.to_numpy(dtype=float))
    positive_mask = finite_mask & (pc_tpm_all.to_numpy(dtype=float) > 0.0)
    pc_tpm_positive = pc_tpm_all.iloc[np.flatnonzero(positive_mask)].copy()
    transcript_factors = signaling_bayes_factors(
        pc_tpm_positive,
        universe,
        q=0.75,
    )

    transcript_output = pd.DataFrame(
        {
            "gene_symbol": transcript_factors.index,
            "pc_median_tpm": pc_tpm_positive.loc[
                transcript_factors.index
            ].to_numpy(),
            "bayes_factor": transcript_factors.to_numpy(),
            "q": transcript_factors.attrs["q"],
            "T_q": transcript_factors.attrs["T_q"],
            "minimum_factor": transcript_factors.attrs["minimum_factor"],
            "background_size": transcript_factors.attrs["background_size"],
            "background_filter": "finite PC median TPM > 0",
        }
    ).sort_values(
        ["bayes_factor", "pc_median_tpm", "gene_symbol"],
        ascending=[False, False, True],
        kind="stable",
    )

    protein_table = pd.read_csv(PROTEIN_RESULTS, sep="\t")
    protein_factors = protein_table.set_index("gene_symbol")["bayes_factor"]

    initial_prior = pd.Series(
        0.5,
        index=universe_index,
        name="initial_prior",
    )
    protein_bayes_factors = protein_factors / 0.5
    transcript_bayes_factors = transcript_factors / 0.5
    after_protein = binary_bayes_factor_update(
        initial_prior, protein_bayes_factors, missing_factor=1.0
    )
    final_posterior = binary_bayes_factor_update(
        after_protein,
        transcript_bayes_factors,
        missing_factor=1.0,
    )

    combined = pd.DataFrame(index=universe_index)
    combined["initial_prior"] = initial_prior
    combined["protein_observed"] = combined.index.isin(protein_factors.index)
    combined["protein_factor"] = protein_factors.reindex(combined.index).fillna(0.5)
    combined["protein_bayes_factor"] = protein_bayes_factors.reindex(
        combined.index
    ).fillna(1.0)
    combined["posterior_after_protein"] = after_protein
    combined["pc_transcript_measurement_available"] = combined.index.isin(
        pc_tpm_all.index
    )
    combined["pc_transcript_zero_excluded"] = (
        pc_tpm_all.reindex(combined.index).eq(0.0).fillna(False)
    )
    combined["pc_transcript_observed"] = combined.index.isin(transcript_factors.index)
    combined["pc_transcript_factor"] = transcript_factors.reindex(combined.index).fillna(0.5)
    combined["pc_transcript_bayes_factor"] = transcript_bayes_factors.reindex(
        combined.index
    ).fillna(1.0)
    combined["posterior_probability"] = final_posterior
    combined["posterior_to_prior_ratio"] = final_posterior / initial_prior
    combined = combined.sort_values(
        ["posterior_probability", "gene_symbol"],
        ascending=[False, True],
        kind="stable",
    ).reset_index()

    TRANSCRIPT_RESULTS.parent.mkdir(parents=True, exist_ok=True)
    transcript_output.to_csv(TRANSCRIPT_RESULTS, sep="\t", index=False)
    combined.to_csv(COMBINED_RESULTS, sep="\t", index=False)
    finite_values = pc_tpm_all.iloc[np.flatnonzero(finite_mask)]
    universe_values = pc_tpm_all.reindex(universe_index)
    summary = {
        "input_workbook": str(WORKBOOK),
        "input_sheet": "Median TPM",
        "input_measurement": "PC (Median TPM,n=74)",
        "background_filter": "finite PC median TPM > 0",
        "finite_background_genes_before_zero_exclusion": int(finite_mask.sum()),
        "zero_background_genes_excluded": int((finite_values == 0.0).sum()),
        "negative_background_genes_excluded": int((finite_values < 0.0).sum()),
        "positive_background_genes_used": transcript_factors.attrs[
            "background_size"
        ],
        "signaling_candidates_with_positive_pc_tpm": transcript_factors.attrs[
            "candidate_size"
        ],
        "signaling_candidates_with_zero_pc_tpm_excluded": int(
            universe_values.eq(0.0).sum()
        ),
        "q": transcript_factors.attrs["q"],
        "T_q_positive_pc_median_tpm": transcript_factors.attrs["T_q"],
        "minimum_factor": transcript_factors.attrs["minimum_factor"],
        "candidates_above_neutral": int(
            (
                transcript_output["bayes_factor"]
                > transcript_factors.attrs["minimum_factor"]
            ).sum()
        ),
        "candidates_at_neutral": int(
            (
                transcript_output["bayes_factor"]
                == transcript_factors.attrs["minimum_factor"]
            ).sum()
        ),
        "node_prior_probability": 0.5,
        "zero_or_missing_pc_values_source_factor": 0.5,
        "zero_or_missing_pc_values_bayes_factor": 1.0,
        "posterior_probabilities_are_independent": True,
        "posterior_minimum": float(final_posterior.min()),
        "posterior_mean": float(final_posterior.mean()),
        "posterior_maximum": float(final_posterior.max()),
        "protein_factor_file": str(PROTEIN_RESULTS),
        "transcript_factor_file": str(TRANSCRIPT_RESULTS),
        "combined_output_file": str(COMBINED_RESULTS),
    }
    SUMMARY_OUTPUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"PC transcript background genes: {transcript_factors.attrs['background_size']}")
    print(f"PC signaling candidates observed: {transcript_factors.attrs['candidate_size']}")
    print(f"PC zero background genes excluded: {summary['zero_background_genes_excluded']}")
    print(f"PC q: {transcript_factors.attrs['q']}")
    print(f"PC T_q: {transcript_factors.attrs['T_q']:.12g}")
    print(f"PC candidates above neutral: {summary['candidates_above_neutral']}")
    print(f"combined posterior nodes: {len(combined)}")
    print("posterior probabilities are independent and are not normalized across nodes")
    print("top five combined nodes:")
    print(
        combined[
            [
                "gene_symbol",
                "protein_factor",
                "pc_transcript_factor",
                "posterior_probability",
                "posterior_to_prior_ratio",
            ]
        ].head(5).to_string(index=False)
    )


if __name__ == "__main__":
    main()
