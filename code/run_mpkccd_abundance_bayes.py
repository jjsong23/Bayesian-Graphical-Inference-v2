"""Apply signaling-universe Bayes factors to the mpkCCD abundance workbook."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from bayes_factors import signaling_bayes_factors


PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = PROJECT_ROOT / "data" / "mpkccd_protein_abundances_raw.xlsx"
UNIVERSE = (
    PROJECT_ROOT
    / "data"
    / "node_selection"
    / "mouse_signaling_symbols_liberal.txt"
)
OUTPUT = PROJECT_ROOT / "results" / "mpkccd_protein_abundance_bayes_factors.tsv"
SUMMARY_OUTPUT = (
    PROJECT_ROOT / "results" / "mpkccd_protein_abundance_bayes_summary.json"
)


def main() -> None:
    table = pd.read_excel(
        WORKBOOK,
        sheet_name="Original Data in webpage",
        usecols=["Gene Symbol", "Log10 Abundance"],
    )
    log10_abundances = pd.to_numeric(
        table.set_index("Gene Symbol")["Log10 Abundance"],
        errors="coerce",
    )
    linear_abundances = pd.Series(
        np.power(10.0, log10_abundances.to_numpy(dtype=float)),
        index=log10_abundances.index,
        name="linear_relative_abundance",
    )
    universe = UNIVERSE.read_text(encoding="utf-8").splitlines()

    factors = signaling_bayes_factors(linear_abundances, universe, q=0.75)
    result = pd.DataFrame(
        {
            "gene_symbol": factors.index,
            "log10_abundance": log10_abundances.loc[factors.index].to_numpy(),
            "linear_relative_abundance": linear_abundances.loc[
                factors.index
            ].to_numpy(),
            "bayes_factor": factors.to_numpy(),
            "q": factors.attrs["q"],
            "T_q": factors.attrs["T_q"],
            "minimum_factor": factors.attrs["minimum_factor"],
            "background_size": factors.attrs["background_size"],
            "input_transform": "10 ** supplied Log10 Abundance",
            "background_filter": "all finite back-transformed abundance values",
        }
    ).sort_values(
        ["bayes_factor", "linear_relative_abundance", "gene_symbol"],
        ascending=[False, False, True],
        kind="stable",
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, sep="\t", index=False)
    effective_cutoff = factors.attrs["T_q"] * math.sqrt(2.0 * math.log(2.0))
    summary = {
        "input_workbook": str(WORKBOOK),
        "input_sheet": "Original Data in webpage",
        "supplied_measurement": "Log10 Abundance",
        "analysis_measurement": "linear relative abundance",
        "input_transform": "10 ** supplied Log10 Abundance",
        "background_filter": "all finite back-transformed abundance values",
        "background_size": factors.attrs["background_size"],
        "signaling_candidates_observed": factors.attrs["candidate_size"],
        "q": factors.attrs["q"],
        "T_q_linear_relative_abundance": factors.attrs["T_q"],
        "effective_non_neutral_cutoff_linear_relative_abundance": effective_cutoff,
        "effective_non_neutral_cutoff_log10_abundance": math.log10(
            effective_cutoff
        ),
        "minimum_factor": factors.attrs["minimum_factor"],
        "candidates_above_neutral": int(
            (result["bayes_factor"] > factors.attrs["minimum_factor"]).sum()
        ),
        "candidates_at_neutral": int(
            (result["bayes_factor"] == factors.attrs["minimum_factor"]).sum()
        ),
        "output_file": str(OUTPUT),
    }
    SUMMARY_OUTPUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"background proteins: {factors.attrs['background_size']}")
    print(f"signaling candidates observed: {factors.attrs['candidate_size']}")
    print(f"q: {factors.attrs['q']}")
    print(f"T_q (linear relative abundance): {factors.attrs['T_q']:.12g}")
    print(f"effective non-neutral cutoff: {effective_cutoff:.12g}")
    print(f"minimum factor: {factors.attrs['minimum_factor']}")
    print(f"candidates above neutral: {summary['candidates_above_neutral']}")
    print(f"output: {OUTPUT}")
    print("top five:")
    print(result.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
