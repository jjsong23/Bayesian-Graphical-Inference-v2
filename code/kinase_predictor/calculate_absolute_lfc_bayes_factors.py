"""Calculate evidence factors from the absolute kinase-level median LFC.

The input vector contains one value for each kinase retained in the top-10-hit
analysis: abs(median raw phosphosite LFC). All 185 values form the empirical
background for T_q, matching the project's existing Bayes-factor-like scoring
function.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = PROJECT / "results" / "kinase_predictor" / "top10_hit_analysis"
OUTPUT_DIR = PROJECT / "results" / "kinase_predictor" / "absolute_lfc_bayes_factors"
Q = 0.75
MINIMUM_FACTOR = 0.5


def main() -> None:
    prior_results = json.loads((ANALYSIS_DIR / "kinase_results.json").read_text(encoding="utf-8"))
    background = pd.Series(
        {row["kinase"]: abs(float(row["median_raw_log2_change"])) for row in prior_results},
        dtype=float,
    )
    t_q = float(background.quantile(Q))
    results = []
    for row in prior_results:
        absolute_lfc = abs(float(row["median_raw_log2_change"]))
        raw_factor = 1.0 - math.exp(-0.5 * (absolute_lfc / t_q) ** 2)
        bayes_factor = max(MINIMUM_FACTOR, raw_factor)
        results.append({
            "kinase": row["kinase"],
            "kinase_family": row["kinase_family"],
            "unique_hit_sites": row["unique_hit_sites"],
            "signed_median_raw_lfc": row["median_raw_log2_change"],
            "absolute_lfc": absolute_lfc,
            "absolute_lfc_over_T_q": absolute_lfc / t_q,
            "unfloored_factor": raw_factor,
            "bayes_factor": bayes_factor,
            "at_minimum_factor": bayes_factor == MINIMUM_FACTOR,
            "signed_profile_classification": row["classification"],
            "signed_profile_bh_fdr": row["bh_fdr"],
        })

    results.sort(key=lambda r: (-r["bayes_factor"], -r["absolute_lfc"], -r["unique_hit_sites"], r["kinase"]))
    for rank, row in enumerate(results, start=1):
        row["bayes_factor_rank"] = rank

    summary = {
        "kinases_scored": len(results),
        "background_size": len(background),
        "q": Q,
        "T_q": t_q,
        "minimum_factor": MINIMUM_FACTOR,
        "kinases_above_minimum_factor": sum(not row["at_minimum_factor"] for row in results),
        "kinases_at_minimum_factor": sum(row["at_minimum_factor"] for row in results),
        "kinase_statistic": "absolute value of each kinase's median raw phosphosite log2(PKA-null / PKA-intact) change",
        "background_source": "all 185 finite kinase-level absolute LFC values",
        "formula": "max(0.5, 1 - exp(-0.5 * (absolute_LFC / T_q)^2))",
        "input_analysis": str(ANALYSIS_DIR / "kinase_results.json"),
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "kinase_absolute_lfc_bayes_factors.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    (OUTPUT_DIR / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.DataFrame(results).to_csv(OUTPUT_DIR / "kinase_absolute_lfc_bayes_factors.tsv", sep="\t", index=False)
    print(json.dumps({"summary": summary, "top_15": results[:15]}, indent=2))


if __name__ == "__main__":
    main()
