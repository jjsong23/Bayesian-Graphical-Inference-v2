#!/usr/bin/env python3
"""Run node selection with all six collecting-duct streams enabled.

The generic GUI defaults remain unchanged for reproducibility. This script
creates the lab-specific Aqp2/collecting-duct profile by adding proteome and RNA
evidence for CCD, OMCD, and IMCD to the four existing node streams.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
GUI_DIR = PROJECT_ROOT / "gui"
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

from workflow_engine import (  # noqa: E402
    default_configuration,
    load_registry,
    normalize_configuration,
    select_nodes,
)


COLLECTING_DUCT_STREAMS = (
    "collecting_duct_proteome_ccd",
    "collecting_duct_proteome_omcd",
    "collecting_duct_proteome_imcd",
    "collecting_duct_rna_ccd",
    "collecting_duct_rna_omcd",
    "collecting_duct_rna_imcd",
)
OUTPUT_DIR = PROJECT_ROOT / "results/collecting_duct_node_selection"


def main() -> None:
    registry = load_registry(PROJECT_ROOT)
    supplied = default_configuration(registry)
    baseline_config = normalize_configuration(supplied, registry)
    baseline_factors, baseline_selected, baseline_summary = select_nodes(
        PROJECT_ROOT, registry, baseline_config
    )

    stage_rows: list[dict[str, object]] = [
        {
            "stage": 0,
            "added_stream": "existing four-stream node selection",
            "selected_proteins": baseline_summary["selected_protein_count"],
            "selected_total_nodes": len(baseline_selected),
            "new_proteins_at_stage": baseline_summary["selected_protein_count"],
        }
    ]
    previous_proteins = int(baseline_summary["selected_protein_count"])
    final_factors = baseline_factors
    final_selected = baseline_selected
    final_summary = baseline_summary
    for stage, stream_id in enumerate(COLLECTING_DUCT_STREAMS, start=1):
        supplied["node_streams"][stream_id]["enabled"] = True
        config = normalize_configuration(supplied, registry)
        final_factors, final_selected, final_summary = select_nodes(
            PROJECT_ROOT, registry, config
        )
        selected_proteins = int(final_summary["selected_protein_count"])
        stage_rows.append(
            {
                "stage": stage,
                "added_stream": stream_id,
                "selected_proteins": selected_proteins,
                "selected_total_nodes": int(len(final_selected)),
                "new_proteins_at_stage": selected_proteins - previous_proteins,
            }
        )
        previous_proteins = selected_proteins

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    final_factors.to_csv(
        OUTPUT_DIR / "node_posteriors_all_collecting_duct_streams.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    final_selected.to_csv(
        OUTPUT_DIR / "selected_node_universe_all_collecting_duct_streams.tsv",
        sep="\t",
        index=False,
    )
    stage_table = pd.DataFrame(stage_rows)
    stage_table.to_csv(
        OUTPUT_DIR / "cumulative_selection_by_stream.tsv", sep="\t", index=False
    )

    baseline_symbols = set(
        baseline_factors.loc[
            baseline_factors["gui_selected_in_graph"], "gene_symbol"
        ].astype(str)
    )
    final_symbols = set(
        final_factors.loc[
            final_factors["gui_selected_in_graph"], "gene_symbol"
        ].astype(str)
    )
    newly_selected = final_factors[
        final_factors["gene_symbol"].isin(final_symbols - baseline_symbols)
    ].copy()
    newly_selected.to_csv(
        OUTPUT_DIR / "newly_selected_proteins_from_collecting_duct.tsv",
        sep="\t",
        index=False,
    )

    final_config = normalize_configuration(supplied, registry)
    (OUTPUT_DIR / "configuration.json").write_text(
        json.dumps(final_config, indent=2), encoding="utf-8"
    )
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis": "existing node evidence plus six collecting-duct streams",
        "collecting_duct_streams": list(COLLECTING_DUCT_STREAMS),
        "baseline_selected_proteins": int(baseline_summary["selected_protein_count"]),
        "baseline_total_nodes_including_second_messengers": int(len(baseline_selected)),
        "final_selected_proteins": int(final_summary["selected_protein_count"]),
        "final_total_nodes_including_second_messengers": int(len(final_selected)),
        "newly_selected_proteins": int(len(final_symbols - baseline_symbols)),
        "prior_probability": final_summary["node_prior_probability"],
        "output_probability_cutoff": final_summary[
            "node_output_probability_cutoff"
        ],
        "posterior_probabilities_are_independent": final_summary[
            "posterior_probabilities_are_independent"
        ],
        "segment_dependence_note": (
            "CCD, OMCD, and IMCD are separate requested streams, but the three "
            "proteome streams share one experiment and the three RNA streams "
            "share one experiment; they are tagged with dependence groups."
        ),
        "cumulative_stages": stage_rows,
        "workflow_summary": final_summary,
    }
    (OUTPUT_DIR / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
