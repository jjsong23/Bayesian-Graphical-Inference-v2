# Data and results handoff

The GitHub repository intentionally omits `data/`, `results/`, `outputs/`, large spreadsheets, compressed source files, and the incremental SQLite cache. They are retained in the complete project ZIP prepared for OneDrive.

## Restore the complete workspace

1. Clone the GitHub repository.
2. Extract the complete project archive.
3. Copy or merge the archive's `data/`, `results/`, `outputs/`, and root-level `.xlsx` files into the repository root.
4. Preserve the relative paths. The evidence registry and analysis scripts resolve inputs from those paths.
5. Start the workbench with `python launch.py` (portable) or
   `gui/run_workbench.ps1` (PowerShell convenience wrapper).

The archive is a point-in-time research record. Do not overwrite an older archive when publishing a revised one; use a dated filename and retain its SHA-256 checksum.

## Data layout

- `data/node_selection/`: signaling-universe definitions, GO terms,
  abundance/transcript evidence, node-level Bayes-factor tables, and the
  selective PKA-Cα/PKA-Cβ knockout source and processed node evidence under
  `pka_subunit_ko/`.
- `data/kinase_predictor/` and `data/pka_ko/`: kinase-model resources and PKA-knockout phosphoproteomic inputs.
- `data/edge_characterization/`: edge evidence, raw-source downloads, compatibility tables, factor catalogs, the mapped OmniPath direction source, and the incremental unordered-pair cache.
- `data/colocalization/`: localization profiles and scientist-reviewable compartment compatibility matrices.
- `data/experimental_ppi/`: STRING/BioGRID/IntAct-derived prior-support resources and evidence tiers.
- `data/phospho_data_original.xlsx`: raw 1/2/5/15-minute dDAVP/vehicle replicate phosphoproteomics used only by the optional temporal path-validation stage.
- `results/`: versioned node, edge, colocalization, AlphaFold-candidate, path, and GUI-run outputs.
- `outputs/`: auxiliary generated output files.

## External sources represented

The working data tree includes derived or downloaded material from resources such as Gene Ontology, UniProt, Human Protein Atlas, STRING, STITCH, OmniPath, COMPARTMENTS, BioGRID, IntAct, PhosphoSitePlus-derived site annotations where available, KinasePredictor, mpkCCD localization/proteomics, principal-cell RNA sequencing, and the PKA-knockout phosphoproteomic study data supplied for this project.

This file is an inventory, not a redistribution license. Before public sharing, verify the current terms, required citations, and redistribution permissions for every external source. In particular, keep restricted datasets in controlled OneDrive storage if their terms do not permit GitHub redistribution.

## Reproducibility convention

Each GUI run creates a new directory under `results/gui_runs/` with its configuration, node posterior table, selected universe, edge matrix, supported edge list, path tables, eligibility audit, and summary. Path runs also write `top_path_network.json`, a compact union of up to 50 ranked paths used by the interactive network view. Directional runs write the partially directed propagation matrix and complete ontology/OmniPath direction audits. Temporal runs additionally write the annotated path table, gene-level response audit, intensity–variance trend, and temporal summary. Treat these directories as immutable records. The SQLite edge cache stores raw pair-level observations for reuse; changes to weights or normalization settings rescore cached observations without recomputing unchanged source lookups.

Tq-aware leave-one-stream-out analyses write under
`results/sensitivity_analysis/`. Their tabular outputs distinguish node losses,
unique undirected edge changes, exact path-route changes, conditional Bernoulli
KL divergence, and Tq/reference-multiplier robustness. Generated sensitivity
results and PNGs belong in the companion archive rather than Git history.
