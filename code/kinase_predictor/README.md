# PKA-KO KinasePredictor annotation

This folder contains the reproducible workflow used to add KinasePredictor v0.8 sequence-motif predictions to the PKA-knockout phosphoproteomic dataset.

## Project files

- `../../data/pka_ko/PKA-KO_database_raw.xlsx` — unchanged input workbook.
- `../../data/kinase_predictor/v0.8/kinase-predictor-v0.8-official.zip` — official ESBL KinasePredictor v0.8 archive.
- `../../data/kinase_predictor/v0.8/official_package/` — extracted official executable, README, and matrices.
- `reference_gui/` — the separately supplied KinasePredictor Python GUI source.
- `score_pka_ko_phosphosites.py` — batch implementation of the v0.8 matrix-scoring method.
- `build_pka_ko_kinase_predictions.mjs` — creates and verifies the annotated Excel workbook.
- `../../results/kinase_predictor/pka_ko_scored_rows_top10.json` — scored intermediate records containing ranks 1–10.
- `../../results/kinase_predictor/PKA-KO_database_with_KinasePredictor_top10.xlsx` — top-10 annotated workbook.
- `../../results/kinase_predictor/qa_top10/` — rendered previews used for visual verification.

## Scoring

After removing the asterisk from a valid centered 13-residue sequence, each residue is used to select one value from each kinase's position-specific matrix. The 13 selected values are summed. Kinases are ranked from greatest to lowest score using the serine/threonine or tyrosine matrix set according to the center residue, and ranks 1–10 are retained.

Rows labeled `Multiple sites` are retained but not scored because they do not identify a unique centered phosphosite sequence.

The raw score is a sequence-motif match, not a probability, posterior probability, or Bayes factor.

## Reproduction order

1. Run `score_pka_ko_phosphosites.py` with Python containing pandas and NumPy.
2. Run `build_pka_ko_kinase_predictions.mjs` with the Codex bundled Node runtime and `@oai/artifact-tool` available.

Official source: https://esbl.nhlbi.nih.gov/Databases/Kinase_Logos/KinasePredictor.html
