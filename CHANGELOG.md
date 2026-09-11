# Changelog

This file summarizes user-visible and scientifically meaningful code changes.
The canonical takeover instructions remain in `CODEX_PROJECT_CONTEXT.md`; dated
implementation detail is retained under `docs/lab_notebook/`.

## 2026-09-10

### Visual evidence and calibration transparency

- Replaced variable network node size and edge width with a fixed geometry and
  a continuous red–gray–green posterior color scale. This separates evidence
  strength from graph topology and makes similarly supported edges easier to
  compare.
- Connected network selections directly to the completed-run evidence
  inspector. Selecting a node or edge now loads its ledger; selecting a ranked
  path highlights the route and initially inspects its lowest-probability edge.
- Added compact log2(BF) distributions for every active evidence stream. Each
  inspector row reports the selected BF's percentile/tie interval, and its
  detail view marks that BF against all modeled node hypotheses or all unique
  undirected graph pairs from the same run.
- Added readable node and edge calibration tables showing each stream's
  starting and fitted weight plus starting, preferred, and fitted Tq/reference
  multiplier. The existing JSON, TSV, and optimizer-trace outputs remain the
  machine-readable audit source.

## 2026-09-09

### Direction-aware path inference

- Added mapped OmniPath mouse core source-target directions to the propagation
  layer while preserving undirected Bayesian edge-existence probabilities.
- Retained ontology precedence: OmniPath orients only pairs that ontology left
  unresolved and cannot reopen or reverse an ontology-disallowed traversal.
- Kept bidirectional or conflicting evidence unresolved and traversable both
  ways.
- Added a separate GUI control for OmniPath directions and expanded the
  direction audit with source-specific counts, resources, references, and
  ontology/OmniPath conflicts.

### Evidence transparency

- Added a completed-run inspector for any modeled node or unordered edge.
- The inspector reports all registered streams, including disabled streams,
  with support/refute/neutral status, scope, BF, weight, weighted log2-odds
  contribution, normalization, and missing-data handling.
- Added an arithmetic reconciliation between the displayed update ledger and
  the posterior stored by the completed run.

### Predicted-network visualization

- Added a merged network view that overlays up to 50 ranked paths and collapses
  repeated relationships.
- Added Top 5/10/25/all controls, deterministic re-layout, neighborhood
  highlighting, exact node/edge details, and accessible keyboard selection.
- Encoded node and edge evidence magnitude on a capped log-odds scale so highly
  supported relationships remain distinguishable. Green supports, red refutes,
  gray is near the prior, and orange marks nodes without a Bayesian node
  posterior.
- Added arrowheads only where the propagation matrix allows exactly one
  direction. The compact backend payload is saved as `top_path_network.json`.

### Sensitivity analysis

- Added Tq-aware leave-one-evidence-stream-out analysis across all 12 registered
  node streams and all 8 registered edge streams.
- Optional streams use matched add-one contexts; mutually exclusive source
  alternatives use matched replacement contexts.
- Node-stream ablation propagates through node selection, edge rebuilding,
  directionality, and path ranking. Edge-stream ablation rebuilds edges and
  paths with the selected nodes fixed.
- Added conditional Bernoulli KL divergence, cutoff flips, unique-edge and path
  Jaccard measures, path-rank changes, Tq robustness ranges, checkpoints, PNG
  summaries, and machine-readable methods metadata.

### Performance and auditability

- Added an exact-pair cache query path for small dynamic expansions, avoiding a
  multi-gigabyte SQLite cache scan when at most 50,000 pairs are requested.
- Added per-stream edge log-odds contribution collection without changing the
  posterior integration equation.
- Expanded regression coverage for OmniPath precedence/conflicts, merged path
  graphs, evidence reconstruction, all-stream ablation contexts, GUI controls,
  and small incremental-pair cache behavior.

## 2026-08-25

- Added all-stream score-distribution and Tq-response figures.
- Added reproducible Windows-safe PNG and runnable-application release archives.
- Renamed the portable application entry point to `launch.py`.
- Added the complete project summary and explicit Git-versus-OneDrive data
  handoff documentation.

## Earlier development

Earlier scientific decisions—including the 0.5 independent priors, optional
negative evidence, collecting-duct evidence streams, site-centric
phosphoproteomic scoring, regularized positive-control calibration, scaffold
closure, target extension, and temporal path validation—are recorded in the
dated lab notebook and detailed methods files under `docs/`.
