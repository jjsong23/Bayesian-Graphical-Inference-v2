# Changelog

This file summarizes user-visible and scientifically meaningful code changes.
The canonical takeover instructions remain in `CODEX_PROJECT_CONTEXT.md`; dated
implementation detail is retained under `docs/lab_notebook/`.

## 2026-09-15

### Compact structural screening and reference-interactome evidence

- Restricted new AlphaPulldown work to one total model prediction per pair:
  one prediction from `model_1_multimer_v3`.
- Moved each GPU task's heavy model output to node-local `$LSCRATCH`; only the
  validated ipTM table, top PDB, small metadata, and task manifest persist.
- Isolated feature/MSA objects by search run and added post-frontier cleanup
  that retains only next-frontier features and records every deletion/byte.
- Added validation and compact SQLite import for the earlier `ppi_screen`.
  Low ipTM values are retained, a nonempty top PDB is required, duplicates are
  audited, and the five-model legacy protocol remains distinctly labeled.
- Added exact-protocol-first structural cache lookup with explicit compatible
  protocol fallbacks and source/protocol provenance in the edge ledger.
- Added optional HuRI binary interaction evidence mapped from human Ensembl to
  mouse symbols. Reported pairs support edges; non-reporting is neutral by
  default and can become only a scoped, configurable weak negative assumption.
- Added Biowulf helpers to fetch HuRI, submit the legacy-screen import, and
  fast-forward an existing checkout from GitHub without touching ignored data.
- Added regression tests for prediction count, feature cleanup, legacy cache
  compatibility/validation, low-score retention, and HuRI non-report scoping.

## 2026-09-14

### Version 2 complete-pipeline workbench

- Connected the existing configurable node-selection engine directly to the
  Version 2 backward frontier search; a Version 2 run now begins with raw
  registered node evidence rather than a manually exported node universe.
- Added an exact handoff for all enabled inexpensive edge streams that preserves
  their Tq/reference controls, weights, and scope-aware negative evidence while
  avoiding a second application of the weight.
- Added target/receptor resolution, directionality export, sequence preparation,
  a Version 2-local persistent pair cache, and saved pipeline-level state.
- Added complete-pipeline CLI commands for configuration, initialization,
  single-stage execution, run-to-checkpoint, Biowulf submission, structural
  collection, status, and trace generation.
- Added a dedicated minimal GUI at `/v2.html` exposing all registered node and
  edge datasets and the complete initialize/step/submit/collect lifecycle.
- Deferred sequence validation until after cheap shortlist selection so an
  unmapped protein that is never submitted cannot block an otherwise valid run.
- Corrected the Biowulf submission manifest to contain real tab separators.
- Replaced repeated multi-gigabyte incremental-cache scans with indexed bulk
  pair lookup and isolated Version 2's writable incremental and structural
  caches under its ignored `runtime/` directory.

### Version 2 stepwise development inference

- Added an optional development state machine that advances exactly one
  auditable decision per command while preserving the production search's
  Bayesian equations, structural transformation, ranking, and pruning rules.
- Added complete candidate-enumeration, cheap-evidence, shortlist,
  structural-cache, combined-evidence, and frontier keep/reject ledgers.
- Added an append-only JSONL event history and a continuously refreshed,
  self-contained HTML trace showing the current frontier, next action, stage
  counts, and linked table previews.
- Added explicit structural import checkpoints: partial results remain waiting,
  and complete imports pause before Bayesian integration so raw scores can be
  inspected first.
- Added CLI controls (`init --development` and `trace`) and synthetic regression
  coverage for every stage through the first frontier transition.

## 2026-09-11

### Version 2 backward frontier inference

- Added a resumable backward beam search that evaluates cheap pair evidence
  only at the active frontier, retains 20 candidates per partial path, and
  globally retains five distinct upstream frontier nodes after structural
  rescoring.
- Added an auditable heuristic conversion from raw ipTM to a bounded structural
  BF, including negative support below a configurable reference score, while
  preserving raw scores separately from weights and transformation settings.
- Added a persistent SQLite pair cache, atomic run state, per-round candidate
  and pair manifests, loop prevention, signaling-role filtering, direction-aware
  backward traversal, and receptor/maximum-depth termination.
- Added portable preparation commands for audited node-symbol FASTA generation
  and for appending an external endpoint's HPA vector to the internal HPA pair
  catalog without reading or constructing a complete square adjacency matrix.
- Added portable NIH Biowulf Slurm scripts for cached AlphaPulldown feature
  generation, one exact pair per one-GPU array task, ipTM collection, and
  pause/import/resume operation without any AI assistant or API dependency.
- Added a complete design/methods plan, example configuration, root CLI entry
  point, and synthetic end-to-end tests of two-round search and cache reuse.

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
