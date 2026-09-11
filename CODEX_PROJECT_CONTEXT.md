# Codex project context

Read this file before changing the project. It is the canonical short handoff for a new Codex session: it defines what the repository is for, where its data live, which scientific conventions are intentional, and which implementation decisions must not be silently reversed. The detailed scientific record is in `docs/full_project_methods.md`, `docs/edge_characterization_methods_891.md`, `docs/path_finding_target_extension.md`, and `docs/lab_notebook/`. Check `CHANGELOG.md` and the newest lab-notebook entry for changes made after the longer methods documents were written.

## One-sentence objective

Construct an auditable Bayesian graph of renal collecting-duct/principal-cell signaling participants, estimate undirected association probabilities between them from selectable evidence streams, and rank plausible signal-propagating paths from a user-selected start node to a target such as Aqp2.

## Current architecture

```text
9,170 candidate mouse proteins
        |
        | node-level Bayes factors
        v
selected node universe + curated second messengers
        |
        | edge-level Bayes factors on unordered pairs
        v
symmetric undirected probability matrix
        |
        | optional external-target row/column
        v
ontology + OmniPath traversal constraints
        |
        v
product-ranked loopless paths
        |
        +--> optional temporal-order audit
        +--> merged interactive top-path network
```

The validated seed graph contains **891 nodes: 871 proteins and 20 curated secondary messengers**. The GUI is no longer capped at 891: newly selected protein nodes are characterized against the active universe and added as needed. Raw observations for new unordered pairs are cached so later runs calculate only previously unseen pairs under the current source signature. The Git repository contains code, tests, configuration, and documentation; the large `data/`, `results/`, `outputs/`, and SQLite cache trees are restored from the companion OneDrive archive and are intentionally ignored by Git.

A lab-specific Aqp2/collecting-duct profile additionally enables proteome and
RNA abundance for CCD, OMCD, and IMCD. All six new streams plus the existing
four select **3,330 proteins**, or **3,350 total nodes** with second messengers.
This large profile is stored separately. The 891-node graph remains the
immutable edge-catalog seed, but the current generic node defaults select
1,051 proteins (1,071 nodes with messengers) after the site-centric
phosphoprotein revision and therefore extend beyond that seed dynamically.

## Bayesian conventions

- Node and edge evidence are represented as positive Bayes factors (BFs).
- Missing node and edge evidence is neutral by default. The GUI exposes a stage-wide fixed-absence rule and a per-dataset continuous alternative. The fixed rule gives an eligible nondetection one configurable `BF < 1` (default 0.5). The continuous rule removes the 0.5 likelihood floor, scores weak measurements below BF 1, and evaluates eligible nondetections as `x=0` at a small positive BF floor. Hypotheses outside the source's scope remain at `BF = 1`.
- For a probability prior `p`, integration is performed in odds space: `posterior_odds = prior_odds * product(BF_i ** weight_i)`, then converted back to probability.
- The common continuous-evidence transformation uses a source-specific threshold `T_q`. The GUI permits independent normalization multipliers for each applicable dataset.
- STRING and STITCH instead use score odds relative to configurable reference scores, not the Gaussian-complement threshold kernel.
- The historical phrase "complement of the minimum Bayes factor" refers to converting complementary tail evidence into positive likelihood/Bayes support while preserving a neutral floor. Confirm the exact implementation in `code/bayes_factors.py` before describing equations.
- Node and edge posteriors are independent Bernoulli updates. Every protein and every edge begins at probability 0.5 by default, and posterior values are not normalized across candidates or pairs. Do not describe ranked whole-path products as calibrated biological probabilities.

## Node selection streams

The selectable streams currently include:

- mpkCCD protein abundance;
- principal-cell RNA expression, using raw abundance with zero-valued observations excluded from the evidence-background calculation;
- inferred kinase activity from PKA-knockout phosphosite changes; and
- protein-level differential phosphorylation after PKA deletion;
- selective PKA-Cα-knockout and PKA-Cβ-knockout protein-level differential
  phosphorylation, registered as two optional streams;
- rat proteome abundance for CCD, OMCD, and IMCD, harmonized to mouse with
  audited Ensembl release 116 orthology; and
- mouse renal-tubule RNA abundance for CCD, OMCD, and IMCD.

For every protein, the active model compares `present in the signaling system`
with `not present`. The default prior assigns probability 0.5 to each
hypothesis. Protein/RNA source scores have a neutral floor of 0.5 and are
divided by 0.5 to obtain a conventional neutral Bayes factor of 1; kinase and
phosphoprotein multipliers are already neutral at 1. Weighted factors multiply
each protein's prior odds independently. A protein is selected when its
posterior is strictly greater than the configured node cutoff, which defaults
to 0.5. There is no across-protein sum-to-one normalization.

The optional `penalize_unobserved` node policy is disabled in the reproducibility default. When enabled, the engine replaces the neutral BF for each eligible nondetection with `unobserved_bayes_factor` (default 0.5) before applying the stream weight. Repeated nondetection BFs multiply the node's odds downward. Eligibility is stream-specific: non-kinases always remain neutral under kinase-activity evidence, whereas protein, RNA, phosphoprotein-response, and collecting-duct abundance streams treat all protein candidates as eligible. The engine writes eligible, observed, penalty-applied, source-BF, effective-BF, and weighted-log-BF audit columns for every active node stream. A strictly positive BF is required; posteriors can approach zero but are not made irreversibly equal to zero.

Every node dataset also has an independent `continuous_negative_evidence` switch. When enabled, `BF=max(epsilon,(1-exp(-0.5*(x/Tq)^2))/0.5)` and eligible nondetections use `x=0`; the fixed penalty is skipped for that stream. The positive stored Tq is retained because zero-padding sparse 9,170-candidate backgrounds can make q75 equal zero. See `docs/continuous_negative_evidence.md`.

With the option disabled, proteins with no detected phosphosite remain neutral under phosphoprotein evidence. Curated secondary messengers are inserted separately because they are not gene products measured by the protein/RNA screens and are never subjected to protein nondetection penalties.

For the six collecting-duct streams, zero is a nondetection and receives BF=1 when the optional penalty is off, or the configured nondetection BF when it is on;
each q75 uses all positive finite source values in that assay segment before
candidate filtering. Rat mapping uses UniProt first and rat-symbol fallback,
prefers high-confidence orthologs, flags low-confidence fallbacks, retains valid
multi-ortholog relationships, and aggregates duplicate contributions by the
maximum abundance. The three segment streams within each assay share a
dependence group. Entry points are
`code/node_selection/build_collecting_duct_evidence.py` and
`code/node_selection/run_collecting_duct_node_selection.py`.

The selective PKA-subunit streams are built by
`code/node_selection/build_pka_subunit_ko_evidence.py`. Duplicate source rows
are collapsed by median signed LFC for each UniProt/site-pattern key. Each
comparison uses its own all-source-site absolute-LFC background. Every site is
scored against the comparison's single q75, and a protein inherits its
strongest site-level factor; there is no adjustment for its number of detected
sites and no accumulation across sites. The double-KO phosphoprotein-response
stream follows the same rule. Pmod is audit-only. Both streams are disabled
in the generic default and share dependence group
`pka_subunit_ko_phosphoproteomics`. Read
`docs/pka_subunit_ko_node_evidence.md` before changing this logic.

The canonical symbol field is `symbol`; duplicate symbols were removed when the original liberal signaling universe was consolidated.

## Edge characterization streams

Bayesian edge characterization remains **undirected and symmetric**. An unordered biological pair is represented once in pair tables and mirrored across the edge-probability matrix; self-edges are excluded. Path inference has a separate conservative propagation-direction layer. Versioned ontology-role rules are applied first. Mapped OmniPath mouse core source-target records can uniquely orient a pair only when ontology left it unresolved; they cannot reopen or reverse an ontology-disallowed traversal. A bidirectional OmniPath record leaves an unresolved pair traversable both ways. The allowed traversal keeps the original undirected edge probability, so directionality does not re-estimate edge existence. Activation/inhibition sign remains audit metadata and is not used to score or orient a path.

Selectable streams include:

- mpkCCD subcellular localization;
- observed-phosphosite-restricted KinasePredictor kinase–protein support;
- STRING protein association;
- Human Protein Atlas localization (primary or high-confidence alternative);
- OmniPath interactions, collapsed to undirected support;
- STITCH and curated secondary-messenger associations; and
- optional scaffold-mediated binary closure.

The optional `penalize_unsupported` edge policy is disabled in the reproducibility default. When enabled, a pair that is eligible for an active source but lacks a non-neutral source record receives `unsupported_bayes_factor` (default 0.5) before the stream weight is applied. Scope is source-specific: both proteins must have the relevant localization profiles for mpkCCD/HPA; KinasePredictor requires a kinase and a protein with an observed scorable phosphosite; STRING requires two mapped proteins; OmniPath considers protein–protein pairs; and STITCH considers curated-messenger-to-mapped-protein pairs. Scaffold closure has no negative-absence rule. The update remains symmetric and is also applied to eligible external-target edges. Database absence is an optional modeling assumption, not proof that a biological interaction is impossible.

Every primary edge dataset also has an independent continuous-negative switch. It retains measured low factors below 1 and fills eligible no-record pairs with the numerical floor; out-of-scope pairs stay neutral. It is implemented for internal, incrementally cached, calibration, and external-target edges. Scaffold closure remains positive-only.

Localization compatibility/adjacency matrices for mpkCCD, HPA, and COMPARTMENTS are stored with the archived data so a scientist can validate every assumed compartment relationship. COMPARTMENTS is used in the dedicated colocalization/AlphaFold candidate-filtering workflow. Experimental PPI tiers assembled from STRING, BioGRID, and IntAct support the prior-knowledge filter used before structural prediction.

The incremental raw-pair cache is `data/edge_characterization/incremental_edge_cache/edge_pair_cache.sqlite3`. It is a performance cache, not the only scientific record: each run also writes explicit configuration and audit outputs.

## Scaffold closure: do not restore the rejected heuristic

The accepted optional rule is deliberately simple:

1. Use only the **pre-closure** graph.
2. A protein is an anchor to an exact `adaptor_scaffold` node only when their edge probability is strictly greater than the selected cutoff (default `0.90`).
3. Two proteins qualify if they share at least one such scaffold.
4. Every qualifying pair receives fixed support likelihood `0.90`, equivalent to `BF = 0.90 / 0.50 = 1.8`; every other pair receives neutral `BF = 1`.
5. Apply the rule once. Closure-derived edges never become new anchors.

There is **no scaffold-degree or “promiscuity” penalty, no anchor-excess score, no noisy-OR aggregation, and no empirical scaffold `T_q`**. Those earlier heuristics were rejected as insufficiently biologically justified. Closure is off by default and must be described as hypothesis-generating proximity/co-complex support, not proof of a direct PPI.

Validation output: `results/gui_runs/scaffold_binary_closure_validation_20260804_v2/`. At defaults, 195 scaffold-tagged proteins yielded 14,308 qualifying protein–scaffold anchors, 205,920 qualifying protein pairs, and 81,323 pairs newly raised above probability 0.5.

## Target extension and path inference

`code/path_finding/build_target_adjacency_vector.py` characterizes an external mouse protein against all active nodes and appends a symmetric row and column while preserving the existing matrix. `code/path_finding/ontology_directionality.py` combines the auditable role catalog in `ontology_direction_rules.json` with mapped OmniPath directions and emits a partially directed propagation matrix without changing allowed edge probabilities. `code/path_finding/find_ranked_paths.py` ranks loopless paths by the product of allowed traversal probabilities, implemented with additive `-log(p)` costs and a hop-limited Yen/Dijkstra search.

Path intermediates are controlled by ontology-derived role classes in the GUI. The sensitivity-oriented interface exposes all current role labels and lets the user select which may propagate signals. Multi-role kinase/scaffold proteins remain eligible whenever they match any selected relay role. The separate strict scaffold override excludes every scaffold-tagged intermediate and should remain optional. Start and target nodes are endpoint exemptions.

## Optional temporal path validation

`code/path_finding/temporal_path_ranking.py` annotates each current ranked path with response-order evidence from the raw 1/2/5/15-minute dDAVP phosphoproteomic workbook at `data/phospho_data_original.xlsx`. Replicate mode uses intensity-moderated site/timepoint variances, selects one site per gene by peak moderated signal-to-noise, applies a within-gene Bonferroni correction by default, and propagates response-time uncertainty with Monte-Carlo draws. The GUI exposes the prior df, alpha, draw count, seed, p-adjustment method, and minimum measured-node coverage. This layer is disabled by default.

Temporal validation is post-path annotation, not another edge Bayes factor. It never changes the symmetric edge-existence matrix, ontology directionality, path probabilities, or primary Bayesian rank. It writes `ranked_paths_temporal.tsv`, a gene-level response audit, the empirical intensity–variance trend, and a dedicated JSON summary. Paths with insufficient measured nodes remain temporally unrankable; that is missing temporal information, not evidence that the path is absent. Detailed methods and reporting cautions are in `docs/temporal_path_validation.md`. `code/path_finding/temporal_sensitivity.py` audits dependence on variance prior df and alpha.

## Workbench entry points

- Portable root launcher: `launch.py`
- Start: `gui/run_workbench.ps1`
- Detached start: `gui/run_workbench.ps1 -Background`
- Detached launcher: `gui/launch_workbench.py`
- Server: `gui/server.py`
- Main workflow: `gui/workflow_engine.py`
- Per-hypothesis evidence reconstruction: `gui/evidence_inspector.py`
- Evidence metadata: `gui/evidence_registry.json`
- Browser client: `gui/web/`
- Bayesian primitives: `code/bayes_factors.py`
- Incremental edge cache: `code/edge_characterization/incremental_edge_cache.py`

Use the detached start when the backend must survive its launching terminal or
Codex execution session. It performs an `/api/config` health check, records
`gui/.workbench.pid`, and reuses an already healthy listener. A visible browser
tab alone does not mean the Python backend is still running.

The registry drives the cards and generic factor combination. A genuinely new raw-data transformation still needs a dedicated, audited preprocessing module that emits the standard factor-table shape.

## Run outputs and current validated records

- GUI runs: `results/gui_runs/<immutable_run_id>/`
- Scaffold-rule validation: `results/gui_runs/scaffold_binary_closure_validation_20260804_v2/`
- Target extensions: `results/path_finding/target_extensions/`
- Ranked paths: `results/path_finding/ranked_paths/`
- Current methods: `docs/full_project_methods.md`
- Current edge methods: `docs/edge_characterization_methods_891.md`
- Collecting-duct node results: `results/collecting_duct_node_selection/`
- Tq-aware end-to-end ablation: `code/sensitivity_analysis/analyze_end_to_end_ablation.py`
- Recent decisions: newest dated entry under `docs/lab_notebook/`

The full data/results tree is not stored in GitHub. Restore it from the dated complete-project archive described in `DATA_AND_RESULTS.md`, preserving paths.

## 2026-09-10 current implementation handoff

The registry currently contains **12 node streams and 8 edge streams**; the
generic profile enables 4 node and 6 edge streams. The following capabilities
postdate the longer August methods summary and are part of the current code:

1. **Tq-aware end-to-end ablation.**
   `code/sensitivity_analysis/analyze_end_to_end_ablation.py` tests every
   registered stream by default, including disabled optional sources in a
   matched add-one comparison. Mutually exclusive alternatives, such as the two
   HPA variants, are compared in replacement contexts. Only the focal
   Tq/reference multiplier is swept; weights, priors, cutoffs, non-focal Tq
   values, and calibration results are held fixed. Node ablation repeats node
   selection and rebuilds edges, directionality, and paths. Edge ablation keeps
   the node set fixed but rebuilds downstream edge/path results. Outputs include
   posterior changes, unique-edge cutoff flips, exact-route rank changes,
   conditional Bernoulli KL divergence, Jaccard similarities, robustness
   ranges, PNG summaries, checkpoints, and a methods README. Calibration is
   deliberately not refit because that would measure optimizer compensation
   rather than the marginal effect of removing the stream.

2. **Per-hypothesis evidence inspection.**
   A completed GUI run can inspect one node or one unordered edge. The ledger
   shows every registered stream—including disabled ones—its support/refute/
   neutral/disabled call, assay scope, applied BF, weight, weighted change in
   log2 odds, normalization control, and missing-data rule. The backend
   reconstructs the posterior from the displayed terms and reports any mismatch
   with the stored posterior. New runs retain a compact log2(BF) distribution
   for every active stream. Each ledger row shows the selected BF's empirical
   percentile/tie interval (or a histogram-based estimate for a continuous edge
   stream), and selecting the row opens a marked distribution plot. Node scopes
   contain all modeled proteins; edge scopes contain all unique unordered graph
   pairs. The API is
   `/api/jobs/<job_id>/evidence/node?symbol=...` or
   `/api/jobs/<job_id>/evidence/edge?node_a=...&node_b=...`. At present, the
   server must still know the completed job in its in-memory job registry; a
   server restart preserves run files but not inspector access to that old job.

3. **Ontology plus OmniPath directionality.**
   Directionality is enabled by default for path search, and OmniPath direction
   use has a separate enabled-by-default checkbox beneath it. Ontology has
   precedence. OmniPath adds a disallowed reverse traversal only for an
   ontology-unresolved, uniquely directed pair. Every run writes the propagation
   matrix, complete edge audit, rules, class-pair catalog, and mapped OmniPath
   direction evidence. The validated default run oriented 35,061 of 224,434
   supported unique edges (15.6%): 33,720 by ontology only, 813 by OmniPath only,
   and 528 by both. These are traversal constraints, not new edge probabilities
   and not activation/inhibition calls.

4. **Merged predicted-network view.**
   Path runs write `top_path_network.json` and expose the same payload in the GUI
   preview. The backend unions up to the top 50 ranked paths, collapses repeated
   unordered relationships, preserves contributing path ranks, and marks an
   arrow only when exactly one traversal is allowed. The browser shows Top 5,
   Top 10, Top 25, or all available visualized paths in one nonlinear layout.
   Node and edge geometry is fixed; a continuous color intensity encodes
   posterior evidence magnitude on a capped log-odds scale. Green supports, red
   refutes, gray is near the 0.5 prior, and orange denotes a curated messenger
   or external target without a node posterior. Selecting a node or edge also
   loads its per-stream evidence ledger. Selecting a ranked path highlights the
   complete route and initially loads its lowest-probability edge. The latest
   local default validation merged 25 paths into 19 nodes and 38 relationships,
   12 uniquely directed and 26 unresolved.

5. **Small incremental-pair performance path.**
   The SQLite cache now uses exact primary-key probes and exact pair-factor
   loads when at most 50,000 dynamic unordered pairs are requested. Larger
   expansions retain the set-based cache scan. This changes query strategy, not
   scientific scoring or cache semantics.

6. **Readable calibration output.**
   When positive-control calibration is enabled, the result page now renders
   separate node and edge tables. For each optimized primary stream they show
   the starting and fitted weight together with the starting, preferred, and
   fitted Tq/reference multiplier. The underlying calibration JSON, TSV tables,
   and optimizer trace are unchanged and remain the definitive downloadable
   audit.

The latest dated implementation record is
`docs/lab_notebook/2026-09-10_gui_evidence_visualization.md`; the preceding
multi-day record is `docs/lab_notebook/2026-08-26_to_2026-09-09.md`.

## First checks in a new Codex session

1. Read this file, `CHANGELOG.md`, and the most recent lab-notebook entry.
2. Inspect `git status` and preserve unrelated user changes.
3. Confirm the companion data archive has been restored before running the GUI.
4. Run the unit tests listed in `README.md` before and after logic changes.
5. For a scientific claim, trace the exact factor column and source metadata through `gui/evidence_registry.json`, the preprocessing script, and the immutable run configuration.
6. Record substantive method changes in both the relevant detailed methods document and a dated lab-notebook entry.

## 2026-08-25 sensitivity and shareable release handoff

The complete current overview is `docs/PROJECT_COMPLETE_SUMMARY_2026-08-25.md`.
The new score-distribution renderer,
`code/sensitivity_analysis/render_evidence_score_distributions.py`, saves one
PNG for every node/edge evidence stream plus consolidated figures and a TSV
summary under the dated sensitivity-analysis result directory. Eligible
nondetections are shown at zero, item-specific Tq populations are shown as a
median plus 10th-90th percentile band, and non-Tq references are labeled
explicitly.

`code/release/build_shareable_packages.py` creates a complete PNG archive and a
runnable lab-sharing application archive under ignored `deliverables/`. The app
archive intentionally omits the 3.35 GB disposable SQLite cache and old outputs,
but includes every source required for the 891-node seed workflow and for
incremental characterization of newly added nodes. Read
`docs/SHAREABLE_RELEASE_2026-08-25.md` before changing the package scope.
The portable root launcher is `launch.py`. Release archive names and internal
paths are deliberately short so Windows Explorer can extract them from the
project's deep `deliverables` directory without crossing the 260-character
legacy path limit.

## Scientific cautions

- Treat evidence-stream dependence explicitly; kinase activity and phosphoprotein response share the PKA-knockout experiment, and scaffold closure is derived from the integrated graph.
- HPA primary and high-confidence variants are mutually exclusive alternatives.
- A database association can represent functional linkage or co-complex support rather than direct physical binding.
- Preserve raw score, mapping, direction/sign, threshold, and source columns in audits even when the active graph simplifies them.
- Do not silently change node identifiers, species mappings, zero handling, threshold backgrounds, neutral defaults, edge directionality, or path eligibility semantics.
- Do not describe temporal validation as updating a Bayes factor or as re-ranking the canonical path table. Its `temporal_evidence_rank` is separate, and the primary Bayesian `rank` is preserved.

## Optional positive-control parameter calibration (2026-08-17)

The GUI now has a disabled-by-default stage 00 that fits node weights/Tq scales
and edge weights/Tq scales independently to user-supplied known-present controls.
Unknown nodes and pairs are unlabeled, not negatives. The objective is mean
positive-control negative log probability plus one standardized quadratic
regularizer. It uses deterministic bounded multi-start SciPy Powell; weights are
bounded 0–3, ordinary Tq/reference multipliers 0.25–4, and phosphoproteomic
multipliers 0.05–1. Preferred weights are 1. Every primary stream has its own
GUI-editable preferred Tq/reference multiplier: ordinary streams default to 1,
while PKA-KO kinase activity, PKA-KO phosphoprotein response, both selective
PKA-subunit phosphoprotein responses, and KinasePredictor default to 0.1.
Preferred values are regularization anchors, not direct scoring controls.
Derived scaffold closure is held fixed. Node fitting precedes node selection; edge fitting follows endpoint
extension and incremental-pair caching. Read
`docs/positive_control_parameter_calibration.md` before changing this logic.
