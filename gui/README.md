# Graphical Bayesian Inference workbench

This local browser application runs the project from node selection through
edge integration, external-target extension, and ranked path finding. It is
local by design: arbitrary target characterization reads the project's raw
mouse datasets and executes the audited Python analysis code on this computer.

## Start

From PowerShell:

```powershell
cd C:\Users\songjj\Documents\Codex\2026-07-21\un\graphical_bayesian_inference
.\gui\run_workbench.ps1
```

The foreground command keeps the server attached to that PowerShell window.
Closing the window stops the backend even though an already-open browser tab
may remain visible. To keep the server running independently, use:

```powershell
.\gui\run_workbench.ps1 -Background
```

The background launcher starts a hidden detached process, waits for the
configuration API to become healthy, writes its PID to `gui/.workbench.pid`,
and reuses an already healthy server rather than creating a duplicate.

The workbench opens at `http://127.0.0.1:8765`. It binds only to the local
loopback interface by default.

If the backend stops, the page now reports **Local analysis server is not
running** rather than misclassifying the browser network error as an invalid
scientific configuration. Restart the server and choose **Reconnect to local
server**.

## Workflow

1. **Select nodes.** Enable any subset of protein abundance, principal-cell
   RNA, kinase activity, double-KO phosphoprotein response, the separately
   selectable PKA-Cα- and PKA-Cβ-KO phosphoprotein-response streams, and the six lab-specific
   collecting-duct streams (proteome and RNA for CCD, OMCD, and IMCD). Each enabled
   factor can be raised to a nonnegative user weight. Every protein is an
   independent present-versus-absent hypothesis with prior probability 0.5 by
   default. Weighted Bayes factors update its odds without normalization across
   proteins. Proteins whose posterior is strictly above the selected node
   cutoff (default 0.5) are retained; curated second messengers can be included
   separately. Candidates outside the validated 891-node seed are inserted
   after incremental edge characterization. Every node dataset also has its
   own `Tq ×` control. An optional global nondetection rule gives each eligible
   but unobserved candidate a configurable BF below 1 in every enabled node
   stream. It is off by default; the conservative default strength is BF 0.5.
   Non-kinases remain out of scope and neutral in the kinase-activity stream.
   The six collecting-duct streams and two selective PKA-subunit-KO streams are
   off in the generic default. The 891-node edge catalog remains the immutable
   seed; the revised site-centric double-KO stream expands the current default
   to 1,051 proteins, or 1,071 nodes with curated messengers.
2. **Characterize edges.** Enable any subset of mpkCCD localization,
   KinasePredictor, STRING, HPA, OmniPath, STITCH, and the optional derived
   scaffold-closure stream. Edge prior odds are multiplied by
   each sparse Bayes factor raised to its selected weight. HPA primary and HPA
   high-confidence are alternatives and cannot be enabled together. Each edge
   dataset has an independent normalization control. An optional global rule
   assigns a configurable BF below 1 to a pair that is eligible for an enabled
   source but has no non-neutral relationship in that source. Out-of-scope
   pairs remain neutral. The option is off by default and its default penalty
   is BF 0.5. Scaffold closure is positive-only and off by default because it
   reuses the selected pre-closure graph rather than adding an independent
   experiment.
3. **Find paths.** Choose a selected start node and any mouse protein target.
   Targets outside the node universe are characterized from the same evidence
   streams and appended as an undirected row/column. Paths are ranked by the
   product of their edge probabilities. The path panel exposes every
   ontology-derived node-role class and the GO root term(s) used to assign it.
   Users select which roles may serve as internal signal relays. The
   sensitivity-oriented default enables ligand, binding, generic
   signaling-process, catalytic, regulatory, and second-messenger classes;
   `adaptor_scaffold` is available but disabled. Multi-role nodes qualify when
   any selected class matches. The optional scaffold override still removes
   every scaffold-tagged intermediate, even if another selected role matches.

   Conservative directionality is enabled by default and combines ontology-role
   rules with mapped OmniPath source-target records. Established ontology
   restrictions are applied first. OmniPath can add a unique direction only for
   a pair ontology left unresolved; it cannot reopen or reverse an
   ontology-disallowed traversal. The original edge probability is retained in
   every allowed direction. Pairs with no direction evidence, and unresolved
   pairs with bidirectional evidence, remain traversable both ways. Source
   disagreements are recorded in the audit. OmniPath direction use has a
   separate checkbox beneath the master directionality control.

   Completed path runs also display a merged pathway map. The browser overlays
   up to the first 50 ranked paths, collapses repeated relationships to one
   edge, and lets the viewer switch among the top 5, 10, 25, or all available
   visualized paths. Node radius and color intensity encode distance of the node
   posterior from the neutral 0.5 prior; green indicates support, red indicates
   counterevidence, and orange identifies curated or external nodes without a
   Bayesian node posterior. Edge width and color use the same encoding for edge
   posteriors. Magnitude is drawn on a capped log-odds scale so relationships
   remain visually distinguishable even when several posteriors are very close
   to one. Arrowheads appear only when the propagation matrix actually
   disallows the reverse traversal. Selecting a mark reveals the exact
   probability, ontology roles, direction status, and contributing path ranks.

4. **Validate temporal order (optional).** Compare the current path order with
   the raw 1/2/5/15-minute dDAVP phosphoproteomic time course. Replicate mode
   stabilizes the three-replicate variance with an empirical intensity trend,
   selects one phosphosite per gene by moderated signal-to-noise, applies a
   within-gene Bonferroni correction by default, and propagates response-time
   uncertainty with Monte-Carlo draws. The panel exposes prior df, alpha, draw
   count, seed, correction method, and the minimum number of scored nodes. This
   stage is off by default and never changes the Bayesian path probability or
   primary rank.

The path count can be set from 1 to 500. All requested paths are shown in the
scrollable result table and written to `ranked_paths.tsv`; the separate
`ranked_path_edges.tsv` retains every constituent edge. The liberal intermediate
policy retains a substantially denser search graph than the earlier conservative
policy, so large path requests can take considerably longer. Runtime varies with
graph size, edge cutoff, hop limit, and endpoint connectivity.

Each run writes a new immutable folder under `results/gui_runs/` containing the
configuration, node posterior table, selected-node universe, adjacency matrix,
supported edge list, path tables, eligibility audit, compact
`top_path_network.json` union graph, and JSON summary.
Directional runs additionally contain `propagation_adjacency_matrix.tsv`,
`ontology_directionality_audit.tsv.gz`, `ontology_class_pair_catalog.tsv`, and
`ontology_direction_rules.json`. Runs using OmniPath directions also contain
`omnipath_direction_evidence.tsv.gz`, including mapped directions, record counts,
resources, references, and consensus-direction counts.
Temporal runs additionally contain `ranked_paths_temporal.tsv`,
`temporal_gene_responses.tsv.gz`, `temporal_variance_trend.tsv`, and
`temporal_validation_summary.json`. The unmodified `ranked_paths.tsv` remains
the canonical Bayesian ranking.
`configuration.json` and `analysis_summary.json` retain the exact selected
ontology-class IDs, while `intermediate_node_eligibility.tsv` reports every
node's complete class list, matching selected roles, endpoint exemption, and
final eligibility decision.

## Per-dataset calibration preferences

Every primary node and edge evidence card exposes a separate `Preferred Tq ×`
control when positive-control calibration is enabled. STRING and STITCH show
the same concept as `Preferred Ref ×`, matching their ordinary `Ref ×` controls.
This preferred value is the center of that stream's regularization penalty; it
does not change an ordinary uncalibrated run and is distinct from the submitted
`Tq ×`/`Ref ×` starting value. Ordinary streams default to 1.0. The five
designated phosphoproteomic streams default to 0.1. Each value is preserved in
`configuration.json` and reported as the `preferred` entry in the calibrated
parameter audit. Derived scaffold closure has no Tq/reference parameter and is
therefore not given a preferred-scale control.

Completed calibrated runs also show the fitted values directly in the result
panel. Separate node and edge tables place each stream's starting and fitted
weight beside its starting, preferred, and fitted Tq/reference multiplier. The
JSON, TSV parameter table, optimizer trace, and bounds remain downloadable for
machine-readable audit, but reading those files is no longer required to see
the fitted model.

## Per-dataset normalization

`Tq × = 1.00` reproduces the finalized analysis for any enabled stream. A value below 1 lowers a
dataset's threshold and increases sensitivity; a value above 1 raises the
threshold and is more stringent. Scalar thresholds (protein, RNA, kinase
activity, all three protein-level phosphoprotein-response streams, the six
collecting-duct abundance streams, and OmniPath) and observation-specific
thresholds (mpkCCD localization, KinasePredictor, and HPA) are all recomputed from their
stored raw statistics. The multiplier is applied to every threshold belonging
to that dataset before its likelihood and Bayes factor are recalculated.

STRING and STITCH are the exceptions: their integrations use score odds rather
than the Gaussian complement kernel. STRING uses the odds of the
STRING combined score relative to a reference score of 0.041, not the Gaussian
complement kernel. STITCH messenger–protein associations use score odds relative
to the database's 0.150 reporting floor. Their cards are labeled `Ref ×`; each
control scales that stream's reference score rather than calling it `Tq`.

Scaffold-mediated closure is calculated in a separate second pass. A protein is
an anchor to an `adaptor_scaffold` node only when its pre-closure edge
probability exceeds the card's `Anchor >` setting (default 0.90). If two
proteins have qualifying edges to at least one shared scaffold, the pair gets
the fixed `Support L` likelihood (default 0.90), equivalent to BF 1.8 relative
to neutral likelihood 0.50. All other pairs remain neutral. Scaffold degree,
anchor strength above the cutoff, and the number of shared scaffolds do not
alter the factor. No empirical Tq is calculated. Closure-derived edges are
never fed back as new anchors. Every qualifying pair, its complete supporting-
scaffold list, and its before/after probability are written to
`scaffold_triadic_closure_audit.tsv.gz` in that run's output folder.

Node factors are recalculated for all 9,170 candidates. The validated 891-node
graph is retained as an immutable seed, while newly selected proteins are
characterized against every selected node and inserted into the graph. Raw
statistics for each new unordered pair are stored in
`data/edge_characterization/incremental_edge_cache/edge_pair_cache.sqlite3`.
Later runs reuse those rows and only characterize pairs that have never been
seen under the current evidence-source signature. Changes to weights or
normalization controls rescore cached raw statistics without rereading the
large source datasets.

New node-specific localization thresholds use the fixed 891-node background,
so adding more nodes does not silently alter an earlier pair. Dynamic
KinasePredictor evidence uses global top-ten model rank for the same reason.
The run summary reports requested pairs, cache hits, newly characterized pairs,
source signature, and database totals.

The collecting-duct extension derives each segment's q75 from positive values
only. Zeros are nondetections: they stay neutral when the optional nondetection
rule is off and receive its configured BF when it is on. Rat proteome rows are mapped to
mouse genes with Ensembl release 116 orthology (UniProt first, rat-symbol
fallback), and every mapping decision is audited. Enabling all six new streams
in addition to the existing four selects 3,330 proteins, or 3,350 nodes after
the 20 curated second messengers. The three segment streams within each assay
share a source experiment and are tagged as dependent even though they remain
separately selectable as requested.

## Posterior distribution plots

Each completed run displays two compact probability histograms: one for all
9,170 modeled protein-node hypotheses and one for every unique undirected edge
hypothesis in the selected graph. The engine computes 50 fixed-width bins over
the full [0, 1] probability range and sends only the bin counts to the browser,
so even multi-million-pair runs do not inflate the job-status response. Bar
height is log10(count + 1), which keeps both the neutral-prior spike and the
smaller supported tails visible. The plot also shows the configured output
cutoff and reports the below-prior, exact-prior, above-cutoff, mean, and range
statistics. The below-prior count makes negative node evidence directly visible.
The same distribution summaries are retained in `analysis_summary.json`.

## Per-hypothesis evidence inspector

After a workflow completes, the result panel can reconstruct the update ledger
for one node or one unordered edge from that exact run. Enter a gene symbol in
the node form, or enter the two endpoint symbols in either order in the edge
form. The inspector reports the stored prior, posterior, cutoff decision, and
one row for every registered evidence stream. Enabled rows are classified as
`supports` when their applied Bayes factor is above 1, `refutes` when it is
below 1, and `neutral` when it equals 1. Disabled streams remain visible for
auditability but contribute zero log-odds.

Each row exposes the applied Bayes factor, fitted/configured weight, and its
weighted change in log2 posterior odds. It also reports where that applied BF
falls in the stream's full run-specific distribution. Node distributions cover
all modeled protein candidates; edge distributions cover every unique
undirected graph pair, including neutral and negative-evidence assignments.
Selecting a row opens a log2(BF) histogram centered on BF=1, marks the selected
hypothesis, and reports its empirical percentile or a compact histogram-based
estimate when the full edge vector was intentionally not serialized. Selecting
a node or edge in the merged network automatically loads the same ledger;
selecting a ranked path highlights its route and initially loads the path's
lowest-probability edge. Selecting a row also adds the raw statistic, source
factor, normalization reference or Tq multiplier, missing-data scope, and
shared-source note when available. The summary independently reconstructs
the posterior from the displayed contributions and compares it with the value
stored in `node_posteriors.tsv.gz` or `edge_adjacency_matrix.tsv`. This
arithmetic check makes stale or incomplete audit data visible instead of
silently presenting a plausible-looking decomposition. Edge inspection is
performed for only the requested unordered pair, so it does not rerun the
full pairwise workflow.

The merged network uses fixed node sizes and fixed edge widths. Color carries
integrated posterior evidence strength on a continuous counterevidence–neutral–
support scale: stronger red is farther below 0.5, gray is near 0.5, and stronger
green is farther above 0.5. Orange still identifies curated messengers or
external endpoints without a Bayesian node posterior. This graph color is the
combined posterior, whereas the inspector's `Weight` column is the distinct
per-stream exponent used in `BF^weight` during updating.

## Temporal validation interpretation

`temporal_soft_precedence` averages the probability that an earlier path node's
sampled response time precedes a later node's. Draw-wise Kendall tau-a is
reported as a mean and 2.5/97.5 percentile Monte-Carlo interval. A path receives
the separate temporal evidence rank only if it reaches the configured measured-
node threshold (default three). Fewer than two scored nodes yields no temporal
order statistic. Insufficient coverage is not treated as evidence against a
path. See `docs/temporal_path_validation.md` for equations, the input schema,
the initial coverage result, and reporting cautions.

## Extending evidence streams

The GUI is driven by `evidence_registry.json`. A new precomputed stream needs:

- a unique ID, label, and description;
- its factor-table path and factor column;
- its neutral value for node evidence, or target-factor column for edge
  evidence;
- default enabled state and weight; and
- normalization handler, label, reference, help text, and default multiplier;
- an optional exclusivity or dependence group.
- for node streams, an observed-status column and eligibility policy when
  nondetection should be usable as negative evidence.

The interface renders registry entries automatically. The engine's combination
logic is generic for positive factors and optional BFs below 1. A stream requiring a new raw-data
transformation should first receive an audited preprocessing module that emits
the standard factor table.

## Scientific safeguards

- Missing sparse edge evidence is neutral (`BF = 1`) by default. When the
  optional edge-negative rule is enabled, an eligible unsupported pair receives
  the configured `BF < 1`; a pair outside the source's assay/database scope
  remains neutral.
- Eligible node nondetections can optionally supply a configurable `BF < 1`;
  out-of-scope nodes remain neutral and the reproducibility default is off.
- Node posteriors are independent Bernoulli probabilities and are not
  normalized across the 9,170 protein candidates.
- Protein/RNA source scores use 0.5 as their neutral likelihood floor and are
  divided by 0.5 before integration, making neutral evidence BF = 1.
- Edge probabilities are independent Bernoulli updates from the selected prior.
- Kinase activity and phosphoprotein response share a source and are visibly
  marked as dependent.
- The selective PKA-Cα- and PKA-Cβ-KO response streams are separately
  selectable but share a dependence-group warning because they come from the
  same TMT experiment.
- HPA alternatives are mutually exclusive.
- Path products are ranking scores, not calibrated whole-path probabilities.
- Every run preserves its exact configuration and complete audit tables.
- Scaffold closure is labeled as dependent proximity/co-complex evidence, not
  proof of a direct binary PPI, and is disabled by default.
- Every active stream records its normalization or fixed-rule parameters in
  the run summary.
