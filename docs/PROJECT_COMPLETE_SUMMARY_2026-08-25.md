# Complete project summary: Bayesian Graphical Inference of renal signaling

Date: 2026-08-25  
Repository: `jjsong23/Bayesian-Graphical-Inference`  
Primary application entry point: `python launch.py`

## 1. Project objective

This project is a hypothesis-generation framework for reconstructing signaling
networks in renal collecting-duct principal cells. It was motivated by the
Bayesian kinase-to-phosphosite mapping strategy described by Deshpande and
colleagues, but it has been expanded beyond kinases to include all candidate
signaling participants and multiple types of relationships.

The framework answers three related questions:

1. Which proteins or small signaling molecules should be considered plausible
   participants in the biological system? This is the **node-selection** stage.
2. Which unordered pairs of selected nodes are plausibly associated? This is
   the **edge-characterization** stage.
3. Which signal-propagating paths best connect a user-selected starting node to
   a target such as Aqp2? This is the **path-inference** stage.

The primary graph remains probabilistic and undirected. A separate ontology
layer can partially orient edge traversal for path finding without changing the
underlying undirected edge-existence probability. An optional temporal layer
can then annotate whether measured phosphoproteomic response times are
consistent with a proposed path order.

The software is intended to make each scientific assumption selectable and
auditable. It is not a clinically validated model, and its posterior values
must not be interpreted as empirically calibrated frequencies of biological
truth.

## 2. Overall architecture

The current workflow is:

```text
9,170 candidate mouse proteins
        |
        | selectable node Bayes factors
        v
selected proteins + 20 curated secondary messengers
        |
        | selectable edge Bayes factors on unique unordered pairs
        v
symmetric undirected probability matrix
        |
        | optional external-target row and column
        | optional ontology traversal restrictions
        v
ranked loopless signal-propagating paths
        |
        | optional time-course annotation
        v
temporally annotated path hypotheses
```

The immutable, extensively audited edge seed contains **891 nodes**: 871
proteins and 20 curated secondary messengers. It represents 396,495 unique
unordered pairs. It is no longer an implementation limit. Node selection can
choose proteins outside this seed; the engine characterizes all new unordered
pairs incident to those nodes and stores their raw evidence in an incremental
SQLite cache.

Under the current generic default node settings, 1,051 proteins exceed the 0.5
node posterior threshold. Adding the 20 curated messengers gives 1,071 graph
nodes, so a default run can extend beyond the seed. A lab-specific profile that
enables all CCD, OMCD, and IMCD proteome/RNA streams selects 3,330 proteins, or
3,350 nodes after adding the messengers.

## 3. Bayesian model

### 3.1 Binary hypotheses and priors

Every node and every unique unordered edge is treated as an independent binary
hypothesis:

- the node or edge is present; or
- the node or edge is not present.

The objective starting prior is 0.5 for presence and 0.5 for absence. Therefore
the prior odds are one:

```text
prior_odds = prior_probability / (1 - prior_probability) = 1.
```

There is no requirement that all node probabilities or all edge probabilities
sum to one. This is not a competition in which support for one candidate
removes probability from another candidate.

### 3.2 Evidence integration

For evidence stream `s`, Bayes factor `BF_is`, and scientist-selected weight
`w_s`, hypothesis `i` is updated as:

```text
posterior_odds_i = prior_odds_i * product_s(BF_is ** w_s)
posterior_i = posterior_odds_i / (1 + posterior_odds_i)
```

Equivalently, the implementation adds weighted log Bayes factors to prior log
odds. A weight of zero removes the stream; weight one uses it as stored; weight
two counts its log evidence twice. Weights are not normalized to sum to one.
The node and edge stages now use the same odds-space update rule.

### 3.3 Quantitative evidence and Tq

Many quantitative streams use the Gaussian-complement kernel:

```text
L(x, Tq) = 1 - exp(-0.5 * (x / Tq)^2)
```

`x` is a nonnegative source-specific score, such as abundance, absolute LFC,
localization similarity, or kinase-model score. `Tq` is an empirical scale,
usually the 75th percentile of the appropriate background score distribution.
The GUI multiplier changes the effective threshold to:

```text
effective_Tq = stored_Tq * user_multiplier.
```

When the historical positive-or-neutral mode is used, likelihood is floored at
0.5 and divided by the neutral likelihood 0.5:

```text
BF = max(0.5, L(x, Tq)) / 0.5.
```

Consequently, a single observation gives BF between one and two. Sparse sources
can yield larger pair-level factors when several independent site contributions
are multiplied for the same pair.

Tq is necessary because raw evidence streams have different numerical scales.
It makes `x/Tq` dimensionless and defines where a source transitions from weak
to strong evidence. It is not a probability normalizer and does not force
posteriors to sum to one.

STRING and STITCH are exceptions. Their confidence scores are converted to odds
and compared with a configurable reference-score odds value. Their GUI control
is therefore a **reference multiplier**, not a Gaussian Tq multiplier.

### 3.4 Negative evidence

Missing or weak evidence is neutral by default. Two optional negative-evidence
models are available.

The first is a simple fixed absence factor. When enabled, a source-eligible
nondetection or unsupported pair receives a configurable BF below one, default
0.5. Repeated negative factors multiply the odds downward. Hypotheses outside
the source's biological and technical scope remain neutral.

The second is continuous negative evidence. It removes the historical 0.5
likelihood floor:

```text
BF = max(epsilon, L(x, Tq) / 0.5).
```

Eligible nondetections are represented as `x=0`, which maps to a small positive
BF floor rather than exact zero. Measured low scores can therefore produce
graded factors below one. The positive stored Tq is retained; a background
quantile is never recomputed after zero-padding all 9,170 candidate proteins.

For nodes, kinase-activity evidence applies only to kinases, whereas the
abundance, transcript, phosphoprotein, and collecting-duct streams treat all
protein candidates as eligible. Secondary messengers are never penalized by
protein-only datasets.

For edges, eligibility is source-specific:

- mpkCCD localization requires profiles for both endpoints;
- HPA localization requires HPA profiles for both endpoints;
- KinasePredictor requires a kinase endpoint and an observed scorable
  phosphosite on the protein endpoint;
- STRING requires two mapped protein endpoints;
- OmniPath applies to protein-protein pairs;
- STITCH applies to curated-messenger/protein pairs with the required mapping;
- scaffold closure has no negative-absence rule.

Database absence is thus penalized only when the pair was technically and
biologically testable under that source. The option is a sensitivity assumption
and must not be described as proof that the relationship is absent.

## 4. Candidate signaling-node universe

The initial universe was assembled from a liberal Gene Ontology-based signaling
screen. The canonical identifier for downstream analysis is the mouse gene
symbol in the `symbol` column. Duplicate symbols were removed and the source
GO-term catalog was retained for audit. The consolidated backend candidate
catalog contains 9,170 proteins.

The liberal philosophy prioritizes sensitivity. Broad signaling participants,
ligands, receptors, enzymes, adaptors, binding proteins, and regulators are
retained rather than imposing a narrow kinase-centric definition. Path
propagation uses a separate selectable ontology policy so that a broad node
universe does not automatically mean every selected protein can act as an
intermediate.

Twenty curated secondary messengers are added after protein node selection.
They are not gene products and are not scored by protein abundance, RNA, or
protein-localization nondetection rules. Their relationship evidence is supplied
through STITCH and curated biochemical rules where available.

## 5. Node-selection evidence streams

There are currently 12 quantitative node streams.

### 5.1 mpkCCD protein abundance

The source is a protein-abundance dataset from mpkCCD cells. Linear relative
abundance is used as the raw score. All finite source abundances define the
background and its q75 Tq before filtering to the signaling candidates. This
stream originally yielded fewer non-neutral candidates than RNA because its
empirical distribution and measurement coverage were substantially narrower.

### 5.2 Principal-cell RNA abundance

Only principal cells from the single-cell dataset are used. The current score
is raw median TPM rather than log10 TPM. Zero values are treated as
nondetections and are excluded from the positive background used to calculate
Tq. Positive finite principal-cell values define the empirical background.

### 5.3 PKA-double-knockout kinase activity

The PKA-knockout phosphoproteome is scored with KinasePredictor. For every
observed phosphosite, the top ten ranked kinase models are retained as hits to
increase sensitivity. The kinase-level raw statistic is based on absolute
phosphosite log2 fold change; no protein-abundance normalization or
phosphorylation-specific normalization is performed. A kinase not present in
the candidate kinase mapping remains neutral.

### 5.4 PKA-double-knockout phosphoprotein response

This stream asks whether a protein has at least one phosphosite that changes in
the PKA knockout. The raw protein score is the strongest absolute phosphosite
LFC. Proteins without a detectable phosphosite remain neutral unless a negative
evidence option is deliberately enabled.

There is **no correction for the number of observed sites on a protein**. The
accepted biological rule is site-centric: a protein is a hit if at least one
particular phosphosite shows a sufficiently strong change. Sites are not
accumulated and proteins are not penalized or advantaged according to how many
sites happened to be detected.

### 5.5 PKA-C-alpha and PKA-C-beta knockout phosphoprotein responses

The separated PKA catalytic-subunit dataset supplies two optional streams:

- PKA-Cα knockout phosphoprotein response; and
- PKA-Cβ knockout phosphoprotein response.

They are kept separate because the two catalytic subunits can have opposing or
nonredundant downstream effects that a double knockout may hide. Within each
comparison, duplicate rows are collapsed by median signed LFC for each
UniProt/site-pattern key. Absolute site LFC is compared with that comparison's
all-source-site q75. A protein inherits its strongest site factor. There is no
site-count correction or site accumulation. Pmod is retained for audit but is
not part of scoring. The two streams share a dependence group and are disabled
in the generic default.

### 5.6 Collecting-duct segment proteomics

Rat proteomics is represented as three independent selectable streams:

- CCD proteome;
- OMCD proteome; and
- IMCD proteome.

Rat identifiers were harmonized to mouse gene symbols using audited Ensembl
release 116 orthology. UniProt mapping is preferred, rat-symbol fallback is
audited, high-confidence orthologs are preferred, valid one-to-many orthologs
are retained, and duplicate contributions are aggregated by maximum abundance.
Each segment uses all positive finite values from that segment to define q75
before candidate filtering.

### 5.7 Collecting-duct segment RNA

Mouse renal-tubule RNA abundance is likewise represented as three streams:

- CCD RNA;
- OMCD RNA; and
- IMCD RNA.

Each segment has its own background and Tq. Zero is a nondetection. The three
RNA streams share a dependence group; the three proteome streams share another
dependence group. All six segment streams are optional so the generic baseline
remains reproducible.

## 6. Edge-characterization evidence streams

Edge characterization is undirected and symmetric. A unique unordered pair is
stored once in pair tables and appears twice only in the square matrix, at
`(i,j)` and `(j,i)`. Self-edges are excluded.

There are seven primary quantitative edge streams plus one derived scaffold
stream.

### 6.1 mpkCCD subcellular localization

Each protein has an abundance profile across the 1K, 4K, 17K, 200K pellet, and
200K supernatant fractions. Pairwise raw evidence is based on the dot product of
the two profiles. Directed threshold evaluations are symmetrized by averaging
the two endpoint-relative likelihoods. Missing profiles are neutral unless the
negative-evidence policy explicitly makes the pair eligible and unsupported.

### 6.2 KinasePredictor kinase-protein support

Potential phosphosites were enumerated from the mouse UniProt reference
proteome, and the final edge analysis is restricted to observed phosphosites.
For each site, the top ten globally ranked kinase models are retained. Only
kinases and protein products in the current node universe are considered. Each
site prediction supplies undirected kinase-protein evidence. Multiple qualifying
sites for the same pair contribute in log-BF space, so their factors multiply.
The graph records only the existence of a kinase-protein relationship, not the
direction of phosphorylation flow.

### 6.3 STRING v12

STRING combined scores are treated as functional-association confidence, not
as proof of a direct physical PPI. Mapped mouse protein scores are converted to
odds and divided by the odds of a configurable reference score. The resulting
factor applies symmetrically. Small molecules are outside protein-only STRING
scope.

### 6.4 Human Protein Atlas localization

Mouse proteins are mapped to human orthologs and HPA subcellular-location
profiles. Localization overlap is represented by cosine similarity. Two
mutually related alternatives are exposed:

- HPA primary localization; and
- HPA high-confidence localization.

They should generally not be enabled together as if independent. Source
profiles, mapping audits, and location-specific q75 values are retained.

### 6.5 OmniPath core

Mouse OmniPath post-translational interactions are collapsed to one unordered
protein pair. Direction, stimulation, inhibition, and source annotations remain
available for audit, but the edge-existence Bayes factor intentionally ignores
direction and sign at this stage. Curation effort is transformed into support
using a Tq-scaled saturating kernel.

### 6.6 STITCH messenger-protein evidence

STITCH supplies protein-chemical association evidence for curated secondary
messengers. Fourteen of the 20 messenger nodes currently have STITCH support;
the other six remain neutral for this stream. Like STRING, STITCH uses score
odds relative to a reference score rather than the Gaussian-complement Tq.

### 6.7 Scaffold-mediated closure

Scaffold closure is a derived, optional, deliberately simple binary rule:

1. Use only the graph before closure.
2. A protein-scaffold anchor qualifies when its probability is strictly above
   the configured threshold, default 0.90.
3. Two proteins qualify when both have such an edge to the same exact
   `adaptor_scaffold` node.
4. The qualifying protein pair receives likelihood 0.90, equivalent to BF 1.8
   against neutral likelihood 0.5.
5. All other pairs receive neutral BF one.
6. Apply closure once; closure-generated edges cannot create new anchors.

There is no scaffold-degree or “promiscuity” penalty, no noisy-OR score, and no
empirical scaffold Tq. Those earlier heuristics were rejected because they were
not biologically justified for this application.

## 7. Colocalization and AlphaFold-pair filtering

A separate colocalization graph was created to reduce the number of protein
pairs requiring structural prediction. It integrates:

- mpkCCD fractionation localization;
- HPA localization; and
- COMPARTMENTS localization.

Scientist-reviewable compatibility matrices encode both exact compartment
identity and plausible adjacency. For example, cytosolic and plasma-membrane
proteins are not identical in location but can interact at the membrane
interface. These `C` matrices are compatibility assumptions, not interaction
probabilities.

The three localization streams are combined in the same Bayesian odds
framework. The resulting matrix is symmetric; unique-pair counts use only one
triangle. In the 848-protein analysis, there were 359,128 possible unique
protein pairs, of which 167,884 had non-neutral adjacency-aware colocalization
support. Thus 191,244 pairs were removed by the strict near-baseline filter.

Known PPI evidence from STRING, BioGRID, and IntAct was assembled into tiers.
Tier A and Tier B relationships were used as an exclusion filter for structural
prediction, not as a new Bayesian edge factor. BioGRID/IntAct absence remained
missing evidence rather than negative evidence.

The final scaffold-focused AlphaFold candidate calculation used the
`adaptor_scaffold` class as a reproducible proxy for scaffold/structural
proteins. Starting from the 167,884 colocalized pairs, Tier A/B known pairs were
removed and at least one scaffold endpoint was required. This produced 66,703
unique unordered pairs:

- 58,265 scaffold-other pairs;
- 8,438 scaffold-scaffold pairs; and
- 80 pairs with ambiguous-only prior evidence retained.

Validation found no reversed duplicates, self-pairs, missing expected pairs, or
remaining Tier A/B overlaps.

Using the supplied empirical runtime model, these 66,703 pairs represented
approximately 27,298 GPU-hours of folding work. The ideal four-GPU wall-time
estimate was approximately 285.8 days for the benchmark configuration, or 58.3
days under the one-prediction-per-model approximation. These are planning
estimates, not guaranteed runtimes.

## 8. Incremental characterization of newly selected nodes

The GUI is not capped at 891 nodes. When node selection adds a protein outside
the seed, the engine computes only unordered pairs incident to newly added
nodes. It reconstructs raw profiles and observations for mpkCCD localization,
KinasePredictor, STRING, HPA, OmniPath, and STITCH, then stores these raw values
in:

`data/edge_characterization/incremental_edge_cache/edge_pair_cache.sqlite3`

Weights and Tq/reference multipliers are not baked into the cache. Cached raw
observations are rescored under the current configuration, so changing a model
parameter does not trigger new database/source scans. A source signature guards
against reuse after evidence inputs change.

The cache is a performance optimization rather than the scientific record. The
shareable application release omits the current 3.35-GB cache because it
contains several historical source signatures and is not portable across ZIP
extraction timestamps. All source files required to rebuild it are included.
The first run that adds nodes can therefore be slower; later runs reuse the new
local cache.

## 9. External targets and path inference

A target such as Aqp2 does not have to be part of the selected signaling-node
universe. The target-extension function computes one new adjacency vector
between the target and every active graph node using the same selectable edge
streams. It appends this vector as a symmetric row and column to the adjacency
matrix. The target can be an endpoint but is not automatically allowed as an
internal signaling relay.

Path search ranks loopless paths using edge probabilities. The primary path
score is the product of the traversed edge probabilities; equivalently, the
negative logarithms of the edge probabilities form additive path costs. The GUI
exposes the start node, target, maximum hops, minimum usable edge probability,
and number of paths to report.

Only selectable signaling-role classes may act as intermediates. Multi-role
kinase/scaffold proteins remain allowed when they have any propagating role.
Only proteins whose permitted role is exclusively scaffold/adaptor are excluded
under the default policy. The user can select which ontology classes qualify as
signal relays.

Whole-path products are ranking scores, not calibrated probabilities that the
complete biological pathway is correct. Paths are hypotheses supported by the
current graph assumptions.

## 10. Ontology-based partial directionality

Edge existence remains undirected. A separate propagation matrix removes a
reverse traversal only when endpoint ontology roles uniquely support one
direction. Unresolved or conflicting pairs remain traversable in both
directions. This conservative fail-open policy avoids inventing causal order.

The current versioned rule catalog includes:

- ligand to receptor;
- receptor regulator to receptor;
- kinase to kinase/phosphatase-binding protein;
- phosphatase to kinase/phosphatase-binding protein;
- GTPase regulator to GTPase;
- cyclase to secondary messenger;
- phospholipase to secondary messenger;
- nitric-oxide synthase to secondary messenger;
- phosphodiesterase to secondary messenger; and
- secondary messenger to second-messenger-binding protein.

All endpoint roles are considered. Rules pointing uniquely one way orient the
traversal. Rules supporting both ways are marked conflicting and retain both
directions. Pairs without a rule also retain both directions.

In the initial 892-node Prkaca-to-Aqp2 validation, 169,769 undirected edges were
above probability 0.5. Of these, 24,679 were uniquely oriented (14.5368%),
144,027 had no matching ontology rule, and 1,063 had conflicting multi-role
rules. Directionality removed only the 24,679 disallowed reverse traversals.

Biochemical activation/inhibition sign remains separate from traversal
direction. OmniPath signs and the phosphodiesterase negative-effect annotation
are retained for audit but do not yet change the path score.

## 11. Optional temporal validation

The temporal module uses the 1-, 2-, 5-, and 15-minute dDAVP/vehicle replicate
phosphoproteome. It is an optional post-path annotation: it does not contribute
a node or edge BF and does not change the canonical Bayesian rank.

For each usable replicate, log2 treated/control response is calculated only
when both reporter intensities are positive. Raw site/timepoint variance is
moderated toward an empirical intensity-dependent prior. Representative sites
are selected by maximum moderated absolute t statistic. The default within-gene
Bonferroni gate adjusts for choosing among multiple sites and timepoints.

Response-time draws use:

```text
tau = sum_t(t * abs(response_t)) / sum_t(abs(response_t)).
```

Pairwise soft precedence estimates `P(tau_i < tau_j)`. Draw-wise Kendall tau-a
produces a conditional Monte-Carlo interval. These are not posterior credible
intervals from a full longitudinal Bayesian model.

The initial default validation loaded 6,755 valid phosphosite trajectories and
assigned representatives to 2,661 genes. Only 237 passed the within-gene gate.
Seven of the 25 paths had at least two scored nodes, but none reached the default
three-node threshold for a temporal evidence rank. This was reported as
insufficient coverage, not negative evidence against those paths.

## 12. Parameter calibration

Node and edge calibration are separate. The user supplies known-present positive
controls for one stage. Bounded Powell optimization adjusts selected evidence
weights and Tq/reference multipliers to increase the probability of those
controls. Unknown hypotheses are not silently labeled absent.

A quadratic regularization term keeps weights and scale multipliers near
scientist-selected preferred values. Every stream exposes its own preferred Tq
or reference multiplier. Phosphoproteomic streams can be anchored near a low
preferred multiplier because they have low sensitivity but high specificity.
Regularization prevents the optimizer from obtaining a trivial fit by pushing
all weights or thresholds to extreme boundaries.

The calibration objective is conditional on the supplied positive controls. It
does not prove that fitted parameters generalize to unknown biology. Held-out
controls remain preferable when enough validated examples are available.

## 13. Sensitivity analyses

The current evidence-information audit evaluates all 12 node streams and all
seven primary edge streams on the fixed 891-node edge seed. Scaffold closure is
excluded from intrinsic-information and Tq rankings because it is a derived
binary rule with no independent score distribution or Tq.

### 13.1 Weight of evidence

For hypothesis `i` and stream `s`, signed weight of evidence is:

```text
w_s * log2(BF_is) bits.
```

Its absolute magnitude measures evidence strength regardless of direction.
This quantity is exactly the change in base-2 log posterior odds supplied by
that stream. It is unbounded and is not the same as uncertainty reduction.

### 13.2 Realized KL information

Starting from prior 0.5, a stream alone yields:

```text
p_is = BF_is ** w_s / (1 + BF_is ** w_s).
```

Binary entropy is:

```text
H(p) = -p*log2(p) - (1-p)*log2(1-p).
```

Because the equal prior has one bit of entropy, realized information gain is:

```text
D_KL(Bernoulli(p) || Bernoulli(0.5)) = 1 - H(p).
```

This is zero at posterior 0.5 and approaches one bit as posterior approaches
zero or one. It measures certainty gained, not whether the evidence is positive
or negative. Total stream information is the sum over eligible hypotheses.

### 13.3 Leave-one-stream-out influence

For every default-enabled stream, the analysis compares the full integrated
posterior with the posterior after removing that stream while holding the prior,
other streams, weights, and scales fixed. It reports posterior shifts, decision
flips across 0.5, Spearman agreement, and conditional KL divergence:

```text
D_KL(Bernoulli(p_full) || Bernoulli(p_without_stream)).
```

### 13.4 Agreement and disagreement

All 66 unique node-stream pairs and all 21 unique edge-stream pairs are tested.
Continuous-negative sensitivity uses BF floor 0.001 and compares signs of
`log2(BF)`. Pairwise tables report joint eligibility, non-neutral overlap, sign
concordance, positive-hit Jaccard overlap, and Spearman correlation. Higher-order
tables record the positive, negative, and non-neutral stream counts for every
hypothesis, sign consensus, and binary sign entropy.

Under the continuous-negative scenario, 3,371 node candidates and 136,615
unique edges had at least one positive and one negative stream. Restricting to
the default-active model yielded 1,034 conflicted nodes and 134,429 conflicted
edges. These values reflect the modeled low/absence assumptions, not proven
biological contradictions.

### 13.5 Tq response curves

Every stream was evaluated at 41 log-spaced scale multipliers over the same
bounds permitted by GUI calibration. The plots report:

- geometric mean effective BF, `geomean(BF ** weight)`, over eligible
  hypotheses; and
- total realized KL information from the 0.5 prior.

All 19 evaluated streams showed monotonic decline in both quantities as the
Tq/reference multiplier increased. All maxima occurred at the smallest allowed
multiplier: 0.05 for the phosphoproteomic/kinase-model streams and 0.25 for the
other streams. There were zero monotonicity violations across 779 curve points.

This is mathematically expected under positive-or-neutral scoring: decreasing
Tq cannot weaken a BF. Therefore maximizing evidence magnitude or realized KL
alone is **not a valid unsupervised Tq optimizer**; it will always prefer a search
boundary. Scientifically defensible scale selection requires external controls,
held-out predictive performance, and regularization.

## 14. Graphical workbench

The local minimal web application exposes the entire workflow end to end. Users
can:

- enable or disable every node and edge evidence stream;
- independently edit each stream's weight and Tq/reference multiplier;
- choose fixed or continuous negative evidence per stream;
- edit stage priors, selection cutoffs, and BF floors;
- enter known-positive controls and calibrate node and edge parameters
  separately;
- edit every stream's preferred calibration Tq/reference multiplier;
- include secondary messengers;
- add nodes beyond the seed and compute only newly required edge pairs;
- enable scaffold-mediated closure;
- choose the start node, target, path count, hop limit, and edge cutoff;
- choose the ontology classes allowed as signaling intermediates;
- enable conservative ontology directionality;
- optionally enable temporal validation;
- inspect node and edge posterior distributions;
- cancel queued or running analyses from the GUI; and
- download each run's explicit output tables and summaries.

The server permits one active analysis at a time. Additional submissions receive
a clear conflict response rather than silently hanging at “Starting.” A
cooperative cancellation flag is checked throughout long-running calculations.

Every run creates a new immutable directory under `results/gui_runs/`. Outputs
include the submitted configuration, node posterior audit, selected universe,
edge probability matrix, supported unique-pair table, ranked paths, ontology
eligibility/directionality audits, plot-ready posterior distributions, and a
machine-readable summary. Temporal runs add response and temporal-ranking
tables.

## 15. Reproducibility and scientific reporting cautions

The project preserves raw inputs where redistribution permits, processed
mapping tables, source-specific factors, pair/site audits, symmetric-matrix
checks, diagonal checks, exact reconstruction checks, and JSON summaries.

The following statements should be avoided:

- Do not describe node or edge posteriors as calibrated biological truth
  frequencies.
- Do not call neutral posterior 0.5 “zero.”
- Do not say missing data proves absence; negative evidence is an optional
  eligibility-scoped model.
- Do not call STRING associations experimentally demonstrated physical PPIs.
- Do not say every kinase-prediction site was observed in the PKA-knockout
  experiment; the phosphosite database includes UniProt-derived annotations.
- Do not claim kinase-activity and PKA phosphoprotein streams are independent;
  they use related source data.
- Do not call compartment compatibility values PPI probabilities.
- Do not equate `adaptor_scaffold` with all structural proteins.
- Do not claim Tier A/B PPI evidence was integrated as another edge BF; it was
  used as an AlphaFold exclusion filter.
- Do not describe the primary edge matrix as directed; only propagation
  traversal is partially oriented.
- Do not describe 891 nodes as a hard limit.
- Do not claim OmniPath direction or sign currently determines the undirected
  edge factor.
- Do not claim temporal annotation updates the Bayesian edge posterior or
  proves causality.
- Do not call temporal percentile intervals fully Bayesian credible intervals.

External-source licenses and citation requirements must be checked before
redistributing data outside controlled lab storage. In particular, a complete
local runtime bundle is different from a public data release.

## 16. Important files

- `launch.py`: portable application launcher and release validator.
- `gui/server.py`: local HTTP server and cancellation/job management.
- `gui/workflow_engine.py`: canonical node, edge, target, path, and output logic.
- `gui/evidence_registry.json`: declarative evidence-stream registry and
  defaults.
- `gui/parameter_calibration.py`: bounded positive-control optimizer.
- `code/bayes_factors.py`: reusable Bayesian transformations.
- `code/edge_characterization/incremental_edge_cache.py`: persistent raw
  pair-evidence cache.
- `code/path_finding/`: target extension, path search, ontology directionality,
  and temporal ranking.
- `data/node_selection/node_universe_combined_nonzero.tsv`: validated 891-node
  seed universe.
- `results/combined_kinase_phosphoprotein_evidence/combined_all_nodes.tsv`:
  complete node-factor catalog.
- `results/backend_bayes_factor_catalogs/edge_factors_891/`: seed edge-factor
  tables.
- `docs/full_project_methods.md`: chronological full methods record.
- `CODEX_PROJECT_CONTEXT.md`: compact handoff for another coding session.

## 17. Running the shareable application

After extracting the application ZIP:

```text
python -m pip install -r requirements.txt
python launch.py --check
python launch.py
```

The default server binds only to `127.0.0.1` and opens
`http://127.0.0.1:8765/`. Use `python launch.py --foreground` to keep it attached
to the terminal and stop it with Ctrl+C. The application does not require an
external web server or a JavaScript build step.

The release contains the code, documentation, seed node/edge factors, raw
quantitative inputs needed to change Tq/reference multipliers, and source files
needed to characterize newly selected nodes. It intentionally excludes Git
history, Python environments, old GUI-run output directories, and the large
disposable incremental SQLite cache.

## 18. Current development status and logical next steps

The end-to-end system is operational. Major remaining scientific tasks include:

1. constructing larger, independently validated positive and negative control
   sets for parameter calibration;
2. validating or revising every ontology propagation rule with domain experts;
3. deciding whether directed/sign-specific OmniPath evidence should become a
   separate causal layer;
4. evaluating dependence-aware integration so correlated datasets are not
   treated as fully independent;
5. improving temporal coverage and cross-species mapping before using time
   order to rank paths; and
6. prospectively testing high-ranking, previously unsupported edges and paths.

Until those validations are complete, the graph and ranked paths should be
reported as auditable probabilistic hypotheses rather than definitive signaling
mechanisms.

## 19. Saved evidence-score and cutoff distributions

The release now includes a complete visual audit of the quantitative inputs to
all 12 node evidence streams and all 8 edge evidence streams. The figures are
under
`results/sensitivity_analysis/evidence_information_2026-08-24/score_distributions/`.
There is one consolidated node figure, one consolidated edge figure, and one
PNG per individual evidence stream. A TSV beside the plots records the number
of scores, number of zeros, score quantiles, cutoff quantiles, population
definition, and an inference-specific note.

The blue histograms show raw quantities entering the evidence transformations.
For streams supporting continuous negative evidence, all eligible hypotheses
are represented and eligible nondetections are entered as zero. A red vertical
line shows a fixed default Tq/reference value. When the inference uses a
node- or phosphosite-specific threshold, the line is the median threshold and a
red band spans its 10th to 90th percentile. Highly skewed nonnegative scores are
displayed on `log10(1 + x)` axes without dropping observations.

Important source-specific distinctions remain explicit:

- STRING and STITCH use odds reference scores, not Gaussian-kernel Tq values.
- Scaffold closure has no empirical Tq; its plot shows the pre-closure edge
  probabilities and the fixed 0.90 anchor cutoff.
- KinasePredictor inference multiplies retained site-level factors. Its
  histogram uses the maximum raw motif score per eligible pair only as a
  readable display statistic and does not replace the per-site calculation.
- mpkCCD and HPA localization use endpoint-specific thresholds followed by
  reciprocal symmetrization; their red bands summarize those endpoint Tq
  populations.

The reproducible renderer is
`code/sensitivity_analysis/render_evidence_score_distributions.py`.

## 20. Shareable release artifacts

Two local release archives are generated under `deliverables/`, which is
intentionally ignored by Git:

- `BGI_Plots_2026-08-25.zip` contains every project-generated PNG using short,
  numbered archive names, plus a manifest that maps each file back to its full
  project-relative path and records dimensions, byte sizes, and SHA-256 hashes.
- `BGI_App_2026-08-25.zip` contains the executable
  workbench, analysis code, tests, documentation, the full node factor catalog,
  all 891-node seed edge catalogs, source values needed for interactive
  Tq/reference rescoring, and raw sources needed for incremental characterization
  of nodes added beyond the seed.

The application archive intentionally omits the 3.35 GB disposable SQLite
incremental-edge cache, historical GUI runs, old outputs, logs, PID files,
browser dependencies, and duplicate exploratory workbooks. A fresh SQLite
cache is created automatically and retains newly computed raw pair evidence.
This omission both reduces package size and avoids treating a cache with a
stale extraction-time evidence signature as primary scientific data.

The application archive was extracted into a clean directory and validated
with `python launch.py --check`. It reported 12 node streams, 8 edge streams,
9,170 node candidates, 4 default node streams, and 6 default edge streams. A
fresh server was then started from the extracted copy; `/api/health` returned
`{"status": "ok"}` and the main page returned HTTP 200. Archive CRC checks and
SHA-256 checksums are generated automatically. The reproducible builder is
`code/release/build_shareable_packages.py`; detailed packaging notes are in
`docs/SHAREABLE_RELEASE_2026-08-25.md`.

The release builder uses short archive and internal root names and flattens the
plot PNGs to unique numbered filenames. The plot manifest retains their full
source provenance. It enforces a maximum hypothetical extraction path of 240
characters when extracting beside the ZIP in `deliverables/`, leaving margin
below the Windows Explorer 260-character legacy limit.
