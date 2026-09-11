# External-target adjacency-vector extension

## Purpose

Path inference may end at a protein that was intentionally excluded from the
signaling-node universe. Aquaporin-2 (`Aqp2`) is the motivating example. The
target-extension function treats such a protein as one additional graph node
and estimates an edge probability between it and every node in the active
universe. It does not rerun or alter any relationship among existing nodes.

The implementation is:

`code/path_finding/build_target_adjacency_vector.py`

Its public entry point is:

```python
build_target_adjacency_vector(
    target,
    *,
    project_root=PROJECT_DEFAULT,
    output_dir=None,
    write_outputs=True,
    write_extended_matrix=True,
    top_n=10,
)
```

`target` can be a mouse gene symbol or a UniProt accession. The target is
resolved against the cached reviewed mouse reference proteome. The function
rejects a target already present in the active universe, because that node's
column already exists in the graph.

## Fixed graph and probabilistic update

Suppose the active universe contains (N) nodes and its current matrix is
(A\in[0,1]^{N\times N}). For an external target (t), the function estimates
the vector

\[
v_t = (P(E_{t1}), P(E_{t2}), \ldots, P(E_{tN})),
\]

where each entry is a separate Bernoulli edge hypothesis. Every hypothesis
starts at prior probability 0.5. For evidence stream (s), a Bayes factor
(BF_{s,ti}) updates the edge odds:

\[
O^{(s)}_{ti} = O^{(s-1)}_{ti}BF_{s,ti},
\qquad
P^{(s)}_{ti}=\frac{O^{(s)}_{ti}}{1+O^{(s)}_{ti}}.
\]

Missing coverage receives (BF=1), so it leaves the current probability
exactly unchanged. The streams are applied in the same order as the current
edge-characterization workflow:

1. mpkCCD differential-fraction localization;
2. observed/curated-site KinasePredictor;
3. STRING v12 functional association;
4. HPA v25.1 primary subcellular localization; and
5. OmniPath core signaling interaction evidence.

The enlarged matrix is then

\[
A'=
\begin{bmatrix}
A & v_t^T\\
v_t & 0
\end{bmatrix}.
\]

This construction makes the added graph unsigned and symmetric. The target
self-edge is zero. The existing (N\times N) block is copied, not recomputed.

## Evidence-stream construction

### mpkCCD localization

The target is looked up in the same subcellular-fractionation workbook used for
the universe. Replicate/source rows are collapsed with the existing pipeline's
rules. Dot products between the target's five-fraction profile and every
covered universe profile are calculated. The target-specific (T_q) is the
75th percentile of its dot products against all covered current nodes. Existing
node-specific (T_q) values are reused from the processed localization table.
The two directed complement-of-minimum likelihoods are averaged to make one
undirected likelihood, then divided by the neutral likelihood 0.5 to obtain the
Bayes factor.

### KinasePredictor

The target-site database is the union of UniProt-annotated target phosphosites
and sites observed for the target in the PKA-knockout phosphoproteomic audit.
Canonical 13-residue sequence windows are preferred; the experimental
dataset's centralized 13-mer is used as a fallback when needed. Every scorable
site is evaluated with the official KinasePredictor matrices. The top ten
distinct predicted kinase genes that are also nodes in the active universe are
retained per site. Each retained kinase-site match contributes a site-specific
Bayes factor, and log Bayes factors are summed when several target sites support
the same kinase-target pair. Non-kinase target-node pairs are neutral for this
stream.

### STRING

The target is mapped to a mouse STRING protein through the same audited
UniProt/preferred-name procedure used for current nodes. All incident records
in the detailed STRING v12 network are scanned, restricted to current-universe
partners, and duplicate target-node records are collapsed by their maximum
combined score. The combined STRING confidence is converted to odds and divided
by the workflow's neutral STRING odds at score 0.041. Absence of a mapped
association is neutral; it is not treated as evidence against an edge.

### HPA localization

The mouse target is mapped to its one-to-one human ortholog through the cached
MGI/HGNC mapping. If the ortholog has HPA Enhanced, Supported, or Approved
primary localization calls, those calls form a binary localization profile.
Cosine overlaps, target- and node-specific 75th-percentile thresholds,
complement-of-minimum scoring, and reciprocal averaging follow the current HPA
edge workflow. A missing HPA subcellular profile is neutral.

### OmniPath

OmniPath core rows incident to the target are restricted to partners in the
active universe and collapsed to unique unordered target-node pairs. Curation
effort is the number of distinct supporting references. The current global
threshold (T_q=6) and complement-based likelihood transform are reused.
Original direction, stimulation/inhibition flags, resources, and references are
retained in audit columns, but direction and sign do not affect the factor.

## Output files

The default output directory is
`results/path_finding/target_extensions/<resolved_symbol>/` and contains:

- `target_adjacency_vector.tsv`: complete canonical-order vector plus all
  stream-level evidence, Bayes factors, and intermediate posteriors;
- `target_adjacency_vector_ranked.tsv`: the same records ordered by final edge
  probability;
- `extended_adjacency_matrix.tsv`: original matrix plus the symmetric target
  row and column;
- `target_phosphosites.tsv`: target-site provenance and sequence windows;
- `target_kinase_predictions.tsv`: retained site-level KinasePredictor hits;
- `analysis_summary.json`: parameters, evidence coverage, validation results,
  and SHA-256 hashes; and
- `README.md`: output-specific interpretation notes.

## Aqp2 validation run

The reference invocation was:

```powershell
python code/path_finding/build_target_adjacency_vector.py Aqp2
```

`Aqp2` resolved to reviewed mouse UniProt accession `P56402`. The active graph
contained 891 nodes, so the function produced 891 target-node probabilities and
a 892 by 892 matrix. Of the target-node edges, 357 were above the exact 0.5
baseline and 534 remained exactly at baseline. Stream-level non-neutral counts
were 243 for mpkCCD localization, 25 for KinasePredictor, 158 for STRING, zero
for HPA, and three for OmniPath. Counts overlap and therefore must not be added
to obtain 357.

The HPA result was neutral because the one-to-one human ortholog `AQP2` was not
present in the downloaded HPA subcellular-location table. This was treated as
missing coverage, not negative localization evidence. Three target
phosphosites were scorable, producing 30 retained site-level kinase predictions
and 25 unique kinase-node hits. The maximum final target-edge probability was
0.998852, for `Prkacb`.

Independent validation confirmed:

- vector length and node ordering matched the universe;
- every probability was between zero and one;
- every pair lacking non-neutral evidence remained exactly 0.5;
- the original 891 by 891 matrix block had maximum absolute difference zero;
- the 892 by 892 matrix was exactly symmetric; and
- the complete diagonal was zero.

## Path-finding handoff

Because an adjacency entry is a probability-like support value rather than a
distance, a shortest-path routine should first convert retained edges to a
nonnegative cost, for example

\[
c_{ij}=-\log(p_{ij}).
\]

An edge-probability cutoff and/or maximum hop count should be selected before
path inference. Without a cutoff, the many 0.5 neutral edges would make the
graph artificially dense and could create biologically uninformative paths.

# Ranked path-finding implementation

## Search space and score

`code/path_finding/find_ranked_paths.py` implements the next stage. The public
function is `find_ranked_paths(start, target, ...)`. The start must resolve to a
node in the active signaling graph. The target may already be in the graph or
may be an external mouse protein. In the latter case, the target-adjacency
extension described above is generated automatically; a cached extension is
used only after its input and output hashes and validation flags are checked.

The default search excludes every edge at the exact 0.5 neutral baseline,
allows at most six edges per path, restricts internal nodes to mechanistic
signal relays, and returns up to 25 loopless paths. For a path
`r = (v_0, v_1, ..., v_h)`, the ranking score is

\[
S(r)=\prod_{i=0}^{h-1}p(v_i,v_{i+1}).
\]

Because logarithms are monotonic, descending (S(r)) is equivalent to
ascending

\[
C(r)=\sum_{i=0}^{h-1}-\log p(v_i,v_{i+1}).
\]

The additive representation permits conventional shortest-path machinery.
The implementation uses Yen's algorithm to rank distinct simple paths. Each
spur-path subproblem is solved by Dijkstra search on a hop-layered graph, which
enforces the maximum-hop setting. No NetworkX dependency is required.

The product naturally penalizes additional edges but does not prohibit a
multi-step path from outranking a direct edge when every edge in the multi-step
path has sufficiently high support. The output also reports hop count,
geometric mean edge probability, and minimum/bottleneck edge probability so
users can inspect why a route ranked highly.

## Signal-propagation constraint

The start and target are always permitted because they are endpoints. The
default policy is sensitivity-oriented: every internal node must have at least
one of these class labels in the active node universe:

- `receptor`
- `receptor_regulator`
- `ligand`
- `kinase`
- `phosphatase`
- `cyclase`
- `phosphodiesterase`
- `phospholipase`
- `nos`
- `gtpase`
- `gtpase_regulator`
- `kinase_phosphatase_binding`
- `second_messenger_binding`
- `signaling_process`
- `signaling_regulation`
- `second_messenger`

This deliberately includes broad ligand, binding, and signaling-process
annotations in order to favor sensitivity over specificity. The only node-
universe class that does not qualify is `adaptor_scaffold`. Consequently, a
protein whose sole class is `adaptor_scaffold` cannot be an intermediate. Any
multi-role scaffold carrying at least one additional signaling annotation
remains eligible. The optional
`--exclude-any-scaffold` flag implements the more literal alternative in which
every scaffold-tagged intermediate is removed regardless of its other roles.
In the active universe, 195 nodes carried `adaptor_scaffold`: 94 had no other
class and were excluded by the default rule, whereas 101 had at least one
additional signaling class and remained eligible.

The graphical workbench makes this policy configurable. Its path panel lists
all ontology-derived role classes, shows the GO root term or terms used to
build each class, and allows any subset to be selected. `adaptor_scaffold` is
therefore available for an explicit sensitivity analysis even though it is off
by default. A node with several class labels is eligible when at least one of
its labels is selected. An empty selection permits no internal nodes (direct
start-to-target paths can still be returned), whereas disabling the
signaling-intermediates-only rule ignores the class selection. The separate
exclude-any-scaffold control takes precedence over the class selector.

The normalized selection is recorded as `allowed_intermediate_classes` in both
`configuration.json` and the path-finding portion of `analysis_summary.json`.
The per-node eligibility audit records the selected class matches and the final
decision, so a scientist can reproduce exactly which ontology policy admitted
each internal node.

For reproducibility, `--allow-all-intermediates` disables the constraint. Each
run writes `intermediate_node_eligibility.tsv`, which records the classes,
matching relay classes, eligibility decision, endpoint exemption, and reason
for every matrix node.

## Prkaca-to-Aqp2 sensitivity-oriented reference run

The current reference run used the GUI workflow defaults with:

- start `Prkaca`;
- target `Aqp2`;
- 25 requested paths;
- at most six hops; and
- edge probability strictly greater than 0.5.

The validated Aqp2 incremental characterization was reused. Of the 892 matrix
nodes, 797 qualified as intermediates and 95 did not; Aqp2 was nevertheless
permitted as the selected endpoint. The baseline-filtered graph contained
169,769 unique undirected edges before the intermediate rule and 139,613 edges
in the permitted search subgraph afterward. The constrained start-node
component contained 792 nodes and Aqp2 had 323 permitted retained neighbors.
All 25 requested paths were found and passed the path validation checks.
Twenty-four of the 25 routes used at least one intermediate that the earlier
conservative policy would have excluded. Newly admitted examples included the
PKA regulatory subunits Prkar1a, Prkar2a, and Prkar2b; the transcriptional
regulators Creb1, Foxo3, Rela, and Crebbp; and the binding/regulatory proteins
Akap1, Bad, Ctnnb1, and Sufu. This confirms that the liberal policy materially
increases sensitivity and also makes the ranked paths more dependent on broad
ontology annotations.

The five highest-ranked paths were:

| Rank | Hops | Product score | Path |
|---:|---:|---:|---|
| 1 | 3 | 0.998812 | Prkaca -> Prkar1a -> Prkacb -> Aqp2 |
| 2 | 3 | 0.998716 | Prkaca -> Prkar2b -> Prkacb -> Aqp2 |
| 3 | 3 | 0.998706 | Prkaca -> Lipe -> Prkacb -> Aqp2 |
| 4 | 3 | 0.998663 | Prkaca -> Foxo3 -> Prkacb -> Aqp2 |
| 5 | 3 | 0.998650 | Prkaca -> Akap1 -> Prkacb -> Aqp2 |

The common-language input `PKA subunit a` is supported as an explicit alias for
mouse `Prkaca`; the analysis summary records when that alias is used. This is a
narrow curated convenience, not a general natural-language gene resolver.

The output directory
`results/gui_runs/Prkaca_to_Aqp2_sensitive_signaling_components/` contains one
record per ranked route, a separate edge-level audit table, the complete node
eligibility audit, a JSON summary with hashes and validation checks, and an
interpretation README. The earlier unconstrained result was preserved in the
separate `Prkaca_to_Aqp2/` directory rather than overwritten.

## Current workbench implementation

The standalone target-vector files and reference counts above are retained as
historical validation artifacts. The current interactive workbench routes an
outside endpoint through the persistent incremental edge-pair database instead
of maintaining a separate one-off target vector. The endpoint is treated as an
added protein, all target-to-selected-node pairs absent from the database are
characterized once, and subsequent runs rescore the cached raw evidence.

For Aqp2 and the 891-node seed, the incremental implementation created exactly
891 unordered rows, a 892 by 892 matrix containing 397,386 unique hypotheses,
and all 25 requested signal-relay-constrained paths. Dynamic KinasePredictor
edges use global top-ten model rank for cache stability, so current path scores
need not be numerically identical to the earlier fixed-seed target-vector run.
The incremental schema and reference policy are documented in
`docs/edge_characterization_methods_891.md`.

## Interpretation limitation

The product score is appropriate for ranking under the model, but it must not
be described as a calibrated probability that the complete biological route
exists. Edge evidence streams can be dependent, and edges within a route need
not be statistically independent. Routes should therefore be compared as
model-supported hypotheses and then reviewed for biological plausibility.
