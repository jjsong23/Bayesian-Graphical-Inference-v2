# Edge characterization of the revised 891-node signaling universe

## Scope and graph definition

The active node universe was regenerated after revising the protein-abundance
and principal-cell RNA preprocessing. Node selection incorporated four streams:
mpkCCD protein abundance, positive-only principal-cell median TPM, PKA-knockout
kinase-activity evidence, and PKA-knockout differential-phosphoprotein evidence.
The two PKA-derived streams contributed 32 kinase-activity-supported proteins
and 250 differential-phosphoprotein-supported proteins; three proteins were
supported by both streams. The retained universe contained 871 gene/protein
nodes plus 20 curated second-messenger molecules, for 891 total nodes.

All edge analyses used an undirected graph. The number of unique unordered,
non-self node pairs was therefore:

$$
\binom{891}{2}=396{,}495.
$$

Every adjacency matrix used the exact row/column order in
`data/node_selection/node_universe_combined_nonzero.tsv`. Matrices were required
to be symmetric, to have a zero diagonal, and to contain probabilities between
zero and one. Evidence missing for a pair was neutral: it preserved the
pre-stream posterior instead of being interpreted as evidence against an edge.

## Bayesian edge update

Each evidence stream supplied an edge-specific likelihood or Bayes factor. For
an existing edge probability \(p\), a likelihood \(L\), and no-edge reference
likelihood \(L_0=0.5\), the Bernoulli update was:

$$
p^\prime =
\frac{pL}{pL+(1-p)L_0}.
$$

Equivalently, for a stream Bayes factor \(BF=L/L_0\):

$$
\operatorname{logit}(p^\prime)
=
\operatorname{logit}(p)+\log(BF).
$$

The first stream began from \(p=0.5\). Subsequent streams updated the posterior
sequentially in the order described below. Only factors greater than one added
positive evidence. No negative factors were introduced.

## Stream 1: mpkCCD fractionation localization

Five-component subcellular-fraction abundance profiles were obtained from
`Proteomics_of_subcellular_fractions.xlsx`. Exact duplicate profiles were
collapsed. When a symbol had multiple distinct profiles, the profile with the
greatest summed abundance was selected.

Usable five-fraction profiles were available for 712 of 891 nodes. For each
pair of profiled nodes, the raw similarity was the dot product of the two
fraction vectors. For each profiled target node \(i\), \(T_{0.75,i}\) was the
75th percentile of its dot products with all other profiled nodes in the
current 891-node universe, excluding self-comparisons.

The directed likelihood calculated using node \(i\)'s threshold was:

$$
L_{j\rightarrow i}
=
\max\left[
0.5,\,
1-\exp\left(
-\frac{1}{2}
\left(\frac{s_{ij}}{T_{0.75,i}}\right)^2
\right)
\right].
$$

Reciprocal directed likelihoods were averaged to produce one undirected
likelihood. Missing profiles were neutral. This stream produced 98,328 unique
pairs with posterior probability greater than 0.5; 298,167 pairs remained
exactly at 0.5.

## Phosphosite database used for KinasePredictor

All 871 protein nodes mapped exactly to the selected accession in the cached
Mus musculus UniProt reference proteome UP000000589. The rebuilt database
contained 4,643 unique experimentally observed or curated phosphosites across
696 proteins. These sites combined PKA-knockout dataset observations and
UniProt phosphosite annotations. Sequence-derived S/T/Y residues that lacked
observed or curated phosphosite support were not used for edge predictions.

Of the 4,643 retained sites, 4,513 had a complete centered 13-residue window
and were scorable. The remaining 130 sites were preserved in the audit but not
scored.

## Stream 2: observed-site KinasePredictor

Each scorable phosphosite was evaluated with the official KinasePredictor v0.8
matrices: 237 serine/threonine matrices or 98 tyrosine matrices, depending on
the central residue. The site-specific \(T_{0.75}\) background was calculated
from the complete applicable matrix-score distribution before restricting
predictions to kinases in the node universe. Thirty sites had a nonpositive raw
75th percentile and used the recorded positive-score 75th-percentile fallback.

Predictor labels mapping to the same mouse kinase gene were collapsed by their
maximum score. For every site, the top ten distinct eligible kinase genes were
retained, where eligibility required a protein node labeled `kinase` in the
current universe. Seventy-five of 127 kinase-class nodes had at least one
mapped predictor model.

Each retained site assignment received:

$$
L_{\mathrm{site}}
=
\max\left[
0.5,\,
1-\exp\left(
-\frac{1}{2}
\left(
\frac{\max(\mathrm{score},0)}{T_{0.75,\mathrm{site}}}
\right)^2
\right)
\right].
$$

Evidence from multiple sites was accumulated on the unordered
\(\{\mathrm{kinase},\mathrm{protein}\}\) pair by summing
\(\log(L_{\mathrm{site}}/0.5)\). Self-predictions were audited but excluded from
the diagonal. Biochemical kinase/target roles were retained for provenance but
did not create direction in the graph.

The analysis recorded 19,301 unique kinase-protein pairs with at least one
top-ten prediction. Of these, 17,846 had a Bayes factor greater than one and
strengthened the localization posterior; 1,455 had only neutral ranked hits and
left the localization posterior unchanged. After this stream, 113,277 of all
396,495 pairs were above 0.5.

## Stream 3: STRING v12 functional associations

Protein nodes were mapped to the Mus musculus STRING v12 network. A total of
870 of 871 proteins mapped; `Cdk3` remained unmapped. The 20 molecule nodes
were ineligible for STRING evidence.

For a reported STRING combined score \(s\), the evidence Bayes factor used
STRING's documented prior \(\pi=0.041\):

$$
BF_{\mathrm{STRING}}
=
\frac{s/(1-s)}{\pi/(1-\pi)}.
$$

The combined score was used once; individual STRING channels were retained for
audit and were not added separately. Absence of a STRING row was neutral.

STRING reported 65,764 unique relationships between mapped proteins in the
current universe. All 65,764 increased their pre-STRING posterior. Of these,
12,324 met the descriptive medium-confidence threshold of 0.4 and 3,665 met
the high-confidence threshold of 0.7. After STRING integration, 147,367 pairs
were above 0.5 and 249,128 remained exactly at 0.5.

## Stream 4: Human Protein Atlas localization

Mouse proteins were mapped to HPA v25.1 through MGI's stringent one-to-one
protein-coding mouse-human orthology report. No symbol-only or paralog fallback
was used. Primary HPA profiles included Enhanced, Supported, and Approved
locations and excluded Uncertain annotations. Each binary 49-location profile
was L2-normalized, so pairwise dot products were cosine similarities rather
than counts of shared annotations.

One-to-one human orthologs were available for 823 protein nodes. Primary HPA
profiles were available for 681 proteins, producing 231,540 pairs with two
profiles. For each profiled node, \(T_{0.75}\) was the 75th percentile of its
cosine similarities to other profiled nodes in the current universe. When this
threshold was zero, the 75th percentile of positive overlaps was used and the
fallback was recorded. A node with no positive overlap remained neutral.

Endpoint-specific likelihoods used the same complement-of-minimum transform as
the mpkCCD localization stream. The two endpoint likelihoods were averaged,
and \(BF_{\mathrm{HPA}}=L_{\mathrm{HPA}}/0.5\) updated the existing posterior.

HPA primary localization strengthened 35,935 unique pairs. Of those, 20,950
had been exactly at 0.5 before HPA and 14,985 were already above 0.5. The final
graph contained:

- 168,317 unique edges above 0.5;
- 228,178 unique pairs exactly at 0.5;
- 396,495 unique pairs in total.

The high-confidence sensitivity analysis used only Enhanced and Supported HPA
annotations. It produced 497 profiled nodes, 123,256 covered pairs, and 14,470
pairs with a Bayes factor greater than one. It is an alternative to the primary
HPA stream and must not be combined simultaneously with it.

## Stream 5: OmniPath core signaling interactions

The OmniPath interaction web service was queried on 2026-07-30 for the curated
core `omnipath` dataset, mouse NCBI taxonomy 10090, post-translational
interactions, protein endpoints, and directed records. The query returned
31,384 directed interaction records. OmniPath is built primarily from human
knowledge and provides the mouse network by orthology translation; this
cross-species provenance was retained as a limitation.

Endpoints were mapped first by exact mouse gene symbol. An exact mouse UniProt
accession from the project protein index was used only when the API gene-symbol
field did not match a node; this recovered four additional unordered pairs.
The mapped data contained 1,242 directed rows connecting in-universe proteins.

Although the source data were directed and could be annotated as stimulating,
inhibitory, or unsigned, the modeled graph remained undirected and unsigned.
All rows mapping to the same unordered pair were collapsed to one pair. The
maximum API `curation_effort` across directional rows was used as a
conservative support score, while all directions, signs, source resources, and
literature references were unioned into the audit table. This produced 1,163
unique supported protein pairs incident to 394 protein nodes. Sixty-one pairs
had records in both directions. The audit contained 457 pairs with at least
one stimulation annotation, 333 with at least one inhibition annotation, 464
with at least one unsigned annotation, and 83 with both stimulation and
inhibition annotations. These categories were not used to orient or sign the
edge.

The background distribution matched the undirected scoring unit. All
non-self interactions in the complete mouse query were collapsed to 29,387
unique unordered gene-symbol pairs, and each pair was assigned its maximum
`curation_effort`. The 75th percentile was:

$$
T_{0.75}=6.
$$

Because an explicit curated OmniPath record is positive evidence even when its
curation effort is low, the sparse-presence likelihood used a zero-anchored
complement transform:

$$
L_{\mathrm{OmniPath}}
=
0.5
+
0.5\left[
1-\exp\left(
-\frac{1}{2}
\left(
\frac{\mathrm{curation\ effort}}{T_{0.75}}
\right)^2
\right)
\right].
$$

Thus, absence from the mapped OmniPath network was exactly neutral
(\(L=0.5\), \(BF=1\)), while every curated positive score produced
\(L>0.5\) and \(BF>1\). The factor and update were:

$$
BF_{\mathrm{OmniPath}}=\frac{L_{\mathrm{OmniPath}}}{0.5},
$$

$$
\operatorname{logit}(P_{\mathrm{final}})
=
\operatorname{logit}(P_{\mathrm{after\ HPA}})
+
\log(BF_{\mathrm{OmniPath}}).
$$

Supported Bayes factors ranged from 1.01379 to 2.0, with a median of 1.11750.
Of the 1,163 supported pairs, 1,083 were already above baseline after the
first four streams and 80 were newly raised above 0.5. The final graph
therefore contained:

- 168,397 unique edges above 0.5;
- 228,098 unique pairs exactly at 0.5;
- 396,495 unique pairs in total.

OmniPath integrates many literature and database resources. Some records can
overlap information already represented in STRING or the
kinase/phosphosite-derived streams. Treating these streams as conditionally
independent may therefore overstate combined evidence. OmniPath was retained
as a separately toggleable sparse factor stream.

## Stream 6: STITCH v5 secondary-messenger associations

STRING is restricted to protein–protein associations, so its companion
protein–chemical database STITCH was used to characterize edges incident to
the 20 curated secondary-messenger nodes. Each messenger was mapped to an
explicitly audited PubChem compound identifier. Generic lipid messengers were
represented by the corresponding PubChem class record and flagged as such in
the mapping audit.

The Mus musculus STITCH v5 bulk file contained 39,116,252 protein–chemical
records. Exact PubChem matching retained 26,430 rows. STITCH merged-compound
(`CIDm`) and stereospecific (`CIDs`) representations, alternate PubChem records
for one curated messenger, and repeated records were collapsed to one
messenger–protein pair by taking the maximum combined score. Mouse Ensembl
protein identifiers were mapped through the existing audited STRING v12 node
mapping. This produced 1,037 unique in-universe messenger–protein pairs across
14 messenger nodes.

STITCH combined scores range from its reporting floor of 0.150 to 0.999 in the
retained data. A score (s) was converted to positive-only evidence relative
to the reporting floor (r=0.150):

$$
BF_{\mathrm{STITCH}}
=
\max\left[
1,
\frac{s/(1-s)}{r/(1-r)}
\right].
$$

The floor prevents an absent or floor-valued association from being treated as
evidence against an edge. The graph update was:

$$
\operatorname{logit}(P_{\mathrm{final}})
=
\operatorname{logit}(P_{\mathrm{after\ OmniPath}})
+
\log(BF_{\mathrm{STITCH}}).
$$

Of the 1,037 mapped pairs, 1,017 had (BF>1); the other 20 were exactly at the
0.150 reference and remained neutral. Because messenger-incident pairs had not
previously received protein-based association evidence, all 1,017 non-neutral
pairs were newly raised above 0.5. The updated graph therefore contained:

- 169,414 unique edges above 0.5;
- 227,081 unique pairs exactly at 0.5;
- 396,495 unique pairs in total.

STITCH combined scores may include curated databases, experiments,
predictions, and text mining. They are functional protein–chemical association
confidence scores, not direct-binding probabilities or directional signaling
claims. STITCH v5 is also an older bulk release, so absence was always treated
as neutral rather than evidence of no relationship.

The OmniPath `small_molecule` dataset was evaluated independently. Exact
PubChem matching found 28 human records for seven messengers, but none of their
protein targets mapped into the current 891-node universe through the existing
one-to-one human–mouse crosswalk. It therefore contributed no second
small-molecule factor. This null coverage audit is preserved under
`data/edge_characterization/omnipath/2026-08-03/`.

## Optional derived stream: scaffold-mediated triadic closure

An optional second-pass stream represents the hypothesis that two proteins
with physically supported interactions to the same scaffold are more likely to
occupy the same signaling complex or local microenvironment. This is a
scaffold-specific proximity inference; it is not interpreted as proof that the
two inferred partner proteins bind each other directly.

The biological rationale is that multivalent scaffold proteins can nucleate
signaling complexes and increase partner proximity
([Good et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3137898/)). This
hypothesis is represented by a deliberately simple binary rule, not by a
degree-corrected network-link predictor.

The operation uses only protein nodes. A common-neighbor node qualifies as a
scaffold when its semicolon-delimited node classes contain the exact
`adaptor_scaffold` token. Multi-role scaffold proteins remain eligible. For
proteins (i) and (j), scaffold (s), and a physical anchor score (A), the pair
qualifies exactly when:

$$
\exists s:\quad A_{is}>c\ \land\ A_{js}>c.
$$

The only permitted anchors are accepted-protocol AlphaPulldown/AlphaFold scores
or reported-positive HuRI interactions. The default AlphaPulldown/AlphaFold
cutoff is (c=0.90), exposed as `Anchor >`; a score exactly at 0.90 does not
qualify. A reported-positive HuRI pair qualifies directly. Localization,
STRING, OmniPath, kinase prediction, STITCH, and the integrated cheap-edge
posterior cannot populate (A). Scaffold degree, the amount by which an anchor
exceeds 0.90, and the number of shared scaffolds do not alter the assigned evidence.
Every qualifying pair receives the support likelihood (L=0.90), exposed as
`Support L`. Relative to the neutral likelihood of 0.50, the stored factor is:

$$
BF_{ij}^{\mathrm{closure}}=\frac{0.90}{0.50}=1.8.
$$

Nonqualifying pairs receive neutral (BF=1) by omission from the sparse factor
table. The factor is appended to the pre-closure log odds with the selected
stream weight. Closure is calculated exactly once: a closure-derived edge is
never reused as an anchor, preventing recursive graph densification. There is
no empirical (T_q), degree penalty, continuous shared-scaffold score, or
noisy-OR aggregation.

The Version 2 frontier recalculates this rule from its persistent physical-pair
cache. The audit retains every qualifying pair's shared-scaffold list, fixed
likelihood/BF, and both anchor sources/protocols. The older validation under
`results/gui_runs/scaffold_binary_closure_validation_20260804_v2/` used the
general integrated edge posterior as (P) and is now historical/superseded; its
counts are not results of the physical-only method.

## Validation and backend storage

All six cumulative adjacency matrices were independently checked against the
active node list. Each matrix was 891 by 891, had matching row and column node
order, was symmetric, had a zero diagonal, and contained only values in
\([0,1]\). The complete final undirected edge table contained exactly 396,495
unique pairs. Pairs absent from STITCH were verified to be unchanged, and each
STITCH-supported messenger–protein pair was represented once as an unordered
edge.

Complete sparse factor tables are stored under:

`results/backend_bayes_factor_catalogs/edge_factors_891/`

Each primary factor file contains every unique edge with \(BF>1\), including
`omnipath_core_bf_gt1.tsv.gz` with all 1,163 collapsed OmniPath pairs and
`stitch_secondary_messenger_bf_gt1.tsv.gz` with all 1,017 non-neutral STITCH
pairs. Absence
from a complete sparse factor file means \(BF=1\) for that stream. The explicit
neutral KinasePredictor file separately preserves evaluated top-ten pairs whose
Bayes factor remained one.

The Excel catalog generated before the OmniPath and STITCH extensions,
`results/backend_bayes_factor_catalogs/edge_characterization_bayes_factor_catalog.xlsx`
contains the node order and the first four evidence streams. The complete
current OmniPath and STITCH factors are stored in the sparse backend tables
described above. The final full adjacency matrix is:

- `results/edge_characterization/localization_kinase_predictor_string_hpa_omnipath_stitch/combined_adjacency_matrix.tsv`

The predecessor 848-protein/868-node outputs were preserved before rerunning
the workflow under:

`results/archive/2026-07-30_edge_characterization_868_nodes/`

## Incremental extension beyond the 891-node seed

The validated 891-node graph and its 396,495 unordered hypotheses are retained
as an immutable seed catalog. The interactive workflow is no longer restricted
to those nodes. If a change in node-evidence selection makes an additional
protein non-neutral, that protein is inserted and every previously unseen
unordered pair incident to it is characterized before edge integration. An
external path endpoint outside the seed, such as Aqp2, follows the same route.

For (n) already selected nodes and (k) newly added nodes, the required new
pair count is

$$
kn + \binom{k}{2}.
$$

For example, adding eight proteins to the 891-node seed required exactly
(8(891)+\binom{8}{2}=7{,}156) new pair calculations. Adding a ninth protein
later required only 899 additional pairs; the prior 7,156 rows were reused.

### Persistent pair database

The running cache is a local SQLite database:

`data/edge_characterization/incremental_edge_cache/edge_pair_cache.sqlite3`

The database stores the 396,495 canonical seed pair keys and, separately, one
raw-evidence row for every characterized incremental pair. Pair keys are
canonical unordered gene-symbol pairs, enforced as unique by the database
primary key. Incremental rows retain:

- mpkCCD localization dot product and both node-specific thresholds;
- the retained KinasePredictor site hits, raw motif scores, and site-specific
  thresholds;
- STRING combined score;
- HPA primary and high-confidence cosine similarities and both reciprocal
  thresholds;
- OmniPath maximum curation effort; and
- STITCH combined score for messenger–protein pairs.

Storing raw statistics rather than only final posteriors permits each GUI
stream weight and (T_q)/reference multiplier to be changed without rescanning
the source datasets. Requested Bayes factors are recalculated from cached raw
values, after which prior odds are multiplied exactly as for seed pairs.
Missing evidence remains (BF=1).

Each cache row is associated with an evidence signature derived from the cache
schema version plus the paths, sizes, and nanosecond modification times of all
raw sources and kinase-model files. Only rows with the current signature are
reused. A source change therefore causes requested pairs to be characterized
under a new signature rather than silently mixing data versions. The source
manifest and scientific policies are written to `cache_manifest.json`.

### Stable incremental reference policy

New node-specific localization and HPA thresholds are calculated against the
covered members of the fixed 891-node background. They are not recalculated
when subsequent nodes are added. This makes a cached pair stable and prevents
the order of node additions from altering earlier evidence.

Dynamic KinasePredictor edges similarly use models whose global rank is within
the top ten for a phosphosite. For a new protein target, its UniProt-annotated
and PKA-KO-observed phosphosites are scored against eligible selected kinase
models. If an added node is itself a kinase, existing target phosphosites are
also scored for that new kinase. The global-top-ten restriction is independent
of the changing selected-kinase set and therefore supports safe caching. This
incremental rule differs from the historical seed calculation, which selected
the top ten among kinases eligible in the fixed seed universe.

STRING detailed mouse links are scanned once for each batch of missing pairs.
HPA, mpkCCD localization, OmniPath, and STITCH sources are likewise loaded only
when at least one requested pair is absent. Once every requested pair is in the
database, a repeat run performs no pair characterization.

### Incremental validation

The first validation expansion selected eight proteins beyond the seed and
produced a 899 by 899 symmetric matrix with exactly 403,651 unique pairs. All
7,156 incremental pairs were computed on the first run; an identical second
run reported 7,156 cache hits and zero new calculations. Expanding to 900 nodes
created exactly 899 additional rows. The resulting matrix contained exactly
404,550 unique hypotheses, was symmetric, had a zero diagonal, and remained
within \([0,1]\).

The external Aqp2 validation created 891 cached target-to-seed pairs, a 892 by
892 matrix with 397,386 unique hypotheses, and 25 valid ranked paths. These
checks confirm that node-selection additions and external targets use the same
incremental implementation.
