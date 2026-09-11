# Graphical Bayesian Inference Project: Detailed Methods

This document was reconstructed from the implemented scripts, intermediate
tables, and stored validation summaries. It reflects the calculations actually
performed, including the distinctions between the node-selection model, the
general edge-characterization model, and the separate colocalization screen
used to prioritize AlphaFold predictions.

## Project objective and analytical structure

The project was designed to construct a graphical model of signaling in renal
principal cells, with two primary phases:

1. **Node selection:** identify proteins likely to participate in signaling.
2. **Edge characterization:** estimate the probability that two selected nodes
   are functionally related or capable of interacting.

A third, related analysis was developed specifically to reduce the number of
protein pairs submitted to structural prediction. This
**colocalization-screening graph** used subcellular localization information to
identify protein pairs with sufficient spatial compatibility for potential
interaction.

These analyses used related Bayesian-style calculations. As of the 2026-08-10
model revision, node and edge selection use the same hypothesis structure:

- Each protein is an independent Bernoulli hypothesis comparing presence in
  the signaling system with absence from it.
- Each unordered edge is an independent Bernoulli hypothesis comparing the
  presence and absence of a relationship.
- Both begin at prior probability 0.5 by default and are updated in odds space.
  Neither probability vector is normalized across proteins or pairs.
- Earlier stored node-posterior artifacts used a sum-to-one relative-weight
  normalization. Those files are historical and should not be interpreted as
  outputs of the revised binary node model.
- The AlphaFold filter used localization, experimental prior knowledge, and
  scaffold classification as screening criteria rather than treating all three
  as calibrated interaction-probability measurements.

# Part I: Node selection

## 1. Construction of the initial signaling-candidate universe

A liberal, high-recall universe of mouse signaling proteins was constructed
from a mouse UniProt Gene Ontology Annotation file and a local Gene Ontology
OBO file.

Sixteen overlapping functional classes were defined using 28 GO root terms:

- receptor
- ligand
- receptor regulator
- kinase
- phosphatase
- GTPase
- GTPase regulator
- cyclase
- phosphodiesterase
- phospholipase
- nitric-oxide synthase
- adaptor/scaffold
- kinase/phosphatase binding
- second-messenger binding
- signaling process
- signaling regulation

Each root was expanded to include all descendants reachable through asserted
`is_a` and `part_of` relationships. No ontology reasoner was used, so inferred
relationships not present in the supplied ontology were not added.

A protein was included if any of its GO annotations matched one of the expanded
term sets. Molecular Function and Biological Process annotations were both
accepted. Rows were excluded only when:

- the annotation was explicitly negated with `NOT`, or
- the evidence code was `ND`, indicating no biological data.

All other GO evidence codes were retained. Experimental evidence was recorded
as a separate flag but was not required for inclusion.

The initial export contained:

- 9,175 UniProt gene-product records
- 9,170 unique mouse gene symbols
- 2,007 GO terms in the expanded inclusion network

Five symbols—`Akap7`, `Calca`, `Cdkn2a`, `Gnas`, and `Nrxn1`—were represented
by more than one UniProt record. The scoring universe was deduplicated by gene
symbol to produce 9,170 candidates. Downstream metadata joins retained the
first symbol-matched record rather than merging all duplicated metadata; this
should be disclosed if isoform-level distinctions become important.

The initial node prior is the objective two-hypothesis prior:

$$
P_i^{(0)}=0.5
$$

for every candidate protein \(i\). This does not constrain the expected number
of present nodes to one; every protein has its own present-versus-absent
hypothesis.

## 2. General node-evidence scoring function

Protein abundance, principal-cell transcript abundance, kinase activity, and
phosphoprotein evidence were converted into positive-or-neutral evidence scores
using the same basic kernel.

For a measurement \(x\) and empirical threshold \(T_q\):

$$
f(x;T_q)=\max\left[0.5,\;1-\exp\left(-\frac{1}{2}
\left(\frac{x}{T_q}\right)^2\right)\right]
$$

The minimum score was fixed at 0.5. Consequently:

- \(f=0.5\) represented neutral evidence.
- \(f>0.5\) represented positive evidence.
- No stream generated a score below 0.5, so none supplied explicit negative
  evidence.

The source score is converted to a conventional Bayes factor by comparing it
with the neutral likelihood 0.5:

$$
BF_i=\frac{f_i}{0.5}.
$$

Thus, a source score of 0.5 becomes neutral \(BF=1\), while a score approaching
1 becomes \(BF\approx2\). Each node is updated independently:

$$
O_i^{\mathrm{new}}
=O_i^{\mathrm{old}}BF_i^w,
\qquad
P_i^{\mathrm{new}}=\frac{O_i^{\mathrm{new}}}{1+O_i^{\mathrm{new}}}.
$$

In the reproducibility default, missing measurements receive \(BF=1\) rather
than being interpreted as evidence that the protein is absent. Posterior
probabilities are not divided by a sum over other proteins.

The stored values between 0.5 and 1 are most precisely described as
**Bayes-factor-like likelihood scores**. The revised integrator explicitly
divides them by their neutral value before updating odds, so the operational
Bayes factors are neutral at 1.

### Optional negative evidence from nondetection

The workbench additionally provides a deliberately simple, optional rule for
using nondetection as evidence against node presence. When `penalize_unobserved`
is enabled, an eligible candidate that is not observed in active stream \(j\)
receives a user-configured nondetection Bayes factor \(B^-\), constrained to

$$
0 < B^- \leq 1.
$$

The default value is \(B^-=0.5\), representing two-to-one evidence against the
presence hypothesis. The effective stream Bayes factor is

$$
BF_{ij}^{\mathrm{effective}}=
\begin{cases}
B^- & \text{if candidate }i\text{ is eligible but unobserved in stream }j,\\
BF_{ij}^{\mathrm{source}} & \text{if it is observed},\\
1 & \text{if stream }j\text{ cannot assess that candidate.}
\end{cases}
$$

It then enters the same weighted odds update used by every other node and edge
factor. From a 0.5 prior, one unweighted nondetection at \(B^-=0.5\) gives a
posterior of \(1/3\); two give 0.20; three give approximately 0.111; and four
give approximately 0.0588. Repeated nondetections therefore drive the posterior
toward zero. An exact zero BF is intentionally prohibited because it would make
the absence conclusion irreversible regardless of later positive evidence.

Eligibility prevents absence from being inferred where the assay is
conceptually inapplicable. All protein candidates are eligible for the protein,
RNA, phosphoprotein-response, and segment-abundance streams. Only proteins with
the `kinase` ontology class are eligible for kinase-activity nondetection;
non-kinases remain neutral. Curated small-molecule second messengers are added
after protein scoring and are not penalized by protein-assay nondetection.

The option is disabled by default so historical node universes remain exactly
reproducible. Each run records the option, \(B^-\), per-stream eligible and
observed counts, every applied penalty, the original source BF, the effective
BF, and the weighted log BF in its configuration, summary, and node audit
table. Nondetection is not proof of biological absence; users should enable the
rule only when they accept the assay-coverage assumption.

The GUI also provides a per-stream continuous alternative. For an enabled
stream, it removes the 0.5 likelihood floor, evaluates eligible nondetections as
`x=0`, and allows low measured values to produce BF below 1. It supersedes the
fixed nondetection factor only for that stream. The positive measured Tq is
retained because zero-padding sparse candidate universes can make q75 equal to
zero. The exact equation, stream coverage, floor, and audit fields are specified
in `docs/continuous_negative_evidence.md`.

## 3. mpkCCD protein-abundance evidence

Protein abundance was obtained from the mpkCCD protein-abundance workbook:

- Sheet: `Original Data in webpage`
- Identifier: `Gene Symbol`
- Measurement: `Log10 Abundance`

The supplied values were first returned to the linear relative-abundance scale:

$$
x_i^{\mathrm{linear}}=10^{x_i^{\mathrm{Log10\ Abundance}}}
$$

This operation reverses the logarithmic representation in the workbook; it does
not recover raw peptide-ion intensities. The empirical background for \(T_q\)
was calculated before restricting the output to signaling candidates and
consisted of all 6,750 finite back-transformed abundance values.

The threshold was:

$$
T_{0.75}=1.2630340049\times 10^9
$$

Because the implemented score reaches the neutral floor until the transformed
survival probability exceeds 0.5, the effective boundary for a factor strictly
greater than 0.5 was \(1.4871088962\times10^9\), equivalent to a supplied
log10 abundance of 9.1723427717.

Of the 9,170 signaling candidates:

- 2,190 had an abundance measurement.
- 467 received a factor strictly greater than 0.5.
- 1,723 measured signaling candidates remained at the 0.5 floor.
- Unmeasured candidates received 0.5 during posterior integration.

Each candidate's posterior is updated independently from 0.5. Neutral or
missing abundance evidence leaves that candidate at its pre-stream
probability; abundance evidence for another protein cannot lower it.

A methodological detail worth reporting is that no additional normalization
was introduced. The only preprocessing step was the inverse log10 transform
described above.

## 4. Principal-cell transcriptomic evidence

The transcriptomic evidence was restricted to principal cells.

Data were taken from:

- Sheet: `Median TPM`
- Identifier: `Gene Symbol`
- Measurement: `PC (Median TPM,n=74)`

The source sheet contained 8,022 finite principal-cell median TPM values. Zero
values were treated as nondetections rather than quantitative observations:
3,900 zeros were excluded, leaving a background of 4,122 strictly positive
values. The background was still formed from all eligible genes in the sheet,
not only genes in the signaling universe. Zero or missing values were assigned
the neutral factor during integration and therefore did not reduce a node's
probability.

The threshold was:

$$
T_{0.75}=7.1225\ \text{TPM}
$$

Of the signaling candidates:

- 1,418 signaling candidates had a strictly positive principal-cell transcript
  measurement and were scored.
- 1,108 signaling-candidate zeros were excluded from scoring.
- 319 scored candidates received a factor above 0.5.
- 1,099 scored candidates remained at the 0.5 floor.
- Zero and missing values received the neutral factor.

Protein-abundance evidence was integrated first, followed by principal-cell
transcript evidence. After these two streams:

- 655 proteins had at least one protein-or-PC factor above the neutral floor.
- 8,515 proteins had both protein and PC factors at the neutral floor.

Under the revised model, a protein with only neutral evidence remains at
posterior 0.5, while positive evidence raises it above 0.5. The project's
historical term “nonzero nodes” should therefore be interpreted as **nodes with
posterior above the neutral baseline**, not nodes with a posterior greater than
zero.

## 5. KinasePredictor annotation of the PKA-knockout phosphoproteome

### 5.1 Phosphosite sequence scoring

KinasePredictor v0.8 was reproduced from the official ESBL matrix package.

For each valid single-site record:

1. The centralized sequence was reduced to a 13-residue window.
2. The central residue was required to be serine, threonine, or tyrosine.
3. Serine/threonine sites were scored against 237 serine/threonine matrices.
4. Tyrosine sites were scored against 98 tyrosine-kinase matrices.
5. For each kinase model, the 13 matrix values corresponding to the 13 amino
   acids were summed.
6. Kinases were ranked from highest to lowest motif score.
7. The top ten kinase models were retained.

The raw KinasePredictor score is a sequence-motif compatibility score. It is
not itself a probability, posterior, or Bayes factor.

Rows labeled as multiple sites, or rows lacking a valid centered 13-residue
sequence, were retained for auditing but were not scored.

The PKA-knockout dataset contained:

- 9,117 total rows
- 7,174 scored single-site rows
- 1,943 unscored rows
- 7,174 unique scored phosphosites
- 71,740 top-ten kinase assignments

No additional protein-abundance normalization or phospho-specific correction
was applied to the raw phosphosite log2 changes, following the selected
analysis plan.

### 5.2 Top-ten hit aggregation

Every kinase appearing in ranks 1–10 was treated as an equal-weight hit for
that phosphosite. Prediction rank affected whether a kinase entered the top ten
but did not supply an additional rank weight during kinase-level aggregation.

Phosphosites were identified by the unordered key:

$$
\text{UniProt accession} \mid \text{site}
$$

Duplicate records would have been collapsed by site key, using median raw log2
change and median site \(P\)-value. No duplicate site keys were present in this
dataset.

A total of 257 distinct KinasePredictor labels received at least one hit.
Kinases were required to have at least ten unique hit sites, leaving 185 kinase
labels for downstream analysis.

For each kinase, the following were calculated:

- number of unique hit sites
- mean and median raw phosphosite log2 change
- first and third quartiles
- positive, negative, and zero site counts
- proportion of sites with a nominal site-level \(P<0.05\)
- two-sided Wilcoxon signed-rank test against zero

The Wilcoxon calculation used a tie-corrected normal approximation. The 185
kinase-level \(P\)-values were adjusted by the Benjamini–Hochberg procedure.

Kinases were classified as:

- **Increased:** FDR < 0.05 and median LFC > 0
- **Decreased:** FDR < 0.05 and median LFC < 0
- **No consistent change:** otherwise

This produced:

- 49 increased kinase profiles
- 14 decreased kinase profiles
- 122 with no consistent signed change

These signed classifications were retained for interpretation but were not the
statistic used for node selection.

### 5.3 Absolute-LFC kinase evidence

For node selection, the kinase statistic was:

$$
x_k=
\left|
\operatorname{median}
\left(
\mathrm{raw\ phosphosite\ log_2(PKA-null/PKA-intact)}
\right)
\right|
$$

The background consisted of all 185 finite kinase-level absolute median LFC
values.

The threshold was:

$$
T_{0.75}=0.17
$$

The usual positive-or-neutral kernel was applied, producing:

- 32 kinase labels above the 0.5 floor
- 153 kinase labels at the floor

The raw score was converted to a neutral-relative multiplier:

$$
m_k=\frac{f_k}{0.5}
$$

Thus, a factor of 0.5 became multiplier 1, while a factor approaching 1 became
multiplier 2.

KinasePredictor labels were mapped to mouse gene symbols using the official
Kinase Logos metadata. Case-normalized symbols and a small explicit alias table
were used where necessary. Of the 185 predictor labels:

- 182 mapped to proteins in the signaling universe.
- These corresponded to 176 unique mouse gene nodes.
- Five genes received more than one predictor label.
- When multiple labels mapped to one gene, the largest evidence factor was
  selected; ties were resolved using absolute LFC, hit count, and label.

Every non-kinase receives multiplier 1. Under the revised independent binary
model, this leaves its posterior numerically unchanged. Evidence supporting a
kinase cannot lower a non-kinase merely by reallocating probability mass.

## 6. Differential-phosphoprotein evidence

A separate node-evidence branch asked whether a protein carried phosphosites
that changed strongly after PKA deletion, regardless of which kinase was
predicted to act on those sites.

The background consisted of the absolute raw LFC values for all 7,174 usable
single phosphosites in the PKA-knockout experiment. The background was not
restricted to signaling-universe genes.

Every observed site was scored independently against the same unadjusted
site-level threshold:

$$
T_{0.75}=Q_{0.75}
\left(\{|\mathrm{LFC}_{\mathrm{site}}|\}_{7174}\right)=0.73.
$$

For each signaling protein, the site with the maximum site-level Bayes factor
was retained. Because the transformation is monotonic in absolute LFC, this is
also the site with the maximum absolute raw LFC. A protein was a positive hit
when at least one site had BF greater than 1. Additional sites did not raise
the threshold and their factors were not multiplied. Site count, signed LFC,
UniProt accession, coordinate, and site-level \(P\)-value were retained only
for auditing.

Results were:

- 1,157 signaling proteins with at least one detected phosphosite
- 8,013 signaling proteins without detected sites
- 465 observed proteins above the neutral floor
- 692 observed proteins at the neutral floor

Proteins without detectable phosphosites received multiplier 1 and therefore
no direct phosphoprotein evidence.

Site-level \(P\)-values were retained for descriptive auditing but were not used
to calculate the node evidence factor. Fifteen site-level \(P\)-values were
invalid or unavailable. This site-centric rule replaced the earlier
`q^(1/n)` maximum-of-sites adjustment on 2026-08-20.

## 6A. Selective PKA-Cα- and PKA-Cβ-knockout extension

To recover signals that could be masked when both PKA catalytic subunits are
deleted together, a later mpkCCD phosphopeptide workbook was added as two
separately selectable node streams: PKA-Cα null versus matched intact and
PKA-Cβ null versus matched intact. The generic default leaves both streams off.

The source contained 4,635 phosphopeptide rows. Rows sharing a UniProt
accession and normalized `Site(s)` pattern were collapsed by the median signed
LFC independently for each comparison, producing 3,805 unique site-pattern
records. The median modified P value was retained for audit but was not used to
filter observations or calculate the evidence factor. Candidate mapping used
the signaling catalog's mouse UniProt accession first and exact supplied symbol
fallback second. This mapped 1,680 unique site patterns to 776 signaling
candidates.

For each comparison, the quantitative background contained the absolute LFCs
of all 3,805 duplicate-collapsed source site patterns before signaling-universe
filtering. The ordinary background q75 values were 0.22603867 for PKA-Cα KO and
0.18371188 for PKA-Cβ KO. Every site pattern was scored against its comparison's
single q75. Each mapped protein inherited the largest site-level factor, so one
positively scoring site pattern was sufficient and additional sites neither
raised its threshold nor accumulated extra evidence. Missing
proteins received neutral BF 1 unless the optional GUI nondetection penalty was
enabled.

At unit weight and `Tq × = 1`, the PKA-Cα comparison gave BF > 1 to 267
signaling candidates and the PKA-Cβ comparison gave BF > 1 to 261; 131 were
shared. With both new streams enabled alongside the four historical defaults,
1,232 proteins were selected and the 20 curated messengers produced 1,252
total nodes.

Absolute LFC was used because node selection asks whether a protein
participates, not whether it is activated or inhibited. Oppositely signed
alpha- and beta-specific responses therefore do not cancel. Their signed
values remain in the audit. The two streams share the dependence group
`pka_subunit_ko_phosphoproteomics` because they come from one TMT experiment;
their multiplied BFs must not be described as fully independent evidence.
These streams describe protein-level differential phosphorylation and do not
silently add new KinasePredictor-derived kinase-activity scores. Full details
and reproducibility paths are in `docs/pka_subunit_ko_node_evidence.md`.

## 7. Combination of kinase and phosphoprotein evidence

Both PKA-derived evidence streams were combined with the
protein-plus-principal-cell posterior:

$$
m_i^{\mathrm{combined}}
=
m_i^{\mathrm{kinase}}
m_i^{\mathrm{phosphoprotein}}
$$

$$
O_i^{\mathrm{final}}
=
O_i^{\mathrm{protein+PC}}
m_i^{\mathrm{combined}},
\qquad
P_i^{\mathrm{final}}
=\frac{O_i^{\mathrm{final}}}{1+O_i^{\mathrm{final}}}.
$$

No across-node normalization is applied.

A node was retained when at least one of the following was above its neutral
floor:

- protein-abundance factor > 0.5
- principal-cell transcript factor > 0.5
- kinase multiplier > 1
- phosphoprotein multiplier > 1

Under the current site-centric phosphoprotein rule, this produced 1,051
selected proteins. The validated pre-revision result contained 871 proteins
and remains the immutable seed for the precomputed edge catalog.

The stream-specific counts were:

| Evidence stream | Proteins above its neutral floor |
|---|---:|
| mpkCCD protein abundance | 467 |
| Principal-cell transcript abundance | 319 |
| Kinase activity | 32 |
| Differential phosphoprotein evidence | 465 |
| Kinase and phosphoprotein overlap | 3 |
| Any of the four streams | **1,051** |

The kinase and phosphoprotein streams were derived from the same PKA-knockout
phosphoproteomic experiment. Their product therefore assumes more independence
than can be guaranteed and may overstate evidence. The combined posterior was
explicitly marked as provisional for this reason.

## 8. Addition of second messengers

Twenty second-messenger molecules were added by a curated rule after protein
node selection:

- 7 canonical messengers
- 11 established noncanonical messengers
- 2 context-dependent messengers

Examples included calcium, cAMP, cGMP, IP3, DAG, PIP2, PIP3, nitric oxide,
reactive oxygen species, phosphatidic acid, sphingolipid messengers, cADPR,
NAADP, and 2′3′-cGAMP.

These molecules:

- were not gene products,
- did not receive protein/transcript/phosphoproteomic scores,
- were not assigned a node-selection posterior,
- were retained to make the signaling graph biologically more complete.

The current generic node-selection result therefore contains:

- 1,051 proteins
- 20 second-messenger molecules
- **1,071 total nodes**

The 1,051 proteins have overlapping functional labels. For example, 563 are
labeled `signaling_process`, 279 `kinase_phosphatase_binding`, 247
`adaptor_scaffold`, and 155 `kinase`. Class-composition plots counted a
multi-class protein in every class assigned to it; the class counts are
therefore not mutually exclusive.

## 8A. Collecting-duct extension for the Aqp2-focused analysis

The validated 891-node edge catalog remains the immutable seed. Under the
current site-centric phosphoprotein rule, the generic node defaults select
1,051 proteins and dynamically extend the graph to 1,071 total nodes after the
20 curated second messengers are added. A separate
lab-specific node-selection profile was created to represent Aqp2 signaling
throughout the collecting duct rather than only the cortical collecting duct.
Six additional evidence streams were defined by crossing two assays with three
anatomical segments:

- rat proteome: CCD, OMCD, and IMCD;
- mouse renal-tubule RNA-seq: CCD, OMCD, and IMCD.

Each assay/segment combination was retained as a distinct, independently
selectable GUI stream. The three segment streams within an assay were also
tagged with a shared dependence group because they originate from one source
experiment. Treating them as separate Bayes-factor updates follows the
requested sensitivity-oriented analysis, but their common experimental source
means that the resulting posterior magnitudes should not be interpreted as if
all six measurements were fully independent experiments.

### Rat-to-mouse identifier harmonization

The proteome contained 7,429 unique rat UniProt accessions and 7,429 unique rat
gene symbols. Rat proteins were harmonized to mouse genes with Ensembl BioMart
release 116 orthology data. Mapping was performed in the following order:

1. Match the rat UniProt accession to an Ensembl rat gene and then follow its
   mouse orthology relationship.
2. If no accession mapping was available, match the supplied official rat gene
   symbol to the Ensembl rat gene symbol.
3. If one or more high-confidence mouse orthologs were available for that rat
   protein, retain only those high-confidence mappings.
4. If no high-confidence mapping existed, retain the available low-confidence
   mapping(s) and label them explicitly as a fallback.
5. Retain all selected one-to-one, one-to-many, or many-to-many orthologs rather
   than inventing a single representative for a valid multi-gene relationship.

This procedure produced at least one named mouse ortholog for 7,107 of 7,429
rat proteins. Of these, 5,892 had a selected high-confidence ortholog and 2,678
mapped into the 9,170-protein signaling-candidate universe. The remaining 322
rat proteins were retained in an unmapped-protein audit. When more than one rat
protein or retained ortholog contributed to the same mouse candidate, the
maximum segment abundance was used. Abundances were not summed across proteins
or orthologs.

### Segment-specific backgrounds and Bayes factors

Zeros were common in both input tables and were interpreted as nondetection,
not evidence that a signaling node is absent. For each of the six streams,
\(T_{75}\) was calculated from all strictly positive finite measurements in
that source segment before restriction to the signaling-candidate universe.
Zero, missing, or nonmapped candidate values were excluded from the background
and received neutral evidence. For positive abundance \(x\), the source
likelihood and neutral-relative Bayes factor were:

$$
L(x)=\max\left(0.5,1-\exp\left[-\frac{1}{2}
\left(\frac{x}{T_{75}}\right)^2\right]\right),
\qquad
BF(x)=\frac{L(x)}{0.5}.
$$

Thus, nondetection gives \(L=0.5\) and \(BF=1\), while increasingly large
positive abundances approach \(BF=2\). No cross-node normalization was applied.
The empirical backgrounds and candidate counts were:

| Stream | Positive background values | Zero source rows | \(T_{75}\) | Observed signaling candidates | Non-neutral candidates |
|---|---:|---:|---:|---:|---:|
| Rat proteome CCD | 6,550 | 879 | 66,565.7131 | 2,391 | 521 |
| Rat proteome OMCD | 4,395 | 3,034 | 70,816.4444 | 1,643 | 325 |
| Rat proteome IMCD | 5,621 | 1,808 | 65,193.9136 | 2,128 | 476 |
| Mouse RNA CCD | 25,161 | 29,371 | 22.9 | 5,449 | 1,933 |
| Mouse RNA OMCD | 38,933 | 15,599 | 9.2 | 6,428 | 3,005 |
| Mouse RNA IMCD | 31,140 | 23,392 | 11.7 | 5,774 | 2,351 |

The mouse RNA file contained 100 duplicated nonempty gene symbols represented
by 203 rows. As with the proteome mappings, the maximum reported abundance for
each mouse symbol and segment was used; duplicated rows were not summed.

### Integration result

The six new factors were added to the same independent present-versus-absent
node model as the four existing streams. All nodes began with prior probability
0.5, source likelihoods were converted to neutral-relative Bayes factors, and
weighted Bayes factors multiplied node-specific prior odds. With all ten node
streams enabled at unit weight and their default \(T_q\) multipliers, cumulative
selection was:

| Stage | Selected proteins | New proteins at stage |
|---|---:|---:|
| Existing four streams | 1,051 | -- |
| + proteome CCD | 1,231 | 180 |
| + proteome OMCD | 1,238 | 7 |
| + proteome IMCD | 1,280 | 42 |
| + RNA CCD | 2,412 | 1,132 |
| + RNA OMCD | 3,261 | 849 |
| + RNA IMCD | **3,330** | 69 |

Adding the 20 curated second messengers produced a lab-specific collecting-duct
universe of 3,350 total nodes. This result is stored separately from the
generic profile; the additional six streams are selectable but off
in the generic GUI default configuration. The preprocessing and execution
entry points are `code/node_selection/build_collecting_duct_evidence.py` and
`code/node_selection/run_collecting_duct_node_selection.py`, respectively.

## 9. Historical alternate protein-plus-PC universe

Before the 2026-07-30 preprocessing revision, an alternate universe was
examined using only protein-abundance and principal-cell transcript evidence.
It contained 621 proteins whose historical protein-plus-PC posterior was
strictly above the exact minimum.

This alternate analysis excluded both:

- kinase-activity evidence, and
- differential-phosphoprotein evidence.

It was not simply a “no kinase activity” analysis.

The 621-protein colocalization graph was generated as a strict induced subgraph
of the previously calculated 848-protein graph. Edge values and the original
localization \(T_q\) backgrounds were preserved rather than recomputed after
subsetting.

The project subsequently returned to the kinase-inclusive universe. The
621-protein version should therefore be described as an archived sensitivity or
computational-reduction analysis, not the current primary universe.

> **Version boundary.** The node-selection revision above produced the current
> 871-protein/891-total-node universe. On 2026-07-30, the complete core
> edge-characterization workflow (mpkCCD localization, observed-site
> KinasePredictor, STRING v12, HPA v25.1, and OmniPath core) was recomputed or
> extended from source data for all 396,495 unique pairs. OmniPath source
> directions and stimulation/inhibition annotations were retained for audit
> but collapsed to unsigned, unordered pairs for the modeled graph. Current
> methods and counts are recorded in
> `docs/edge_characterization_methods_891.md`. The detailed narrative below
> retains the historical 848-protein/868-node counts for provenance and should
> not be used as the numerical description of the current graph.

> **Secondary-messenger edge extension (2026-08-03).** STRING itself is
> protein-only. The current workflow therefore adds STITCH v5 mouse
> protein–chemical associations as a sixth edge stream for the 20 curated
> secondary messengers. Exact, audited PubChem matching and duplicate collapse
> produced 1,017 non-neutral messenger–protein factors across 14 messengers.
> The resulting symmetric 891-node graph contains 169,414 unique edges above
> the 0.5 baseline. Full scoring details and coverage limitations are recorded
> in `docs/edge_characterization_methods_891.md`.

> **Incremental graph extension (2026-08-03).** The 891-node result is now an
> immutable seed rather than a hard maximum. Any additional non-neutral protein
> is inserted after all new unordered pairs are characterized. Raw pair
> evidence is retained in a versioned SQLite cache, so repeat runs rescore the
> stored evidence and compute only pairs never previously requested under the
> current source signature. New localization/HPA thresholds use the fixed
> 891-node background, and dynamic kinase edges require global top-ten motif
> rank. The full policy, schema, and validation counts are documented in
> `docs/edge_characterization_methods_891.md`.

# Part II: Phosphosite database for edge characterization

## 10. Protein-centered phosphosite database

A protein-centered phosphosite database was created for all 848 protein nodes
using the UniProt mouse reference proteome UP000000589.

All 848 proteins mapped to the selected UniProt accession exactly. The 20
small-molecule nodes were excluded because they have no amino-acid sequence.

For each protein, the database contained:

- selected and mapped UniProt accession
- canonical sequence and length
- UniProt annotated phosphorylated residues
- PKA-knockout experimentally observed phosphosites
- a union table of observed-or-annotated sites
- every sequence-derived S/T/Y residue
- centered 13-residue windows when available
- KinasePredictor scorable status
- a per-protein human-readable text file

The database contained:

- 4,729 unique observed-or-annotated phosphosites
- 3,537 UniProt-only annotated sites
- 317 PKA-knockout-only sites
- 875 present in both sources
- 679 proteins with at least one observed-or-annotated site
- 169 proteins with none
- 94,473 sequence-derived S/T/Y candidates
- 92,821 candidates with complete 13-residue windows

For subsequent kinase–protein edge inference, only the 4,729
observed-or-annotated sites were considered. The 94,473 sequence-derived
candidates were cataloged but excluded from edge prediction.

The phrase “observed phosphosites” should therefore be used carefully: the
edge-analysis table was the union of experiment-observed sites and
UniProt-annotated sites, not solely sites observed in the PKA-knockout
experiment.

# Part III: General edge-characterization graph

## 11. Graph structure and edge priors

The graph contained all 868 nodes. The number of possible unique unordered
pairs was:

$$
\binom{868}{2}=376{,}278
$$

Edges were modeled as undirected. Each pair appeared once in edge lists but
twice, symmetrically, in square adjacency matrices. Matrix diagonals were set to
zero by convention.

Every off-diagonal pair began with:

$$
P(\mathrm{edge})=0.5
$$

Evidence streams supplied Bayes factors or neutral values. In the validated
reproducibility default, an unsupported pair remained at its existing posterior
rather than being penalized. The interactive framework later added a separate,
disabled-by-default negative-edge sensitivity option described below.

For an evidence likelihood \(L\), the independent binary update was:

$$
P(E_{ij}\mid D)
=
\frac{
P(E_{ij})L
}{
P(E_{ij})L+
[1-P(E_{ij})](0.5)
}
$$

Equivalently:

$$
\operatorname{logit}(P_{\mathrm{new}})
=
\operatorname{logit}(P_{\mathrm{old}})
+
\log(BF)
$$

where:

$$
BF=\frac{L}{0.5}
$$

With the negative-edge option disabled, the stored default evidence streams are
positive or neutral and off-diagonal probabilities do not fall below 0.5.

### 11.1 Optional negative evidence for unsupported eligible pairs

The interactive workflow can let evidence counter the edge-present hypothesis.
When `penalize_unsupported` is enabled, every active primary edge stream first
defines the unordered pairs it could meaningfully assess. If an eligible pair
has no non-neutral factor record, that stream supplies a common configurable
factor (B^-_E), constrained to (0 < B^-_E \le 1). The default is 0.5. The
stream weight (w_k) is applied in the same odds equation used for positive
edge evidence:

$$
O(E_{ij}\mid D_1,\ldots,D_K)
=
O(E_{ij})
\prod_{k=1}^{K} BF_{ijk}^{w_k}.
$$

From a prior probability of 0.5, one unweighted unsupported-pair factor of 0.5
gives posterior (1/3); two give 0.2; three give approximately 0.111. Positive
factors can counter these penalties because all factors enter the same product.
The factor is strictly positive so evidence can make a posterior arbitrarily
small without creating an irreversible mathematical zero.

Eligibility was defined separately for each evidence source:

- mpkCCD localization: both proteins had valid fractionation profiles and
  positive node-specific thresholds;
- HPA primary or high-confidence localization: both proteins had valid profiles
  in the selected HPA policy and valid thresholds;
- KinasePredictor: one endpoint carried the exact `kinase` class and the other
  had at least one observed-or-UniProt-annotated scorable phosphosite;
- STRING: both endpoints were proteins with valid mouse STRING identifiers;
- OmniPath: both endpoints were proteins in the selected graph; and
- STITCH: one endpoint was a curated messenger represented in STITCH and the
  other was a STRING-mapped protein.

Pairs outside these definitions remain neutral (`BF=1`) because the source
could not assess the hypothesis. For dynamically added nodes, cached
localization and STRING identifiers extend the corresponding eligibility rules.
The current incremental cache does not persist a no-hit KinasePredictor
scorable-site flag, so newly added proteins without a validated site-scope flag
remain neutral under negative KinasePredictor evidence. External-target edges
use the same available source-scope checks.

Scaffold-mediated closure is excluded from negative absence. Failure to share a
strong scaffold is not an independent experimental nondetection. The rule is
also not a claim that omission from STRING, OmniPath, or STITCH proves no
interaction; it is explicitly an optional database-absence sensitivity model.

With the validated 891-node default graph, the six default primary edge streams,
all stream weights equal to one, prior 0.5, and (B^-_E=0.5), 332,215 of
396,495 unique pairs fell below the prior, 5,456 remained at the prior, and
58,824 remained above the 0.5 output cutoff. There were 1,112,985 stream–pair
penalty applications because one pair can be eligible and unsupported in more
than one source. These counts describe a sensitivity run, not a recommended
biological cutoff.

The GUI additionally exposes continuous negative evidence independently on
each primary edge stream. It removes positive-only factor floors, preserves
reported low quantitative scores as BF below 1, and represents an eligible
no-record pair as `x=0` at a small positive BF floor. It supersedes the fixed
unsupported-pair factor only for the selected stream; all source eligibility
rules above still apply. See `docs/continuous_negative_evidence.md` for the
equations and source-specific behavior.

### 11.2 Optional scaffold-mediated triadic closure

The interactive workflow includes an optional derived edge stream for
scaffold-mediated proximity. It is evaluated only after all selected primary
edge streams have been integrated. Protein nodes carrying the exact
`adaptor_scaffold` class token serve as possible common scaffold nodes; small
molecules never serve as closure endpoints or common scaffolds.

For a pre-closure protein-scaffold probability (P_{is}) and anchor cutoff
(c), proteins (i) and (j) qualify when at least one annotated scaffold (s)
satisfies (P_{is}>c) and (P_{js}>c). The default cutoff is 0.90 and is
exclusive. Every qualifying pair receives the same support likelihood 0.90,
which corresponds to (BF=0.90/0.50=1.8) relative to the workflow's neutral
likelihood. All nonqualifying pairs receive neutral (BF=1). Scaffold degree,
the number of shared scaffolds, and the amount by which an anchor exceeds the
cutoff do not alter the factor. There is no empirical (T_q), degree penalty,
continuous shared-scaffold score, or noisy-OR operation.

No inferred closure edge is fed back as a new scaffold anchor. This one-pass
restriction prevents recursive densification. The stream is disabled by
default and must be described as dependent proximity or co-complex evidence,
not as an independent experiment or direct-binding proof. With the default
891-node graph and settings, 205,920 qualifying pairs received BF 1.8, 81,323
pairs were newly raised above 0.5, and the supported-pair total increased from
169,414 to 250,737. The revised validation run and its pair-level supporting-
scaffold audit are in
`results/gui_runs/scaffold_binary_closure_validation_20260804_v2/`.

## 12. mpkCCD subcellular-fraction localization evidence

The first edge evidence stream used abundance profiles from five mpkCCD
differential-centrifugation fractions:

- 1K
- 4K
- 17K
- 200K pellet
- 200K supernatant

Repeated source entries with identical profiles were collapsed. When one gene
had distinct profiles, the profile with the greatest summed fraction abundance
was selected instead of summing isoforms.

A total of 615 of the 868 nodes had usable localization profiles.

For two observed proteins \(i\) and \(j\), similarity was the raw
five-dimensional dot product:

$$
s_{ij}=\mathbf{x}_i^\mathsf{T}\mathbf{x}_j
$$

The five-fraction vectors were not L2-normalized in this analysis.
Consequently, the score reflected both fraction-pattern similarity and overall
profile magnitude.

For each observed target node \(i\), \(T_i\) was the 75th percentile of its dot
products with every other observed node, excluding itself.

A directed localization likelihood was calculated:

$$
L_{j\rightarrow i}
=
\max\left[
0.5,\;
1-\exp\left(
-\frac{1}{2}
\left(\frac{s_{ij}}{T_i}\right)^2
\right)
\right]
$$

The two reciprocal values were averaged:

$$
L_{ij}^{\mathrm{loc}}
=
\frac{
L_{i\rightarrow j}+L_{j\rightarrow i}
}{2}
$$

This arithmetic mean made the evidence symmetric.

Pairs missing one or both profiles received the neutral likelihood 0.5.

Results:

- 188,805 pairs had profiles for both endpoints.
- 75,601 pairs were strengthened beyond 0.5.
- The maximum posterior after localization alone was \(2/3\).
- 300,677 pairs remained at 0.5.

## 13. KinasePredictor-derived kinase–protein edges

Observed-or-annotated phosphosites were scored using KinasePredictor.

Only kinase-class nodes in the 868-node universe were eligible as kinase
endpoints. The universe contained 132 kinase-class proteins, of which:

- 76 mapped to at least one KinasePredictor model.
- 79 eligible models mapped to those nodes because several genes had multiple
  model labels.

For each phosphosite:

1. The sequence was scored against all applicable official matrices:
   - 237 S/T models or
   - 98 tyrosine models.
2. \(T_q\) was the 75th percentile of the complete applicable score
   distribution before restriction to universe kinases.
3. If the raw 75th percentile was nonpositive, the 75th percentile of positive
   scores was used and flagged.
4. The top ten distinct eligible mouse kinase genes were retained.
5. Duplicate predictor labels mapping to the same gene were collapsed by
   maximum score.
6. Self-predictions were recorded for auditing but did not create self-edges.

For a raw site score \(r\):

$$
L_{\mathrm{site}}
=
\max\left[
0.5,\;
1-\exp\left(
-\frac{1}{2}
\left(
\frac{\max(r,0)}{T_q}
\right)^2
\right)
\right]
$$

$$
BF_{\mathrm{site}}=\frac{L_{\mathrm{site}}}{0.5}
$$

When multiple phosphosites supported the same unordered kinase–protein pair,
their log Bayes factors were summed:

$$
\operatorname{logit}(P_{ij}^{\mathrm{combined}})
=
\operatorname{logit}(P_{ij}^{\mathrm{localization}})
+
\sum_{\text{site hits}}
\log(BF_{\mathrm{site}})
$$

Results:

- 4,729 sites were considered.
- 4,611 were successfully scored.
- 118 were unscored.
- 45,026 site-level predictions were recorded.
- 74 self-predictions were excluded from edge construction.
- 44,952 predictions were used as edge evidence.
- 19,761 unique undirected pairs received prediction evidence.
- 18,694 pairs received non-neutral kinase evidence.
- 1,067 pairs had only neutral-floor top-ten hits.
- 3,403 supported pairs were kinase–kinase.
- 16,358 were kinase–non-kinase.

After localization and KinasePredictor integration, 91,199 of the 376,278 edges
were above 0.5.

Absence of a kinase prediction was treated as missing evidence, not evidence
against an interaction.

## 14. STRING functional-association evidence

STRING v12.0 mouse data were incorporated as an undirected
functional-association stream.

Of the 848 protein nodes:

- 847 mapped to STRING.
- `Cdk3` was the only unmapped protein.
- Small-molecule nodes were excluded from STRING evidence.

The STRING detailed network contained two orientations for each relationship.
These were collapsed to unique unordered pairs.

This produced:

- 54,922 unique STRING-supported pairs
- 10,501 pairs at or above STRING medium confidence, \(s\geq0.4\)
- 2,999 pairs at or above high confidence, \(s\geq0.7\)

STRING’s combined score \(s\) was converted to a Bayes factor using STRING’s
prior probability of 0.041:

$$
BF_{\mathrm{STRING}}
=
\frac{s/(1-s)}
{0.041/(1-0.041)}
$$

The edge posterior was updated by:

$$
\operatorname{logit}(P_{\mathrm{new}})
=
\operatorname{logit}(P_{\mathrm{previous}})
+
\log(BF_{\mathrm{STRING}})
$$

Only the STRING combined score was integrated. Individual evidence channels
were retained for auditing but were not added separately because they already
contribute to the combined score.

Missing STRING relationships were neutral.

After STRING integration:

- 122,233 edges were above 0.5.
- 254,045 remained at 0.5.

STRING relationships should be described as **functional associations**, not
necessarily direct physical PPIs.

## 15. Human Protein Atlas localization evidence

Human Protein Atlas subcellular localization data, version 25.1 and Ensembl
109, were added as another independent edge evidence stream.

Mouse proteins were mapped to human proteins using the MGI stringent one-to-one
protein-coding orthology report. No symbol-only paralog fallback was used.

Mapping results were:

- 806 mouse proteins with a one-to-one human ortholog
- 696 proteins with a unique corresponding HPA row
- 677 proteins with a primary HPA localization profile
- 482 proteins with a high-confidence profile

Primary profiles included locations assigned as:

- Enhanced
- Supported
- Approved

`Uncertain` locations were excluded from positive evidence. The separate HPA
`Extracellular location` field was retained for auditing but was not included
in the primary profile.

Each protein was represented as a binary vector across 49 HPA locations and
divided by its L2 norm. Therefore, the dot product of two profiles equaled
cosine similarity and did not automatically favor proteins with more
annotations.

For each HPA-profiled node, \(T_i\) was the 75th percentile of its cosine
similarities with every other HPA-profiled node. If this percentile was zero,
the 75th percentile of positive overlaps was used. A node with no positive
overlaps remained neutral.

Endpoint-specific likelihoods were calculated with the standard kernel and
averaged to produce an undirected HPA likelihood. The HPA Bayes factor was:

$$
BF_{\mathrm{HPA}}
=
\frac{L_{\mathrm{HPA}}}{0.5}
$$

Primary HPA evidence strengthened 31,253 pairs:

- 19,109 had previously been at 0.5.
- 12,144 were already above 0.5.

After HPA integration:

- 141,342 edges were above 0.5.
- 234,936 remained at 0.5.

An Enhanced+Supported-only profile was calculated as a sensitivity analysis. It
strengthened 14,034 pairs but was not added on top of the primary HPA stream,
which would have double-counted HPA.

The current general edge-characterization model can be written as:

$$
\operatorname{logit}(P_{ij}^{\mathrm{final}})
=
\log(BF_{ij}^{\mathrm{mpkCCD}})
+
\sum_{\mathrm{site\ hits}}\log(BF_{ij}^{\mathrm{KP}})
+
\log(BF_{ij}^{\mathrm{STRING}})
+
\log(BF_{ij}^{\mathrm{HPA}})
$$

because the initial logit of a 0.5 prior is zero.

# Part IV: Separate colocalization graph for AlphaFold screening

## 16. Purpose and distinction from the general edge graph

A separate graph was constructed specifically to ask:

> Are two proteins sufficiently colocalized—or located in physically adjacent
> compartments—to justify an AlphaFold interaction prediction?

This graph did not use KinasePredictor or STRING because those streams describe
functional relationships rather than spatial compatibility.

The colocalization graph combined:

1. mpkCCD five-fraction localization
2. HPA localization
3. mouse COMPARTMENTS localization

Only the 848 proteins were considered for AlphaFold pairs:

$$
\binom{848}{2}=359{,}128
$$

The 20 small molecules were retained in an all-node audit matrix but excluded
from PPI prediction.

## 17. COMPARTMENTS data processing

The filtered mouse COMPARTMENTS knowledge channel was used.

For each gene:

1. Annotations with confidence score 4 or 5 were retained.
2. Duplicate evidence and protein isoforms were collapsed to the maximum
   confidence score per GO cellular-component term.
3. Scores were converted to weights by dividing by 5.
4. The weighted GO-term profile was L2-normalized.
5. Pairwise cosine similarity was calculated.

A total of 747 protein nodes had a primary COMPARTMENTS profile.

As for HPA, node-specific \(T_q\) values were calculated from similarities to
all other profiled nodes. A positive-overlap fallback was used when the
all-node 75th percentile was zero.

The primary COMPARTMENTS analysis used scores 4–5. A score-3-or-higher
sensitivity analysis was also produced but not simultaneously integrated.

## 18. Initial three-source colocalization integration

For each source:

$$
BF_{\mathrm{source}}
=
\frac{L_{\mathrm{source}}}{0.5}
$$

The three source Bayes factors were combined:

$$
\operatorname{logit}
(P_{ij}^{\mathrm{colocalized}})
=
\log(BF_{ij}^{\mathrm{mpkCCD}})
+
\log(BF_{ij}^{\mathrm{HPA}})
+
\log(BF_{ij}^{\mathrm{COMPARTMENTS}})
$$

Before compartment-adjacency information was introduced:

- 132,053 protein pairs had non-neutral colocalization evidence.
- 227,075 remained at posterior 0.5.

A conservative AlphaFold filter was also defined. A pair was excluded only
when:

- the integrated colocalization posterior was neutral,
- neither HPA nor COMPARTMENTS supplied a shared discrete location, and
- at least one of those two sources had profiles for both proteins that were
  explicitly disjoint.

Unknown pairs were retained rather than interpreted as separated. Under this
conservative rule:

- 213,129 pairs were retained.
- 145,999 were excluded.

## 19. Compartment compatibility matrices

Exact compartment overlap was recognized as overly restrictive because
proteins in adjacent physical zones can interact across an interface.
Scientist-reviewable compatibility matrices, denoted \(C\), were therefore
constructed.

A canonical topology of 30 physical zones was defined, including:

- cytosol
- cytoskeleton
- plasma membrane
- extracellular space
- nucleus and nuclear envelope
- ER lumen and membrane
- Golgi lumen and membrane
- endosome lumen and membrane
- lysosome lumen and membrane
- vesicle lumen and membrane
- mitochondrial matrix, membranes, and intermembrane space
- peroxisomal matrix and membrane
- ciliary interior and membrane
- lipid-droplet surface
- broad or unresolved cytoplasmic and membrane categories

Matrix properties were:

- symmetric
- values bounded from 0 to 1
- diagonal equal to 1
- same resolved compartment given maximum compatibility
- direct physical interfaces generally assigned 0.7
- coarse or unresolved compatibility assigned lower weights
- different native labels mapping to the same resolved zone capped at 0.6

Examples included:

- cytosol–plasma membrane
- extracellular space–plasma membrane
- cytosol–ER membrane
- ER lumen–ER membrane
- cytosol–endosomal membrane
- cytosol–mitochondrial outer membrane
- mitochondrial inner membrane–matrix
- cytosol–lipid-droplet surface

Separate native-location matrices were produced for:

- 49 HPA locations
- 528 COMPARTMENTS GO cellular-component terms
- five mpkCCD fractions

The primary mpkCCD matrix was the identity matrix because centrifugation
fractions are measurement bins rather than literal biological compartments. A
separate fraction-proximity matrix was supplied only as a sensitivity analysis.

All \(C\) values were marked provisional and intended for scientist review.
They were not interpreted as probabilities of interaction.

## 20. Adjacency-aware colocalization scoring

For HPA and COMPARTMENTS, two similarities were calculated.

Exact overlap:

$$
S_{ij}^{\mathrm{exact}}
=
\cos(\mathbf{x}_i,\mathbf{x}_j)
$$

Expected compartment compatibility:

$$
S_{ij}^{C}
=
\mathbf{p}_i^\mathsf{T}C\mathbf{p}_j
$$

where \(\mathbf{p}\) was the L1-normalized location profile.

The adjacency-aware similarity was:

$$
S_{ij}^{\mathrm{adj}}
=
\max
\left(
S_{ij}^{\mathrm{exact}},
S_{ij}^{C}
\right)
$$

The maximum was used so that adding compartment adjacency could not erase
prior exact-overlap evidence.

Node-specific \(T_q\) values were recomputed from the adjacency-aware
similarities. The final source likelihood was the maximum of:

- the previously calculated exact-location likelihood, and
- the newly calculated adjacency-aware likelihood.

This made the added evidence monotonic.

The primary mpkCCD analysis continued to use the identity \(C\) matrix. The
technical fraction-proximity matrix was reported separately.

After adjacency integration:

- 167,884 pairs had non-neutral colocalization evidence.
- 191,244 remained at the neutral baseline.
- 35,831 previously neutral pairs were rescued.
- The neutral-pair count fell by 15.8%.
- 337,799 pairs passed the broader conservative AlphaFold-retention rule.
- 21,329 were conservatively excluded.

The primary supported set used later for the scaffold-restricted analysis was
the stricter **167,884 non-neutral pair set**, not the broader 337,799
conservative retention set.

# Part V: Experimental PPI prior knowledge

## 21. Experimental PPI sources and mapping

Experimentally supported mouse PPIs were collected from:

- BioGRID 5.0.259, compiled June 25, 2026
- IntAct mouse PSI-MITAB release dated January 14, 2026

Only mouse–mouse protein pairs whose endpoints mapped to the 848-protein
universe were retained. Self-interactions were excluded, and all pairs were
canonicalized as unordered pairs.

BioGRID processing retained physical interactions. IntAct processing excluded:

- negative interaction records
- non-protein endpoints
- records not containing two mouse proteins

Endpoints were mapped primarily by UniProt accession. Recorded gene aliases or
official symbols were used as controlled fallbacks.

## 22. PPI evidence tiers

Experimental records were classified into three levels.

### Tier A: direct interaction

Tier A required at least one direct or structural interaction record, defined
as:

- IntAct interaction type MI:0407, or
- BioGRID:
  - Biochemical Activity
  - Co-crystal Structure
  - Far Western
  - Protein-peptide

### Tier B: physical association

Tier B included physical association, proximity, co-complex, affinity-capture,
co-fractionation, two-hybrid, FRET, proximity-labeling, and related physical
methods when no Tier A record existed.

### Ambiguous

Records that did not satisfy Tier A or Tier B were retained as ambiguous.
Examples included association-only annotations, colocalization-only evidence,
and RNA-focused assays.

Tier A took precedence at the pair level. Therefore, “Tier B-only” meant at
least one Tier B record and no Tier A record.

Across all 359,128 protein pairs:

- 51 had Tier A evidence.
- 1,423 had Tier B-only evidence.
- 299 had ambiguous-only evidence.
- 1,474 had Tier A or Tier B evidence.

Within the 167,884 non-neutral colocalization-supported pairs:

- 33 were Tier A.
- 1,073 were Tier B-only.
- 234 were ambiguous-only.
- 1,106 had Tier A or Tier B evidence.

Absence from BioGRID or IntAct was treated as missing evidence, not evidence
that an interaction did not exist.

These tiers were used as a computational exclusion filter for AlphaFold. They
were not added to the Bayesian edge posterior as another Bayes factor.

# Part VI: Scaffold-restricted AlphaFold candidate set

## 23. Definition of scaffold proteins

The node taxonomy did not contain a literal `structural_protein` category. The
reproducible proxy used was the existing `adaptor_scaffold` class.

Among the 848 protein nodes:

- 198 were classified as `adaptor_scaffold`.
- 650 were not.

This class is a signaling-role category. It should not be presented as
identical to all cytoskeletal, extracellular-matrix, or structural proteins.

## 24. Final pair filtering

The final AlphaFold candidate list was constructed using the following
sequence:

1. Begin with all 359,128 unique protein pairs.
2. Retain the 167,884 pairs with non-neutral adjacency-aware colocalization
   evidence.
3. Remove Tier A and Tier B experimentally supported PPIs.
4. Require at least one endpoint to have the `adaptor_scaffold` class.
5. Retain ambiguous-only prior records but mark them.

The calculations were:

| Filtering stage | Unique undirected pairs |
|---|---:|
| All 848-protein pairs | 359,128 |
| Non-neutral adjacency-aware colocalization | 167,884 |
| After removing all colocalized Tier A/B pairs | 166,778 |
| Scaffold-involving before Tier A/B removal | 67,300 |
| Tier A/B scaffold pairs removed | 597 |
| Final pairs with at least one scaffold | **66,703** |

The final 66,703 pairs comprised:

- 58,265 scaffold–other pairs
- 8,438 scaffold–scaffold pairs
- 80 pairs with ambiguous-only prior evidence retained

Independent validation confirmed:

- 66,703 unique unordered pair identifiers
- no reversed duplicates
- no self-pairs
- every pair had at least one scaffold endpoint
- no remaining Tier A/B overlap
- no missing expected pair
- no unexpected pair

The pair list is in the `Pair Data` sheet of
`outputs/scaffold_pairs_848/scaffold_involving_pairs_66703.xlsx`.

## 25. AlphaFold runtime projection

Runtime projections used the empirically supplied model:

$$
\text{tokens}=1.2\times
(\text{combined amino-acid length})
$$

$$
t_{\mathrm{folding}}
=
1.29\times
\text{tokens}^{0.95}
$$

A startup cost of 73 seconds per submitted job was added. The primary batching
projection used ten protein pairs per job.

For the 66,703 scaffold-involving pairs:

- estimated folding work: 27,298 GPU-hours
- serial benchmark-configuration time: approximately 1,143 days
- ideal four-GPU wall time: approximately 285.8 days

A one-prediction-per-model approximation divided the folding component by five
while retaining startup overhead:

- serial estimate: approximately 233.1 days
- ideal four-GPU estimate: approximately 58.3 days

These are planning estimates, not measured completion times. The model fit had
\(r^2=0.89\), assumes ideal scheduling, and may not capture memory limitations,
queueing, failures, or unusually long protein pairs.

# Part VII: Ontology-based partial directionality

## 26. Separation of edge existence from signal-flow direction

Bayesian edge characterization continued to estimate one undirected association
probability for each unordered node pair. Directionality was added only as a
path-traversal layer and did not change that probability. For an edge with
probability (p_{ij}), a uniquely supported ontology direction retained
(p_{ij}) in the allowed direction and assigned zero to the disallowed reverse
traversal. Unresolved pairs retained (p_{ij}) in both directions.

This partial orientation prevents a path search from crossing a resolved
collider (A\rightarrow B\leftarrow C) from (A) through (B) to (C), or a
resolved fork (A\leftarrow B\rightarrow C) from (A) through (B) to (C).
It cannot eliminate these motifs when one or both constituent edges remain
unresolved.

## 27. Conservative ontology rule catalog

Ten rules were defined in the versioned file
`code/path_finding/ontology_direction_rules.json`:

- ligand to receptor;
- receptor regulator to receptor;
- kinase to kinase/phosphatase-binding protein;
- phosphatase to kinase/phosphatase-binding protein;
- GTPase regulator to GTPase;
- cyclase to second messenger;
- phospholipase to second messenger;
- nitric-oxide synthase to second messenger;
- phosphodiesterase to second messenger; and
- second messenger to second-messenger-binding protein.

The kinase-binding rule implements the explicit project policy that a protein
annotated only as a kinase/phosphatase binder should not propagate toward a
connected kinase. It is a modeling heuristic, not a causal relationship
entailed by the GO binding annotation. The phosphatase rule applies the same
policy symmetrically. Broad `signaling_process`, `signaling_regulation`, and
`adaptor_scaffold` labels did not orient an edge by themselves.

Every class on both endpoints was considered. One or more rules supporting
only one direction yielded a uniquely oriented edge. If multi-role annotations
supported both directions, the edge was marked
`unresolved_conflicting_rules` and retained in both directions. Pairs with no
matching rule were marked `unresolved_no_matching_rule` and also retained in
both directions. This fail-open policy avoided inventing causal order when the
ontology labels were insufficient or contradictory.

The implementation wrote four additional artifacts per directional run:

- `propagation_adjacency_matrix.tsv`;
- `ontology_directionality_audit.tsv.gz`;
- `ontology_class_pair_catalog.tsv`; and
- `ontology_direction_rules.json`.

The complete class-pair catalog contained all 153 unordered combinations with
replacement of the 17 ontology classes. The edge-level audit recorded both
nodes, their complete class sets, the matched rules in each direction, the two
allowed-traversal flags, and the final resolution status.

## 28. Initial directionality coverage

The initial validation used the default 891 selected nodes plus external target
Aqp2, giving 892 nodes and 169,769 undirected edges strictly above probability
0.5. Of these:

- 24,679 edges were uniquely oriented (14.5368%);
- 144,027 had no matching rule (84.8370%); and
- 1,063 had conflicting multi-role rules (0.6261%).

Partial orientation therefore removed 24,679 reverse traversals and left
314,859 allowed directed traversals. The largest contributor was the
kinase-to-kinase/phosphatase-binding rule. All 25 validation paths followed
allowed transitions; no path edge used a direction marked as disallowed.

Biochemical sign remains separate. The phosphodiesterase rule records a
negative-effect interpretation in the catalog, but activation and inhibition
do not yet affect path scoring. Directed OmniPath annotations are likewise not
yet combined with this ontology-only layer.

# Part VIII: Temporal validation of ranked paths

## 29. Scope and input

Temporal validation was implemented as an optional post-path annotation. It did
not contribute a node or edge Bayes factor and did not change the primary path
rank. The GUI source was `data/phospho_data_original.xlsx`, a rat inner-medullary
collecting-duct dDAVP/vehicle phosphoproteomic time course with raw sheets at 1,
2, 5, and 15 minutes and three biological replicates. For each replicate, a
log2 treated/control ratio was calculated only when both reporter intensities
were positive. Missing measurements remained missing.

Sheets were merged by UniProt accession and phosphosite position. Records with
no usable gene symbol or protein identifier were excluded. Path-symbol matching
was case-insensitive. The temporal layer did not perform a new rat-to-mouse
orthology mapping and therefore retained that limitation explicitly.

## 30. Replicate variance moderation and representative sites

Raw replicate variance at a site/timepoint was stabilized toward an empirical
intensity-dependent prior. Median variance was calculated in quantile bins of
log10 median reporter intensity and linearly interpolated. With raw variance
\(s^2\), \(n\) replicates, residual degrees of freedom \(d=n-1\), prior variance
\(s_0^2(I)\), and configured prior degrees of freedom \(d_0\), the moderated
variance was:

\[
\widetilde{s}^2 = \frac{d_0s_0^2(I)+ds^2}{d_0+d}.
\]

The moderated standard error was
\(\sqrt{\widetilde{s}^2/n}\). A two-sided Student-t p-value used
\(d_0+d\) degrees of freedom. Each gene was represented by the phosphosite with
the largest peak moderated absolute t statistic across its timepoints.

Because selecting among multiple sites and timepoints inflates an unadjusted
minimum p-value, the GUI default used a within-gene Bonferroni correction:

\[
p_{adjusted}=\min(1,m_gp_{raw}),
\]

where \(m_g\) was the number of tested site/timepoint combinations for that
gene. A gene entered path scoring when \(p_{adjusted}<\alpha\), with default
\(\alpha=0.05\). An unadjusted option was retained only to reproduce the
contributed prototype. The correction addresses the within-gene search but is
not a study-wide path or gene false-discovery-rate procedure.

## 31. Response-time uncertainty and path scores

At each timepoint, response draws were sampled from a normal distribution with
the observed mean and moderated standard error. For every draw, absolute LFC
response time was calculated as:

\[
\tau=\frac{\sum_t t|x_t|}{\sum_t |x_t|}.
\]

The default used 2,000 draws and random seed 0. For every ordered pair of scored
nodes \(i<j\), soft precedence was \(P(\tau_i<\tau_j)\). The path-level soft-
precedence score averaged these probabilities. Draw-wise Kendall tau-a supplied
a mean and 2.5/97.5 percentile interval. These were described as parametric
Monte-Carlo uncertainty intervals, not as posterior credible intervals from a
full longitudinal model.

Paths with fewer than two significant measured nodes had no temporal order
statistic. A separate `temporal_evidence_rank` was assigned only at or above the
configured coverage threshold (default three scored nodes). That secondary rank
used the lower Kendall bound, mean Kendall, soft precedence, and Bayesian path
probability as successive sort keys. The original Bayesian row order and `rank`
column were preserved.

Temporal runs wrote `ranked_paths_temporal.tsv`,
`temporal_gene_responses.tsv.gz`, `temporal_variance_trend.tsv`, and
`temporal_validation_summary.json`. Detailed methods and limitations are in
`docs/temporal_path_validation.md`.

## 32. Initial temporal coverage

The initial integration validation used the current ontology-directed default
Prkaca-to-Aqp2 graph. It loaded 6,755 valid phosphosite trajectories and
assigned representative trajectories to 2,661 genes; 237 genes passed the
within-gene-corrected 0.05 gate. All 25 Bayesian paths and ranks were preserved.
Seven paths had at least two scored nodes, but none reached the default three-
node coverage threshold. Therefore no path received an interpretable temporal
evidence rank. This was treated as insufficient temporal coverage rather than
negative evidence against those paths.

# Part IX: Reproducibility and quality control

The project preserved:

- unchanged raw workbooks
- downloaded reference datasets
- processed mapping tables
- site-level and pair-level audits
- source-specific likelihood matrices
- log-Bayes-factor matrices
- final posterior matrices
- SHA-256 hashes for major inputs and outputs
- JSON analysis summaries
- symmetric matrix checks
- diagonal checks
- probability-bound checks
- exact reconstruction checks

All reported edge counts refer to unique unordered pairs. Square adjacency
matrices contain symmetric entries for both \((i,j)\) and \((j,i)\), but pair
totals do not double-count these entries.

Two backend workbooks were also constructed to support a future interactive
application:

- a node-selection evidence-factor catalog
- an edge-characterization evidence-factor catalog

The node catalog covers 9,170 protein candidates plus 20 curated molecules and
stores the source scores needed to reproduce node evidence. Historical catalog
posterior columns may reflect the former normalized-relative-weight model; the
current workbench recomputes independent binary posteriors from the source
scores and their declared neutral values.

The edge catalog covered all 376,278 unordered pairs and stored the separate
localization, KinasePredictor, STRING, and HPA evidence streams. It
reconstructed the final edge posterior with maximum absolute error
approximately \(1.5\times10^{-9}\).

# Reporting cautions

For scientific accuracy:

- Current node posteriors are independent binary-model probabilities. Do not
  claim they are empirically calibrated frequencies; the 0.5 prior is an
  objective symmetry choice and the evidence mappings remain model choices.
- Do not say “nonzero posterior.” Neutral candidates remain at 0.5; the default
  selection criterion is posterior strictly greater than 0.5.
- Do not say missing data proves absence. Missing node or edge evidence is
  neutral in the reproducibility default. The optional negative-evidence modes
  penalize only source-eligible nondetections/unsupported pairs and must be
  reported as sensitivity assumptions.
- Do not describe STRING relationships as experimentally demonstrated physical
  PPIs. STRING supplies functional associations.
- Do not describe every phosphosite used for kinase–protein edges as observed
  in the PKA-knockout experiment. The table included UniProt-annotated sites as
  well.
- Do not claim the kinase and phosphoprotein streams were independent. They
  came from the same PKA-knockout phosphoproteomic dataset.
- Do not call the \(C\)-matrix values interaction probabilities. They were
  provisional compartment-compatibility weights.
- Do not call the `adaptor_scaffold` set a comprehensive structural-protein
  set.
- Do not say the 66,703 list came from the broader 337,799 conservative
  AlphaFold set. It came from the stricter 167,884 non-neutral colocalization
  set.
- Do not claim that Tier A/B prior knowledge was incorporated as Bayesian edge
  evidence. It was used to remove already-supported pairs from the
  structural-prediction queue.
- Do not describe the Bayesian edge probabilities as directed. The edge matrix
  and pair list remain undirected; only the separate ontology propagation
  matrix is partially directed. Unresolved pairs remain traversable both ways.
- Do not describe 891 nodes as a current implementation limit. It is the
  validated immutable seed; later non-neutral proteins and external targets
  are characterized incrementally and cached.
- Do not claim incremental KinasePredictor pairs use the historical
  seed-eligible top-ten rule. For cache stability, incremental pairs use global
  top-ten model rank for each phosphosite.
- Do not describe the final OmniPath stream as directed or signed. OmniPath
  direction and stimulation/inhibition annotations were retained in audit
  columns but intentionally ignored when assigning the undirected factor.
- Do not imply that small-molecule nodes were scored by the protein-only STRING
  network or protein expression/localization streams. They were added by
  curated rule; 14 of the 20 now receive a separate STITCH protein–chemical
  association stream. The other six remain neutral for that stream.
- Do not say temporal validation updated an edge Bayes factor, changed the
  canonical Bayesian rank, or proved causal transmission. Unrankable paths had
  insufficient significant measured nodes and were not treated as negative.
- Do not call the temporal Kendall percentile range a fully Bayesian credible
  interval. It is a conditional parametric Monte-Carlo interval under the
  moderated normal approximation.

# Key project artifacts

- Current node universe:
  `data/node_selection/node_universe_combined_nonzero.tsv`
- Combined node-selection summary:
  `results/combined_kinase_phosphoprotein_evidence/analysis_summary.json`
- Current sequential edge-characterization summary:
  `results/edge_characterization/localization_kinase_predictor_string_hpa_omnipath_stitch/analysis_summary.json`
- External-target adjacency-vector method:
  `docs/path_finding_target_extension.md`
- Aqp2 target-extension validation summary:
  `results/path_finding/target_extensions/Aqp2/analysis_summary.json`
- Signal-relay-constrained Prkaca-to-Aqp2 path summary:
  `results/path_finding/ranked_paths/Prkaca_to_Aqp2_signal_relay/analysis_summary.json`
- Ontology-directionality method and policy:
  `docs/ontology_directionality.md`
- Versioned ontology-direction rule catalog:
  `code/path_finding/ontology_direction_rules.json`
- Initial ontology-direction validation summary:
  `results/gui_runs/ontology_directionality_validation_20260812/analysis_summary.json`
- Temporal path-validation method:
  `docs/temporal_path_validation.md`
- Initial temporal integration summary:
  `results/gui_runs/temporal_integration_validation_20260814_v3/temporal_validation_summary.json`
- Configurable local analysis workbench:
  `gui/README.md`
- Extensible evidence-stream registry:
  `gui/evidence_registry.json`
- Persistent incremental edge-pair database:
  `data/edge_characterization/incremental_edge_cache/edge_pair_cache.sqlite3`
- Incremental cache source manifest:
  `data/edge_characterization/incremental_edge_cache/cache_manifest.json`
- Adjacency-aware colocalization summary:
  `results/colocalization_adjacency_aware/analysis_summary.json`
- Experimental PPI tier summary:
  `results/experimental_ppi_tiers_848/analysis_summary.json`
- Scaffold pair-filter summary:
  `results/alphafold_scaffold_pair_filter_848/analysis_summary.json`
- Final scaffold-involving pair workbook:
  `outputs/scaffold_pairs_848/scaffold_involving_pairs_66703.xlsx`
