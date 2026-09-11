# Selective PKA catalytic-subunit knockout node evidence

## Purpose

The existing node-selection workflow contains kinase-activity and
phosphoprotein-response evidence derived from a PKA double-knockout experiment.
If PKA-Cα and PKA-Cβ have opposing or partially compensatory effects, deleting
both catalytic subunits can conceal responses that are visible when either
subunit is deleted separately. The workbook `Phosphopeptides-PKAsKO.xlsx`
provides separate PKA-Cα-null/intact and PKA-Cβ-null/intact comparisons from
mouse mpkCCD cells. They are therefore represented as two independently
selectable node-evidence streams.

The new streams are disabled in the generic GUI default. This preserves the
validated 871-protein/20-messenger baseline and prevents the shared experiment
from being silently counted in every historical run. A user can enable either
comparison alone or both comparisons together.

## Source data

The workbook contains 4,635 phosphopeptide rows. The relevant columns are:

- UniProt accession;
- mouse gene symbol;
- protein annotation;
- phosphosite coordinate or coordinate set;
- annotated phosphopeptide sequence;
- log2(PKA-Cα null / PKA-Cα intact) and its modified P value;
- log2(PKA-Cβ null / PKA-Cβ intact) and its modified P value.

The PKA-Cα-null/PKA-Cβ-null comparison in the final workbook column is not used
as a node-evidence stream because it does not compare a knockout with its
matched intact control.

The readme sheet states that CRISPR-Cas9 was used to mutate `Prkaca` or
`Prkacb` in cultured mouse mpkCCD cells and that TMT phosphoproteomics used
three biological replicates.

## Duplicate-site collapse

Several workbook rows describe the same UniProt accession and site-pattern
key, often using overlapping peptide sequences. Treating these rows as
independent sites would overstate the number of opportunities for an extreme
protein-level value. Rows are therefore grouped by:

```text
UniProt accession | normalized Site(s) value
```

Within each key and comparison, the signed log2 fold changes are collapsed by
their median. Pmod values are also collapsed by their median for audit. This
reduced 4,635 source rows to 3,805 unique site-pattern records and collapsed
830 duplicate rows.

Pmod is not used to filter records or calculate the node Bayes factor. It is
retained for scientific review and used only as a deterministic tie-breaker
when two site patterns have the same absolute LFC.

## Mapping to the signaling-candidate universe

Mapping is performed independently of the quantitative score:

1. Match the source UniProt accession to the UniProt accession stored in the
   9,170-protein mouse signaling-candidate catalog.
2. If no accession match exists, split the supplied gene-symbol field on
   semicolons and retain exact symbols present in the candidate catalog.
3. Preserve every valid mapped candidate rather than inventing a mapping for an
   unresolved record.

Of the 3,805 unique source site patterns, 1,680 mapped to the signaling
candidate universe: 1,664 by UniProt and 16 by exact gene-symbol fallback.
The remaining 2,125 records remain in the source/mapping audit but do not
contribute candidate-specific node evidence.

## Quantitative background

Each knockout comparison receives its own empirical background. The
background is formed from the absolute median LFC of all 3,805 finite,
duplicate-collapsed source site patterns before signaling-universe filtering.
Consequently, the threshold is not estimated only from signaling proteins.

The ordinary 75th-percentile background values are:

| Comparison | Background records | Absolute-LFC q75 |
|---|---:|---:|
| PKA-Cα null / intact | 3,805 | 0.22603867 |
| PKA-Cβ null / intact | 3,805 | 0.18371188 |

The two backgrounds are not pooled. This lets each comparison retain its own
measurement scale.

## Site-level scoring and protein hit rule

For each mapped signaling protein and each comparison:

1. Assign every duplicate-collapsed site pattern the comparison's single
   source-background q75, denoted by \(T_q\).
2. Calculate its complement-of-minimum likelihood:

   $$
   L_i=\max\left[0.5,
   1-\exp\left(-\frac{1}{2}
   \left(\frac{x_{site}}{T_q}\right)^2\right)\right].
   $$

3. Convert the site likelihood to a neutral-relative Bayes factor:

   $$
   BF_{site}=\frac{L_{site}}{0.5}.
   $$

4. Assign each protein the maximum site-level BF among its detected site
   patterns. Therefore at least one site with BF greater than 1 is sufficient
   to make the protein a positive hit. Additional sites do not raise its
   threshold and their BFs are not multiplied. The detected-site count remains
   in the audit but has no role in scoring.

An observed low-response protein has BF 1 rather than negative evidence. A
protein without a mapped site pattern also has BF 1 unless the GUI's global
optional nondetection penalty is enabled. The stream-specific `Tq ×` control
multiplies the one site-level \(T_q\) before every site is rescored. A lower multiplier is more
sensitive; a larger multiplier is more conservative.

## Evidence counts

Both comparisons observed the same 776 signaling candidates because they are
columns in the same phosphopeptide table. Their quantitative evidence differs:

| Stream | Observed candidates | BF > 1 | Observed at BF = 1 |
|---|---:|---:|---:|
| PKA-Cα KO phosphoprotein response | 776 | 267 | 509 |
| PKA-Cβ KO phosphoprotein response | 776 | 261 | 515 |

One hundred thirty-one candidates have BF greater than 1 in both streams. The remaining
positively supported candidates are comparison-specific, which is the reason
the two columns are not collapsed into a single absolute-LFC statistic.

With the four historical default node streams plus both new streams enabled at
unit weights and `Tq × = 1`, the node-selection result contains 1,232 proteins.
Adding the 20 curated secondary messengers gives 1,252 nodes. With the two new
streams off, the revised site-centric double-KO default contains 1,051 proteins
plus 20 messengers.

## Interpretation of opposing effects

The node-selection hypothesis is whether a protein participates in the
signaling system, not whether its phosphorylation increases or decreases.
Accordingly, each new stream scores absolute LFC. A large increase and a large
decrease both support node inclusion. Because the streams are separate, an
alpha-specific decrease is not numerically canceled by a beta-specific
increase.

The signed selected LFC remains available in the audit. The current node model
does not convert that sign into activation, inhibition, or causal edge
direction. It also does not infer a new kinase-activity profile from these two
workbook columns; these additions are protein-level differential-
phosphorylation streams. Any future KinasePredictor-derived alpha- and
beta-specific kinase-activity streams should be registered separately so their
statistical meaning remains explicit.

## Dependence and limitations

The alpha- and beta-knockout comparisons share the same TMT study, cell model,
platform, and phosphopeptide-detection process. They are tagged with the shared
dependence group `pka_subunit_ko_phosphoproteomics`. Multiplying their Bayes
factors is allowed for a sensitivity-oriented analysis, but the resulting
posterior should not be described as if the two streams were independent
experiments.

Additional limitations are:

- multi-site phosphopeptides are retained as site-pattern records and are not
  interpreted as unambiguous single residues;
- the maximum-absolute-LFC statistic favors the strongest detected response;
- missingness can reflect assay sensitivity rather than true absence;
- Pmod is audited but does not control a false-discovery rate in this evidence
  transformation; and
- the BF scale is a project-specific evidential transformation, not an
  externally calibrated frequency of true signaling participation.

## Reproducibility files

Preprocessing entry point:

```text
code/node_selection/build_pka_subunit_ko_evidence.py
```

Raw source:

```text
data/node_selection/pka_subunit_ko/raw/Phosphopeptides-PKAsKO.xlsx
```

Generated outputs:

```text
data/node_selection/pka_subunit_ko/processed/pka_subunit_ko_node_factors.tsv.gz
data/node_selection/pka_subunit_ko/processed/pka_subunit_ko_site_audit.tsv.gz
data/node_selection/pka_subunit_ko/processed/pka_subunit_ko_site_level_thresholds.tsv
data/node_selection/pka_subunit_ko/processed/pka_subunit_ko_analysis_summary.json
```

The GUI stream definitions and exact factor-column mappings are stored in
`gui/evidence_registry.json`. Every completed GUI run also records the active
streams, weights, `Tq` multipliers, observation masks, source BFs, effective
BFs, weighted log-BFs, and posterior probabilities.
