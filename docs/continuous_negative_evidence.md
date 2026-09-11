# Continuous negative evidence

## Purpose

The GUI now offers an optional **Allow continuous negative evidence** switch on
every quantitative node stream and every primary quantitative edge stream. It
is disabled by default so existing analyses remain reproducible. It can be
enabled independently for individual datasets.

The historical positive-or-neutral scoring rule was:

$$
L_+(x)=\max\left(0.5,1-\exp\left[-\frac12(x/T_q)^2\right]\right).
$$

It made weak measurements and nondetections neutral. The continuous option
removes that 0.5 lower bound:

$$
L_c(x)=1-\exp\left[-\frac12(x/T_q)^2\right]
$$

and calculates the operational Bayes factor relative to the same neutral
likelihood of 0.5:

$$
BF_c(x)=\max\left(\epsilon_{BF},\frac{L_c(x)}{0.5}\right).
$$

Eligible nondetections or eligible no-record pairs are inserted as `x = 0`.
They therefore evaluate to the user-set numerical floor
$\epsilon_{BF}$ rather than receiving a separate categorical nondetection
factor. A small positive floor is necessary because the integration is carried
out in log-odds space and `log(0)` is undefined. The GUI exposes one floor for
node streams and one for edge streams; both default to `1e-6`.

This kernel has three important properties:

- `BF < 1` is evidence against the node/edge-present hypothesis;
- `BF = 1` occurs when $x=T_q\sqrt{2\ln 2}\approx1.177T_q$; and
- `BF` approaches 2 for very large Gaussian-kernel measurements.

Thus, $T_q$ is a scale parameter, not the exact BF-neutral crossing point.

## Treatment of Tq

The positive empirical Tq already stored for each dataset remains the
denominator. Nondetections are included as `x = 0` in the evaluated evidence
vector and hence in the resulting BF/posterior distribution, but they are not
used to recompute Tq. This is deliberate. For sparse streams, zero-padding the
entire 9,170-candidate universe makes the 75th percentile exactly zero, which
would make `x/Tq` undefined. Retaining the measured positive background keeps
the scale interpretable while still scoring nondetection through the same
continuous kernel as low detection.

The existing per-dataset **Tq ×** (or **Ref ×**) control continues to adjust
that positive reference. The continuous option therefore changes the lower
tail of the evidence model, not the definition of the source's measured
background.

## Node-stream scope

The option is available for all current node streams:

- mpkCCD protein abundance;
- principal-cell RNA;
- PKA double-KO kinase activity;
- PKA double-KO phosphoprotein response;
- selective PKA-Cα and PKA-Cβ KO phosphoprotein responses; and
- CCD, OMCD, and IMCD proteome and RNA streams.

For abundance and phosphoprotein-response streams, every modeled protein is in
scope and a missing/zero value is evaluated as `x = 0`. For kinase activity,
only proteins carrying the exact `kinase` ontology class are eligible;
non-kinases remain at `BF = 1`. Curated small molecules are appended after
protein node scoring and are not assigned protein-assay negative evidence.

When the older stage-wide fixed nondetection rule and a per-stream continuous
switch are both enabled, the continuous rule supersedes the fixed rule only for
that stream. Other active streams may still use the fixed rule.

## Edge-stream scope

The option is available for every primary edge stream:

- mpkCCD localization;
- KinasePredictor;
- STRING;
- HPA primary or HPA high-confidence localization;
- OmniPath core; and
- STITCH messenger–protein associations.

Measured low source values are retained as factors below 1. An eligible pair
with no source record is assigned `x = 0` and receives the BF floor. Eligibility
is still source-specific: both required localization profiles, a kinase and a
scorable phosphoprotein, two STRING-mapped proteins, two protein endpoints for
OmniPath, or a represented messenger and STRING-mapped protein for STITCH.
Missing coverage and out-of-scope endpoint types remain neutral at `BF = 1`.

STRING and STITCH use score odds relative to their positive reference rather
than the Gaussian kernel. Their continuous option removes the old `BF >= 1`
clamp for reported scores and uses the same configured floor for eligible
zero/no-record pairs. KinasePredictor still combines all retained site-level
top-ten hits for a kinase–protein pair; the final stream BF is bounded below by
the configured floor.

Scaffold-mediated closure is intentionally excluded. It is a derived binary
rule, and failure to share a strong scaffold is not an independent quantitative
measurement of noninteraction.

## Bayesian integration and audit fields

No new updating equation is introduced. Each effective stream factor enters
the existing independent binary update:

$$
\operatorname{logit}(P(H\mid D))=
\operatorname{logit}(P(H))+\sum_k w_k\log(BF_k).
$$

The node audit records each source BF, whether the node was eligible and
observed, and whether continuous negative evidence was applied. Node and edge
run summaries list the continuous-enabled stream IDs, the numerical floor, the
number of factors below 1, and the number of eligible zero/no-record
imputations. This makes the sensitivity setting recoverable from every run.

