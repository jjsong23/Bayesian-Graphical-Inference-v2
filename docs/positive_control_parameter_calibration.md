# Positive-control calibration of evidence weights and Tq scales

## Purpose and interpretation

The workbench can optionally fit evidence-stream weights and Tq/reference
multipliers to user-supplied nodes and edges that are accepted as known present.
Node parameters and edge parameters are fitted in two independent stages. This
is target fitting, not supervised binary classification: hypotheses omitted from
the supplied lists are unlabeled, not negative controls, and do not enter the
fit loss.

Calibration is disabled by default. When disabled, the workflow is numerically
identical to the manual configuration. When enabled, node calibration is run
first, the complete candidate catalog is rescored with its fitted node settings,
and node selection is repeated. Known-edge endpoints are then resolved and added
to the graph if necessary, raw evidence for new pairs is characterized once and
cached, and edge calibration is run on the selected/extended graph. The fitted
edge settings are finally used for ordinary full-matrix characterization.

## Bayesian model

For a known node *i*, the workbench evaluates

```text
eta_i = logit(node prior) + sum_k weight_k * log(BF_ik(Tq multiplier_k))
P_i   = inverse_logit(eta_i)
```

The edge stage uses the same expression for each known unordered pair, with the
edge prior and edge-stream parameters. Priors, factor definitions, optional
scope-aware factors below one, and eligibility rules are exactly the same as in
the ordinary workflow. Calibration changes only enabled primary-stream weights
and their Tq/reference multipliers.

## Objective and one regularizer

For either stage, the positive-control fit term is the mean negative log
probability:

```text
fit loss = -mean(log(P of each supplied positive control))
```

One quadratic regularizer acts as a rubber band around scientifically preferred
settings:

```text
R = lambda * mean(
      ((weight - preferred weight) / 0.5)^2,
      ((log multiplier - log preferred multiplier) / log(2))^2
    )
total loss = fit loss + R
```

There is one user-facing regularization strength, lambda, because the
standardizing denominators already put weight deviations and multiplicative Tq
deviations on interpretable scales. A weight movement of 0.5 and a two-fold Tq
movement each count as one standardized unit. The default is lambda = 0.1.

Preferred weights are 1. Every primary evidence stream has its own editable
preferred Tq/reference multiplier on its evidence card in the GUI. Ordinary
streams default to 1. The five designated phosphoproteomic streams—PKA double-KO
kinase activity, PKA double-KO phosphoprotein response, selective PKA-Cα-KO
response, selective PKA-Cβ-KO response, and edge-level KinasePredictor—default to
0.1, reflecting their intended low-sensitivity/high-specificity use. A preferred
multiplier is a regularization anchor used only during positive-control
calibration; changing it does not directly rescore an uncalibrated run. A value
of zero is not allowed because Tq appears in a denominator and because the
multiplier is optimized on a logarithmic scale.

## Bounds and optimizer

The implementation uses deterministic bounded multi-start
`scipy.optimize.minimize(method="Powell")`. Powell was selected because the
objective is low dimensional, may not be smoothly differentiable where source
records cross a neutral threshold, and does not require analytical gradients.
The current bounds are:

- evidence weights: 0 to 3;
- ordinary Tq/reference multipliers and their preferred values: 0.25 to 4;
- phosphoproteomic Tq multipliers and their preferred values: 0.05 to 1.

Tq/reference multipliers are optimized on the natural-log scale. The GUI
requests two starts: the submitted configuration and the preferred anchor
profile. Identical start vectors are collapsed, so one distinct start is run
when those profiles coincide. The run with the lowest complete regularized
objective is retained. The process is deterministic and cooperative
cancellation is checked during objective evaluations.

Derived scaffold-mediated closure is excluded from optimization. It is a
deterministic transformation of the already integrated pre-closure graph, not
an independent primary dataset. Its manual anchor and support settings remain
fixed.

## Target resolution and caching

Known node controls must resolve case-insensitively to modeled signaling
candidates in the node factor catalog. Known edges are treated as unique,
unordered pairs. An edge endpoint absent from the selected node graph is
resolved with the same mouse-universe/UniProt procedure used for external path
targets and is appended with the explicit reason `calibration_edge_endpoint`.
Any missing pair evidence is characterized and added to the persistent SQLite
cache before edge fitting. During optimization, only the named cached pairs are
queried, so the full incremental cache is not rescanned at every Powell step.

## Audit outputs

Each calibrated run writes:

- `submitted_configuration.json`: settings before optimization;
- `configuration.json`: effective fitted settings used for inference;
- `calibration_summary.json`: objective, optimizer, target, and resolution audit;
- `node_calibrated_parameters.tsv` and/or `edge_calibrated_parameters.tsv`;
- target-level probabilities before and after fitting;
- compressed objective-evaluation traces for each fitted stage.

The parameter tables contain current, independently configured preferred,
fitted, lower/upper-bound, regularization-scale, and optimization-scale fields
for every optimized value. Run configurations created before per-stream controls
were introduced remain readable: the former global phosphoproteomic preference
is migrated to any designated phosphoproteomic stream that does not already have
an explicit per-stream value.

## Limitation

A higher probability for the supplied positive controls establishes only that
the chosen settings fit those controls under this model. Without held-out
positives and appropriate negative or positive-unlabeled validation, it does not
establish specificity, population-level probability calibration, or general
predictive performance. Results must therefore be described as regularized
positive-control target fitting, not as a fully trained or validated classifier.
