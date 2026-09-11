# Lab notebook — GUI evidence visualization and calibration transparency

Date: 2026-09-10

## Objective

Make the completed-run interface easier to interpret without changing node,
edge, directionality, or path calculations. The requested changes were to use
color instead of line width for network evidence strength, connect graph
selections to the existing evidence-transparency ledger, show where each
selected stream score falls within that stream's run-specific distribution,
and display calibrated weights directly in the GUI rather than only in JSON.

## Network encoding

The merged top-path graph still visualizes the integrated node posterior and
edge posterior, not an individual source weight. Node radii and edge widths are
now fixed. Posterior distance from the neutral 0.5 prior is converted to the
existing capped absolute-log-odds visual strength and represented as a
continuous color mixture:

- increasingly green means stronger posterior support above 0.5;
- increasingly red means stronger counterevidence below 0.5;
- gray means at or very near 0.5; and
- orange means that a curated messenger or external endpoint has no Bayesian
  node-selection posterior.

Direction arrows remain separate traversal constraints. Their color follows
the edge they annotate; they do not encode activation or inhibition.

Selecting a graph node automatically fills and runs the node evidence
inspector. Selecting an edge does the same for the unordered edge hypothesis.
Selecting a ranked path highlights all graph marks used by that route and
initially opens the ledger for the route's lowest-probability edge. A path does
not receive a fabricated per-stream ledger of its own: its score remains the
product of constituent edge probabilities, so component evidence is inspected
at the node/edge level.

## Per-stream Bayes-factor distributions

Raw datasets use incompatible units. The common quantity actually entering the
Bayesian update is the applied Bayes factor, so distribution placement is
defined on BF rather than on raw abundance, localization, or database-score
units. For each active node stream, all modeled protein candidates are
included. For each active edge stream, all unique unordered graph pairs are
included. This deliberately retains BF=1 neutral assignments and any BF<1
negative-evidence assignments, so the plot represents what the stream actually
contributed to that run.

Each distribution is summarized on a symmetric log2(BF) axis with 41 bins. An
odd bin count gives BF=1 a true center bin. Counts are displayed with a
log10(count + 1) height transformation. The stored summary also contains 101
quantiles and support/neutral/refutation counts. Node percentiles are calculated
exactly from the saved per-node factor column. Edge vectors are not serialized
again: exact percentile/tie intervals are available for compact discrete
streams, while continuous edge streams use the compact quantile/histogram
summary and are explicitly labeled as estimates. Older completed runs remain
inspectable but state that a rerun is needed to add distributions.

The evidence table continues to display the applied BF, configured or fitted
stream weight, and weighted change in log2 posterior odds. The new percentile
column and selected-row histogram supplement these values; they do not replace
the arithmetic posterior reconciliation.

## Calibration presentation

The optimization and saved outputs are unchanged. When positive-control
calibration is enabled, the result page now pivots the existing parameter
records into separate node and edge tables. Each primary stream displays:

- starting weight and fitted weight;
- starting Tq/reference multiplier;
- preferred regularization multiplier; and
- fitted Tq/reference multiplier.

Fitted values indicate increases or decreases with text color, and hovering a
fitted value exposes its allowed bounds. `calibration_summary.json`, the
calibrated-parameter TSV files, and optimizer traces remain the definitive
machine-readable audit.

## Scientific invariants

No evidence score, Tq transformation, stream weight, prior, posterior cutoff,
direction rule, path constraint, or path-ranking equation was changed. The new
factor summaries are read-only audit metadata. Network color represents the
combined posterior; the inspector's `Weight` column remains the distinct
per-stream exponent in `posterior odds = prior odds × product(BF^weight)`.

## Verification

The complete GUI test suite passed after implementation: 54 tests passed and
one optional SciPy-dependent test was skipped in the bundled runtime. Focused
syntax and inspector tests were then rerun after the final presentation edits.
A live API smoke run used one node stream and one edge stream and completed for
467 selected nodes (108,811 unique unordered pairs). Its node and edge
inspectors returned the selected Bayes factor together with its run-specific
distribution and percentile position. The live `/api/config` response was also
checked after the smoke run to confirm that the normal default configuration
was not altered.
