# Version 2 complete, resumable pipeline

## Scope

Version 2 no longer requires a manually prepared node universe from Version 1.
A single saved run now owns the complete sequence:

1. Normalize the submitted Bayesian configuration.
2. Calculate node posteriors from every enabled node evidence stream.
3. Select nodes at the configured posterior cutoff and optionally add second
   messengers.
4. Resolve the requested target and receptor, appending either endpoint when
   it is outside the selected signaling universe.
5. Calculate the configured inexpensive edge evidence over those graph nodes.
6. Export the exact Bayes factor contributed by every enabled edge stream.
7. Build ontology and optional OmniPath traversal constraints.
8. Map graph proteins to UniProt sequences.
9. Initialize the bounded backward frontier at the target.
10. Enumerate candidate upstream parents, update their edge probabilities with
    cheap evidence, and retain the configured shortlist.
11. Reuse cached AlphaPulldown scores, accept verified HuRI positives as
    nonredundant structural substitutes when configured, and package only the
    remaining new shortlisted protein pairs for Biowulf.
12. Import completed ipTM scores, update edge probabilities, prune the frontier,
    and repeat until the receptor is reached or maximum depth is exhausted.

The node posterior is a selection gate; it is not multiplied into an edge or
path score. Edge hypotheses start independently at the configured edge prior.
A path score remains the product of its final edge posteriors (stored as a sum
of logarithms for numerical stability).

## Bayesian handoff

The workbench's normal evidence engine performs node and cheap-edge updates:

```text
posterior log odds = prior log odds + sum_s weight_s * log(BF_s)
```

This preserves each selected stream's Tq/reference multiplier, weight,
scope-aware missingness rule, and optional continuous negative-evidence rule.
For the frontier engine, the complete-pipeline initializer writes one sparse
table per enabled cheap edge stream. It derives the unweighted BF from the
exact weighted log contribution produced by the workbench, then records the
original stream weight in the search configuration. Applying the weight once
inside the frontier engine therefore reconstructs the same contribution.
Pairs omitted from those exported tables are neutral because all applicable
negative factors have already been materialized during the full cheap-edge
calculation.

AlphaPulldown is different: it is never evaluated over the full graph. Its raw
ipTM score is converted to a bounded odds ratio relative to the configurable
structural reference score and is applied only to shortlisted frontier pairs.
The raw score, transformation settings, BF, weighted contribution, and final
edge probability are retained in the round ledgers.

HuRI substitution is deliberately non-additive. For a reported mapped HuRI
positive, the HuRI BF enters the cheap-evidence sum once. When
`substitute_for_structural` is enabled and the HuRI weight is positive, that
pair is omitted from `pairs.tsv`; its structural status is
`substituted_by:huri_binary_interaction` and its structural BF is 1. HuRI
nonreporting or a scoped weak-negative factor does not substitute for an
AlphaPulldown result.

Scaffold-mediated closure is evaluated separately from the full cheap-edge
matrix. A candidate pair receives closure support only if both endpoints have
qualifying physical interactions with the same exact `adaptor_scaffold` node:
an accepted cached AlphaPulldown/AlphaFold score strictly above `Anchor >`, or
a reported-positive HuRI pair. General edge posteriors and nonphysical streams
are never anchors. The fixed closure likelihood is converted to BF relative to
0.5 (default `0.90 / 0.50 = 1.8`). The pair-level audit retains the common
scaffold and source/protocol of both physical anchors.

## GUI

Launch the repository and open `/v2.html`:

```text
python launch.py
```

The page exposes all registered node and edge streams, per-stream weights,
per-stream Tq/reference multipliers, preferred Tq values used by calibration,
continuous negative-evidence switches, global priors/cutoffs, allowed relay
ontologies, target/receptor, frontier sizes, and structural parameters.

`Initialize complete run` performs steps 1–9. `Run to checkpoint` performs all
currently possible frontier work and pauses only when AlphaPulldown input is
needed or the search is complete. `One decision` is the slower development
control. On Biowulf, `Submit AlphaPulldown` submits the current feature and GPU
arrays. When those jobs finish, `Collect + continue` imports their scores and
runs to the next structural checkpoint or completion. Every operation updates
the same atomic run state. Restarting the web server does not invalidate a run;
the run can be recovered by its directory/identifier.

The trace, selected node list, submitted configuration, and pipeline state are
downloadable from the run panel. Complete run directories are written under
`runs/` and are excluded from Git. Writable incremental and structural pair
caches are stored under Version 2's ignored `runtime/` directory, so a linked
Version 1 evidence archive is not used as Version 2's cache destination.

During development runs, the run panel polls a read-only live frontier endpoint.
It draws the top 20 ranked cheap candidates connected to each current source,
shows the exact BF, weight, and weighted log-BF contributions on selection,
and reads completed ipTM CSVs even before formal collection. The displayed
combined posterior is provisional until collection. After frontier selection,
the same graph labels candidates as retained or culled using the final
combined-evidence decision.

The current Version 2 evidence defaults match Version 1 except for HPA
localization. Both HPA alternatives default to `Tq × = 0.5`; lowering the
threshold is more sensitive under the complement kernel because the same
similarity produces a larger `x/Tq` ratio and more pairs receive BF above 1.
Their calibration-preferred multiplier is also 0.5.

## CLI equivalents

The GUI is optional. The exact same engine is available as:

```text
python dynamic_search.py pipeline-config --output configs/my_pipeline.json
python dynamic_search.py pipeline-init --config configs/my_pipeline.json --run-dir runs/my_run
python dynamic_search.py pipeline-step --run-dir runs/my_run --until-checkpoint
python dynamic_search.py pipeline-submit --run-dir runs/my_run
python dynamic_search.py pipeline-collect --run-dir runs/my_run
python dynamic_search.py pipeline-status --run-dir runs/my_run
python dynamic_search.py pipeline-trace --run-dir runs/my_run
```

The first two commands replace the former requirement to run Version 1 node
selection manually. The submit command intentionally fails outside a Slurm
host; local development can still inspect the exact round package.

## Data layout and Biowulf portability

The source repository does not duplicate the multi-gigabyte companion data
archive. For local development, Version 2 can use directory links to the
existing Version 1 `data/` and `results/` trees. These links are ignored by Git;
all Version 2 run outputs and its structural pair cache remain under Version 2.

For Biowulf, copy or extract the companion `data/` and `results/` directories
inside the repository. No directory link and no AI service is required. The
web GUI is a standard-library HTTP server and the batch scripts are ordinary
deterministic Bash/Python/Slurm programs. If a browser cannot connect directly
to a compute/login node, use an NIH Open OnDemand interactive session or an SSH
port forward to the host and port printed by `python launch.py --foreground
--no-browser`.

## Important checkpoint behavior

- Missing sequence mappings do not reject candidates during the cheap screen.
  Sequence availability is validated only after shortlist selection and only
  for pairs that actually require AlphaPulldown.
- A missing AlphaPulldown result is never imputed. Collection remains incomplete
  until every pair in the immutable round manifest has a valid analyzed ipTM
  result.
- The structural cache key includes the unordered pair, protocol identifier,
  and metric, so changing a protocol does not silently reuse an incompatible
  score.
- Directionality restricts traversal only. Bayesian edge existence remains an
  undirected hypothesis.
