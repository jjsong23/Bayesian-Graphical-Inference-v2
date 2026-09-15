# Version 2 plan: hierarchical backward frontier inference

## Objective

Version 1 estimates a broad node universe, characterizes a large undirected
graph, and then searches it for paths. Version 2 adds an alternative for costly
structural evidence: begin at one biological endpoint, expand backward toward a
specified receptor, and calculate only the edge hypotheses that can enter the
current search frontier. This is a hypothesis-prioritization workflow, not a
claim that unexamined pairs do not interact.

The initial use case is a collecting-duct endpoint such as Aqp2 or a
transcription factor and a receptor such as Avpr2. The implementation is
endpoint-agnostic.

## Search state and direction

A partial path is stored in backward order:

```text
target <- candidate parent <- candidate grandparent <- ...
```

When the current frontier node is `b`, a candidate `a` is admissible only if
the direction layer permits forward propagation `a -> b`. An undirected or
directionally unresolved association remains admissible in either orientation;
an explicitly disallowed traversal is removed before scoring. Structural
prediction contributes association strength only and never supplies arrow
direction or activation/inhibition sign.

Candidates are restricted to the user-selected signaling-relay ontology
classes. The target is exempt because it need not itself be a signaling node,
and the receptor is always admissible as the desired terminal point of the
backward search. A node already present in the same partial path cannot be
reused, so every retained path is loopless.

## Evidence hierarchy for one round

For each partial path in the current frontier:

1. Pair its current node with every eligible candidate parent.
2. Apply all configured cheap evidence streams. The first intended stream is
   HPA localization, but any Version 1 factor catalog with `node_a`, `node_b`,
   and a positive Bayes-factor column can be used.
3. Start every edge at the configured prior (default 0.5) and update its odds:

   ```text
   log posterior odds = log prior odds
                        + sum(weight_s * log(BF_s))
   ```

   Missing cheap records use each stream's explicit `missing_bayes_factor`;
   the default is neutral BF=1. This avoids treating incomplete HPA coverage as
   proof against an edge.
4. Rank candidates by cheap-evidence posterior and retain the top 20 **for each
   current frontier path**. The receptor is appended if necessary so the cheap
   screen cannot make the intended endpoint unreachable.
5. Look up these retained unordered pairs in the persistent structural cache.
   Only pairs without a score under the selected AlphaPulldown protocol and
   metric are written to a Biowulf round package.
6. Pause in `waiting_for_structural`. No score is guessed or imputed.
7. After structural results are imported, convert the raw interaction score to
   an explicitly heuristic, bounded odds ratio:

   ```text
   structural BF = odds(score) / odds(reference score)
   ```

   The default ipTM reference is 0.6, BF floor is 0.1, and ceiling is 10. Thus
   scores below 0.6 can refute the edge hypothesis, 0.6 is neutral, and scores
   above 0.6 support it. ipTM is not a calibrated probability, so both the raw
   score and this transformation are retained and the parameters must be
   sensitivity-tested.
8. Add the weighted structural log-BF to the cheap log-odds to obtain the final
   edge posterior for the round.

## Frontier pruning and path ranking

The score of a partial path is the product of its edge posteriors, implemented
as a sum of log probabilities for numerical stability:

```text
log path score = sum(log(edge posterior))
```

All extensions from all current paths are compared together. The next frontier
contains the five highest-scoring extensions with **distinct new frontier
nodes**. Global pruning to five avoids exponential `5^depth` growth while still
letting multiple current nodes compete on the same probability scale.

The default search stops when any retained extension reaches the receptor or
after six edges. The state file records all completed receptor paths seen in
that round. A future extension can continue after the first receptor hit to
collect a fixed number of alternative complete paths.

This is beam search/dynamic programming in the practical sense: pair evidence
is memoized and each round retains only a bounded set of best partial states.
It is not an exact all-paths algorithm, and a true path can be lost during the
top-20 or top-five pruning steps. Both widths therefore remain explicit
sensitivity parameters.

## Persistent records

Each run contains:

- `state.json`: atomic, resumable search state and exact configuration snapshot;
- `round_NNN/cheap_candidates.tsv`: all candidates surviving the cheap stage;
- `round_NNN/pairs.tsv`: unique uncached pairs sent for structural prediction;
- `round_NNN/sequences.fasta`: sequences needed by that batch;
- `round_NNN/scored_candidates.tsv`: cheap, structural, edge, and path scores;
- `solutions.json`: completed forward and backward receptor-target paths.

When initialized with `--development`, the same inference is divided into six
one-command stages. The run additionally contains the complete excluded-node,
per-stream update, cache lookup, and pruning-decision ledgers described in
`V2_STEPWISE_DEVELOPMENT_MODE.md`, plus an append-only event log and a refreshed
`development_trace.html`. Development mode changes observability and pause
points only; it does not use a separate scoring algorithm.

The cross-run cache is
`runtime/pair_cache.sqlite3`. Cheap factors are keyed by unordered
pair and evidence ID. Structural results are keyed by unordered pair, protocol
ID, and metric, so changing the BF reference/weight rescales the saved raw score
without rerunning AlphaPulldown, whereas changing the prediction protocol can
request a separate cached result.

## Biowulf execution boundary

The repository prepares and consumes plain files. On Biowulf, only deterministic
Python/Bash/SQLite code, Slurm, the installed AlphaPulldown module, and its
analysis helper are used. There is no call to Codex, ChatGPT, an LLM API, or an
AI assistant.

For every pending round:

1. `submit_round.sh` launches a CPU feature-generation array with at most ten
   concurrent tasks, using `--skip_existing=True` and a repository-level shared
   feature directory.
2. A dependent GPU array launches one AlphaPulldown pulldown job for each exact
   manifest pair. This avoids the unintended Cartesian products that would
   arise from one multi-bait/multi-candidate list.
3. Each GPU task requests one A100 GPU. The array concurrency is configurable
   and defaults to ten.
4. `run_get_good_pae.sh` produces an ipTM table for each pair.
5. `dynamic_search.py collect` refuses to continue if any required pair lacks a
   numeric ipTM score. Successful rows enter the persistent cache.

The exact operational commands are in `biowulf/README.md`.

## Input contract

The search configuration is JSON so it is portable and auditable. Required
inputs are:

- a node table containing a unique gene-symbol column and, when ontology
  filtering is used, semicolon-delimited classes;
- a FASTA whose first header token is the same safe symbol used in the node
  table, including the external target and receptor;
- one or more sparse cheap-evidence tables with pair columns and positive BFs;
- optionally, a direction table with `source`, `target`, and `allowed` columns.

`directionality_mode="disallowed_rows"` uses only rows explicitly marked false.
`directionality_mode="listed_pair_is_directed"` treats a listed source-target
record as the only allowed orientation for that unordered pair unless both
directions are present.

## Commands

```bash
python dynamic_search.py build-fasta \
  --nodes data/node_selection/node_universe_combined_nonzero.tsv \
  --uniprot data/edge_characterization/kinase_predictor/phosphosite_database/raw/uniprot_mouse_reference_proteome.tsv.gz \
  --output data/backward_search/node_sequences.fasta

# Required when the endpoint is outside the internal HPA pair catalog. This
# calculates only the endpoint-to-universe vector, not a square graph.
python dynamic_search.py build-target-evidence \
  --project-root . \
  --target Aqp2 \
  --base-evidence results/backend_bayes_factor_catalogs/edge_factors_891/hpa_primary_bf_gt1.tsv.gz \
  --output data/backward_search/hpa_primary_with_target.tsv

python dynamic_search.py init --config configs/backward_search.example.json --run-dir runs/example
python dynamic_search.py step --run-dir runs/example

# When status is waiting_for_structural, copy/run the repository on Biowulf:
bash biowulf/submit_round.sh runs/example/round_000

# After the jobs finish:
python dynamic_search.py collect --run-dir runs/example
python dynamic_search.py step --run-dir runs/example
python dynamic_search.py status --run-dir runs/example

# Development alternative: each step advances one decision, not one round.
python dynamic_search.py init --development --config configs/backward_search.example.json --run-dir runs/inspect
python dynamic_search.py step --run-dir runs/inspect
python dynamic_search.py trace --run-dir runs/inspect
```

## Validation and interpretation requirements

- Treat the top-20 and beam-width choices as sensitivity parameters and rerun
  with wider values before making a biological exclusion claim.
- Report AlphaPulldown protocol ID, model settings, template-date cutoff,
  interaction metric, reference score, BF bounds, and weight.
- Inspect low-scoring structural models for sequence length, disorder, membrane
  topology, and MSA limitations. A weak model is evidence under the chosen
  heuristic, not proof of no physical interaction.
- Retain missing coverage as missing; never substitute a failed batch job with a
  negative score.
- Compare recovered paths with held-out known pathways before calibrating or
  interpreting the beam as a sensitivity-optimized biological network.

## Implemented scope and next GUI step

This version implements the complete portable engine, state/cache layer,
AlphaPulldown round packaging, Biowulf batch scripts, result collection, and
tests. It is intentionally CLI-first because Biowulf is batch-oriented. The
stepwise development mode supplies a self-contained HTML inspection view while
keeping state transitions command-driven. The existing Version 1 full-graph
GUI remains intact; a later Version 2 control panel can call the same engine and
present `state.json`/round outputs without changing the scientific algorithm.
