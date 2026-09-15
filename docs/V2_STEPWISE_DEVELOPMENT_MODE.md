# Version 2 stepwise development mode

## Purpose

Development mode runs the same bounded backward-frontier inference as the
ordinary Version 2 command-line workflow, but one scientifically meaningful
decision at a time. It is intended for algorithm development, debugging, and
methods validation. It does not change the prior, Bayes factors, weights,
structural-score transformation, path score, pruning rules, or stopping rules.

The complete audit is written as machine-readable TSV/JSON files. A
self-contained `development_trace.html` summarizes the current frontier,
chronological decisions, row counts, and previews of all completed stages. The
HTML is a viewing aid; the TSV ledgers and frozen `state.json` are authoritative.

## Start and advance a traced run

From the repository root:

```powershell
python dynamic_search.py init `
  --development `
  --config configs/backward_search.example.json `
  --run-dir runs/aqp2_development

python dynamic_search.py status --run-dir runs/aqp2_development
python dynamic_search.py step --run-dir runs/aqp2_development
```

In development mode, **each `step` command performs exactly one stage**. Run
`status` after any step to see the next stage. The same mode can be enabled in
an auditable configuration file with:

```json
"development_mode": true
```

To rebuild the HTML view at any point:

```powershell
python dynamic_search.py trace --run-dir runs/aqp2_development
```

## One search round, stage by stage

### 1. Enumerate candidates

The current node in each partial backward path is paired with every symbol in
the node table. No edge probability is calculated yet. The engine records
whether each proposed parent is admissible and, when it is not, one of these
explicit reasons:

- ontology or configured eligibility filter;
- self-edge;
- reuse of a node already in that partial path, which would make a cycle; or
- disallowed forward propagation from the proposed parent to the current node.

Output: `round_NNN/candidate_enumeration.tsv`.

### 2. Apply cheap evidence

Every admissible candidate begins at the configured edge prior. Cheap evidence
streams are applied in their configuration order:

```text
log posterior odds = log prior odds + sum(weight_s * log(BF_s))
```

For every candidate-stream update, the ledger stores whether a pair record was
found, the BF actually used, the configured missing BF, the weight, the
weighted log-BF contribution, and the probability before/after the update.
Consequently, a missing BF below 1 is visibly negative evidence, a BF of 1 is
neutral, and a BF above 1 is positive evidence.

Outputs:

- `round_NNN/cheap_evidence_ledger.tsv`: long-form arithmetic audit;
- `round_NNN/cheap_scored_all.tsv`: one row per admissible candidate and its
  cheap-stage rank within that frontier path.

### 3. Select the cheap shortlist

The engine retains the configured `cheap_top_n_per_frontier` candidates from
each current frontier path. If `always_include_receptor` is enabled, the
receptor is appended when it did not make the numerical top N. All rejected
rows remain in the audit with the reason `outside_top_n`; forced receptor rows
are explicitly labeled.

Outputs:

- `round_NNN/cheap_selection.tsv`: every scored candidate plus the decision;
- `round_NNN/cheap_candidates.tsv`: the retained shortlist.

### 4. Check the structural cache and prepare the batch

For each shortlisted protein-protein pair, the engine looks for a raw score
with the exact configured protocol ID and metric. Cache hits are reused.
Structurally inapplicable pairs require no prediction. Only unique cache misses
are exported to the AlphaPulldown package.

Outputs:

- `round_NNN/structural_cache_lookup.tsv`;
- `round_NNN/pairs.tsv`, `sequences.fasta`, and `round_config.json` when there
  are cache misses.

When misses exist, status becomes `waiting_for_structural`. Copy the repository
to Biowulf and follow `biowulf/README.md`. Importing results is itself recorded,
but it does not silently integrate them:

```powershell
python dynamic_search.py import `
  --run-dir runs/aqp2_development `
  --scores path/to/structural_scores.tsv
```

Partial imports are allowed. The run remains waiting until every required pair
has a numeric result. Failed or missing predictions are never imputed as weak
binding.

### 5. Integrate all edge evidence

The imported raw interaction score is converted to the configured bounded,
heuristic structural BF. This structural term is added to the same edge log
odds already updated by cheap evidence. The edge posterior is then multiplied
into the probability of the partial path (implemented in log space).

Outputs:

- `round_NNN/edge_evidence_ledger.tsv`: every cheap and structural contribution
  in exact update order;
- `round_NNN/scored_candidates.tsv`: final edge posterior, partial-path score,
  and pre-beam global rank.

### 6. Select the next frontier

All extensions from all current frontier paths compete globally by partial-path
score. The engine keeps at most `beam_width` extensions and requires their new
frontier nodes to be distinct. Every rejected row is labeled either as a
duplicate proposed frontier node or as outside the global beam width. The
engine then checks receptor and maximum-depth stopping rules. If the search
continues, the retained paths become the next round's frontier.

Outputs:

- `round_NNN/frontier_selection.tsv`;
- updated `round_NNN/scored_candidates.tsv` with the retention flag;
- updated `state.json` and `solutions.json`.

## Run-level audit files

- `state.json`: atomic resumable state and frozen configuration snapshot;
- `development_events.jsonl`: append-only chronological stage history;
- `development_trace.html`: self-contained human-readable trace;
- `solutions.json`: completed receptor-to-target paths.

The state file also records `development_stage`, which is the next action,
`development_step_count`, and the most recent event. Re-running `status` or
`trace` does not advance or alter the inference.

## Interpretation boundary

Stepwise mode explains exactly why the implemented algorithm retained a
candidate. It does not make a beam search exhaustive: a biologically correct
path can still be lost at the cheap top-N or global beam pruning stage. Use the
retained/rejected ledgers to test wider settings and to distinguish “not
examined further” from “evidence against the interaction.” The structural BF
remains a configurable heuristic based on an AlphaPulldown score, not a
calibrated probability of physical binding.
