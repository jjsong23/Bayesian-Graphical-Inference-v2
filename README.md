# Bayesian Graphical Inference of Renal Signaling

This research codebase builds and explores a probabilistic signaling graph for renal principal cells, with a current biological focus on vasopressin/PKA regulation of Aqp2 throughout the collecting duct. Bayesian node and edge inference first produces an undirected association graph; an auditable ontology-plus-OmniPath layer can then partially orient supported edges for signal-propagation path searches. It combines heterogeneous evidence in four stages:

1. **Node selection** estimates which candidate signaling participants are relevant to the biological system.
2. **Edge characterization** estimates the probability that two selected nodes are associated.
3. **Path inference** ranks plausible paths from a chosen signaling receptor or regulator to a target protein.
4. **Temporal validation** optionally tests whether measured dDAVP phosphoproteomic response times are consistent with the proposed path order.

The local workbench exposes evidence selection, per-dataset normalization controls, evidence weights, optional scope-aware Bayes factors below 1 for node nondetections and unsupported edge pairs, partial directionality, path constraints, external-target insertion, and optional uncertainty-aware temporal annotations. Completed runs include posterior-distribution plots, a per-node/per-edge evidence ledger that reconstructs the Bayesian update, and an interactive merged network of the highest-ranked paths. Negative evidence can use either a stage-wide fixed absence factor or a per-dataset continuous mode that removes the positive-only floor, treats eligible nondetections as `x=0`, and lets weak observations produce BF below 1. The optional regularized positive-control calibration stage independently fits node and edge weights/Tq scales to supplied known-present controls. The current Version 2 defaults mirror Version 1: all 12 node streams are enabled, eligible unobserved nodes receive the configured negative-evidence treatment, and scaffold-mediated closure is enabled. Unsupported-edge penalties, calibration, and temporal validation remain disabled by default.

For the lab-specific Aqp2 analysis, node selection additionally exposes rat
proteome and mouse RNA abundance for CCD, OMCD, and IMCD as six separate
streams. These six streams are enabled in the current Version 1/Version 2
default profile; their mapping audits are under `results/collecting_duct_node_selection/` and
`data/node_selection/collecting_duct/` in the companion data archive.

## Repository and data archive

GitHub contains the source code, tests, interface, documentation, and configuration. Large datasets and generated results are intentionally excluded from Git history and distributed separately in the complete project archive. Extract that archive so that `data/`, `results/`, and `outputs/` sit beside `code/` and `gui/`.

Read [CODEX_PROJECT_CONTEXT.md](CODEX_PROJECT_CONTEXT.md) first when taking over the project in a new Codex session. It records the current scientific conventions, important entry points, canonical artifacts, and decisions that must not be silently reversed. See [DATA_AND_RESULTS.md](DATA_AND_RESULTS.md) for the data handoff layout.

For the exact changes in the current code state, read [CHANGELOG.md](CHANGELOG.md)
and the latest entry under [`docs/lab_notebook/`](docs/lab_notebook/).

## Architecture at a glance

```text
candidate proteins --node evidence--> selected signaling nodes + messengers
       |                                      |
       |                                      v
       |                         unordered node-pair hypotheses
       |                                      |
       |                              edge evidence integration
       |                                      v
       +--------------------------> symmetric edge probabilities
                                              |
                              ontology + OmniPath traversal constraints
                                              v
                           ranked paths to an internal or external target
                                              |
                         optional temporal audit + interactive merged network
```

Node and edge hypotheses begin at an independent probability of 0.5 by
default. Evidence is integrated in odds space as
`posterior odds = prior odds × product(BF_i ^ weight_i)`. The graph is not a
competition in which probabilities sum to one. Edge existence remains an
undirected hypothesis even when the separate propagation layer disallows one
traversal direction.

## Version 2: backward frontier search

Version 2 now provides a complete GUI-and-CLI pipeline: it performs configurable
node selection itself, materializes the enabled inexpensive edge evidence, and
then launches the bounded backward search. Open `http://127.0.0.1:8765/v2.html`
after running `python launch.py`; no separate Version 1 node-selection run is
required. See [`docs/V2_END_TO_END_PIPELINE.md`](docs/V2_END_TO_END_PIPELINE.md)
for the stage boundaries, equations, GUI controls, and Biowulf checkpoint flow.

The frontier stage avoids structurally scoring the complete graph. The current
end-to-end initialization still materializes the selected nodes' inexpensive
edge-evidence matrix for audit and compatibility with the Version 1 engine;
AlphaPulldown itself is strictly frontier-bounded. The search begins at a selected endpoint,
scores candidate upstream parents with cheap evidence, retains 20 candidates
per current frontier path for AlphaPulldown, and keeps the five best distinct
upstream frontier nodes after structural scores are returned. The search is
resumable and caches raw pair evidence so previously predicted pairs are not
recomputed.

The Biowulf structural protocol now runs exactly one prediction per pair
(`model_1_multimer_v3`, one prediction), performs heavy model work on
node-local `$LSCRATCH`, and returns only compact score/structure/provenance
artifacts. Run-private monomer features are deleted after the frontier moves
away from them, with a byte-level cleanup audit. Existing valid results under
`/data/$USER/ppi_screen` can be imported without copying the model tree, and
HuRI binary interactions are available as optional mapped experimental edge
evidence. See
[`docs/V2_STORAGE_CACHE_HURI_2026-09-15.md`](docs/V2_STORAGE_CACHE_HURI_2026-09-15.md).
When HuRI is enabled, a verified reported HuRI positive replaces the
AlphaPulldown request for that pair by default. The HuRI BF is applied once as
cheap evidence and the structural term remains neutral, so the same interaction
is not double-counted. HuRI non-reporting never suppresses AlphaPulldown.

Scaffold-mediated closure is restricted to physical evidence. For two proteins
to receive the closure BF, each must have either an accepted cached
AlphaPulldown/AlphaFold score strictly above `Anchor >` (default 0.90) or a
reported-positive HuRI interaction with the same exact `adaptor_scaffold`
protein. Localization, STRING, OmniPath, kinase prediction, STITCH, and the
integrated cheap-edge posterior cannot create closure anchors. The default
support likelihood remains 0.90 (`BF = 1.8`), and nonqualifying pairs remain
neutral. Closure is recalculated from the persistent physical-result cache at
each Version 2 frontier so newly collected results can support later rounds.

```powershell
python dynamic_search.py init --config configs/backward_search.example.json --run-dir runs/example
python dynamic_search.py step --run-dir runs/example
python dynamic_search.py status --run-dir runs/example
```

The complete workflow has parallel GUI/CLI controls:

```powershell
python dynamic_search.py pipeline-config --output configs/my_pipeline.json
python dynamic_search.py pipeline-init --config configs/my_pipeline.json --run-dir runs/my_run
python dynamic_search.py pipeline-step --run-dir runs/my_run --until-checkpoint
```

For development and scientific inspection, add `--development` to `init`.
Every subsequent `step` then performs exactly one visible action—candidate
enumeration, cheap-evidence updating, shortlist selection, structural-cache
and batch preparation, final edge updating, or frontier selection—and pauses.
The run writes complete keep/reject ledgers plus a continuously refreshed
`development_trace.html`:

```powershell
python dynamic_search.py init --development --config configs/backward_search.example.json --run-dir runs/inspect
python dynamic_search.py step --run-dir runs/inspect
python dynamic_search.py status --run-dir runs/inspect
python dynamic_search.py trace --run-dir runs/inspect
```

See [`docs/V2_STEPWISE_DEVELOPMENT_MODE.md`](docs/V2_STEPWISE_DEVELOPMENT_MODE.md)
for every pause point, output table, and equation.

The Version 2 run panel also polls a live frontier view. For each active source
node it draws the top 20 cheap-evidence candidates, shows every BF/weight/log-BF
contribution on selection, detects finished ipTM files before collection, and
changes each candidate to retained or culled after combined-evidence pruning.

When a round reports `waiting_for_structural`, the repository can be copied to
NIH Biowulf and submitted using `biowulf/submit_round.sh`. The batch side uses
only deterministic scripts, Slurm, SQLite, and the installed AlphaPulldown
module—no Codex, ChatGPT, LLM API, or AI assistant. See
[`docs/V2_BACKWARD_FRONTIER_PLAN.md`](docs/V2_BACKWARD_FRONTIER_PLAN.md) for the
scientific plan and [`biowulf/README.md`](biowulf/README.md) for exact commands.

## Quick start

Python 3.10 or newer is recommended.

```powershell
python -m pip install -r requirements.txt
python launch.py
```

The workbench binds to `127.0.0.1` and opens at `http://127.0.0.1:8765/`.
Use `python launch.py --foreground` when you prefer a server that stops with
Ctrl+C, or `python launch.py --check` to validate a freshly extracted release.
The PowerShell launcher remains available at `gui/run_workbench.ps1`.

For a complete scientific and software overview, read
[`docs/PROJECT_COMPLETE_SUMMARY_2026-08-25.md`](docs/PROJECT_COMPLETE_SUMMARY_2026-08-25.md).
The reproducible lab-sharing release procedure is documented in
[`docs/SHAREABLE_RELEASE_2026-08-25.md`](docs/SHAREABLE_RELEASE_2026-08-25.md).

## Main directories

- `code/`: Bayesian utilities and reproducible analysis modules.
- `code/backward_search/`: bounded backward beam search and structural-score cache.
- `biowulf/`: portable AlphaPulldown Slurm jobs and submission instructions.
- `configs/`: auditable Version 2 configuration examples.
- `gui/`: minimal local web workbench, evidence registry, and workflow engine.
- `docs/`: detailed methods, target-extension documentation, and lab notebook.
- `notebooks/`: reserved for exploratory notebooks.
- `data/`: source and precomputed evidence tables; supplied in the data archive.
- `results/`: immutable analysis outputs; supplied in the data archive.
- `outputs/`: additional generated artifacts; supplied in the data archive.
- `deliverables/`: ignored release archives produced for lab sharing.

The Tq-aware leave-one-stream-out analysis is implemented in
`code/sensitivity_analysis/analyze_end_to_end_ablation.py`. By default it tests
all registered streams, including optional streams in matched add-one contexts,
while holding non-focal settings fixed and rebuilding downstream nodes, edges,
directionality, and paths as appropriate.

## Tests

From the repository root:

```powershell
python -m unittest discover -s code -p "test_bayes_factors.py"
python -m unittest discover -s code -p "test_engine.py"
python -m unittest discover -s code\path_finding -p "test_*.py"
python -m unittest discover -s code\sensitivity_analysis -p "test_*.py"
python -m unittest discover -s gui -p "test_*.py"
```

## Scientific status

This is an evolving hypothesis-generation workflow, not a clinically validated model. Path scores rank graph-supported hypotheses; they are not calibrated probabilities that an entire biological pathway is correct. Source-specific licenses and redistribution restrictions remain applicable to the data in the companion archive.
