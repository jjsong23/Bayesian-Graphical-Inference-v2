# Bayesian Graphical Inference of Renal Signaling

This research codebase builds and explores a probabilistic signaling graph for renal principal cells, with a current biological focus on vasopressin/PKA regulation of Aqp2 throughout the collecting duct. Bayesian node and edge inference first produces an undirected association graph; an auditable ontology-plus-OmniPath layer can then partially orient supported edges for signal-propagation path searches. It combines heterogeneous evidence in four stages:

1. **Node selection** estimates which candidate signaling participants are relevant to the biological system.
2. **Edge characterization** estimates the probability that two selected nodes are associated.
3. **Path inference** ranks plausible paths from a chosen signaling receptor or regulator to a target protein.
4. **Temporal validation** optionally tests whether measured dDAVP phosphoproteomic response times are consistent with the proposed path order.

The local workbench exposes evidence selection, per-dataset normalization controls, evidence weights, optional scope-aware Bayes factors below 1 for node nondetections and unsupported edge pairs, partial directionality, path constraints, external-target insertion, and optional uncertainty-aware temporal annotations. Completed runs include posterior-distribution plots, a per-node/per-edge evidence ledger that reconstructs the Bayesian update, and an interactive merged network of the highest-ranked paths. Negative evidence can use either a stage-wide fixed absence factor or a per-dataset continuous mode that removes the positive-only floor, treats eligible nondetections as `x=0`, and lets weak observations produce BF below 1. The optional regularized positive-control calibration stage independently fits node and edge weights/Tq scales to supplied known-present controls. Negative evidence, calibration, scaffold closure, and temporal validation are disabled by default.

For the lab-specific Aqp2 analysis, node selection additionally exposes rat
proteome and mouse RNA abundance for CCD, OMCD, and IMCD as six separate
streams. They are off in the generic default profile so the validated 891-node
baseline remains reproducible; the all-collecting-duct profile and its mapping
audits are under `results/collecting_duct_node_selection/` and
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
python -m unittest discover -s code\path_finding -p "test_*.py"
python -m unittest discover -s code\sensitivity_analysis -p "test_*.py"
python -m unittest discover -s gui -p "test_*.py"
```

## Scientific status

This is an evolving hypothesis-generation workflow, not a clinically validated model. Path scores rank graph-supported hypotheses; they are not calibrated probabilities that an entire biological pathway is correct. Source-specific licenses and redistribution restrictions remain applicable to the data in the companion archive.
