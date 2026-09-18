# Running Version 2 on NIH Biowulf

These files contain ordinary deterministic Python, Bash, SQLite, and Slurm
operations. They do not invoke Codex, ChatGPT, an API, or another AI assistant.
AlphaPulldown itself performs the explicitly requested structural prediction.

## 1. Transfer a portable copy

Copy the complete repository, including real `data/` and `results/`
directories, to a Biowulf-accessible project directory.  The Windows
development copy may show these two directories as links to Version 1; those
links themselves are not portable to Linux.  Replace them with the actual
directories before transfer.  NIH recommends Globus for transfers larger than
about 10 GB.  `scp`, `sftp`, and `rsync` transfers must target `helix.nih.gov`,
not `biowulf.nih.gov`.

After transfer, log into Biowulf and work from `/data`, not a transient local
Windows path:

```bash
ssh USER@biowulf.nih.gov
cd /data/$USER/graphical_bayesian_inference_v2
```

After the GitHub remote is configured, do not recopy the project for ordinary
code updates. From the existing Biowulf checkout run:

```bash
cd /data/$USER/graphical_bayesian_inference_v2
bash biowulf/update_from_github.sh
```

This performs `git pull --ff-only origin main`, refuses to overwrite local
tracked edits, and leaves ignored `data/`, `results/`, `runs/`, `runtime/`, and
`.venv/` content untouched.

## 2. Create the Python environment once

Do this in an interactive allocation (or an Open OnDemand graphical session),
not while consuming the shared login node:

```bash
sinteractive --cpus-per-task=4 --mem=16g --time=2:00:00
module load python
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
chmod u+x biowulf/*.sh biowulf/*.sbatch
python launch.py --check
exit
```

`launch.py --check` must report `"status": "ok"`.  A missing seed-universe or
factor-table error usually means the transferred repository still contains a
Windows link instead of the actual companion `data/` or `results/` directory.

## 3. Run the complete pipeline from the GUI

The cleanest Windows workflow uses Biowulf's built-in tunnel support.  Start
the login session inside `tmux` so that an accidental terminal disconnect does
not destroy the controller:

```bash
module load tmux
tmux new -s gbi-v2
sinteractive --cpus-per-task=8 --mem=32g --time=1-12 --tunnel
cd /data/$USER/graphical_bayesian_inference_v2
source .venv/bin/activate
python launch.py --check

# PORT1 is assigned by sinteractive; do not substitute a fixed 8766 here.
python launch.py --foreground --no-browser --host 127.0.0.1 --port "$PORT1"
```

When `sinteractive` starts, it prints an exact workstation command such as:

```text
ssh -L 45000:localhost:45000 USER@biowulf.nih.gov
```

Run the printed command in a second **Windows PowerShell** window and leave it
open.  Then visit the matching address in the Windows browser:

```text
http://localhost:45000/v2.html
```

Use the actual port printed for your job.  In the Version 2 page:

1. Set node evidence streams, weights, Tq/reference multipliers, priors, and
   cutoffs.
2. Set the endpoint, receptor, cheap edge evidence, direction rules, frontier
   widths, and AlphaPulldown transformation.
3. Select **Initialize complete run**.
4. Select **Run to checkpoint**.  Node selection, graph preparation, cheap
   evidence integration, directionality, and backward search run automatically.
   In development mode, the run panel draws the top 20 candidates connected to
   each frontier source. Select a candidate to inspect its exact evidence
   factors and interpretation.
5. When the state says that structural evidence is required, select **Submit
   AlphaPulldown**.  The controller submits the prepared CPU-feature and GPU
   Slurm arrays; it does not run those predictions inside the web server.
6. Check the queue with `squeue -u $USER`.  After the arrays finish, select
   **Collect + continue**.  Repeat submit/collect only when the next frontier
   creates genuinely new pairs.  Cached unordered pairs are not recomputed.
   The live graph detects compact ipTM files as tasks finish and previews which
   candidates will strengthen or weaken before formal collection; final
   retained/culled labels appear after collection and frontier selection.

The GUI controller may be stopped while Slurm jobs run.  The run state is
written atomically under `runs/`; restart the server and reopen the saved run
to continue.  Do not delete `runtime/pair_cache.sqlite3`, because that is the
cross-round structural-score cache.

An NIH Open OnDemand Graphical Session is an alternative to the PowerShell
tunnel: open its terminal, activate the environment, run
`python launch.py --foreground --no-browser --host 127.0.0.1 --port 8766`, and
open `http://127.0.0.1:8766/v2.html` in the session's Firefox browser.  The
explicit `sinteractive --tunnel` route is preferable when you want more than
the Graphical Session's fixed resources.

## 4. Equivalent command-line controller

The web page is a front end to the following complete-pipeline operations; it
is not necessary to submit them manually when using the GUI:

```bash
python dynamic_search.py pipeline-config --output configs/my_pipeline.json
python dynamic_search.py pipeline-init --config configs/my_pipeline.json --run-dir runs/aqp2_complete
python dynamic_search.py pipeline-step --run-dir runs/aqp2_complete --until-checkpoint
python dynamic_search.py pipeline-submit --run-dir runs/aqp2_complete
# After Slurm completes:
python dynamic_search.py pipeline-collect --run-dir runs/aqp2_complete
python dynamic_search.py pipeline-step --run-dir runs/aqp2_complete --until-checkpoint
python dynamic_search.py pipeline-status --run-dir runs/aqp2_complete
python dynamic_search.py pipeline-trace --run-dir runs/aqp2_complete
```

The older low-level commands below remain useful for debugging a frontier that
has already been initialized, but are not the recommended full-pipeline entry
point.

```bash

# Run once if data/backward_search/node_sequences.fasta has not been prepared.
python dynamic_search.py build-fasta \
  --nodes data/node_selection/node_universe_combined_nonzero.tsv \
  --uniprot data/edge_characterization/kinase_predictor/phosphosite_database/raw/uniprot_mouse_reference_proteome.tsv.gz \
  --output data/backward_search/node_sequences.fasta

# If the endpoint is outside the signaling universe, append its HPA vector to
# the internal HPA pair catalog without constructing a square matrix.
python dynamic_search.py build-target-evidence \
  --project-root . \
  --target Aqp2 \
  --base-evidence results/backend_bayes_factor_catalogs/edge_factors_891/hpa_primary_bf_gt1.tsv.gz \
  --output data/backward_search/hpa_primary_with_target.tsv

python dynamic_search.py init \
  --config configs/backward_search.example.json \
  --run-dir runs/aqp2_backward

python dynamic_search.py step --run-dir runs/aqp2_backward
bash biowulf/submit_round.sh runs/aqp2_backward/round_000
```

After the prediction array finishes:

```bash
source .venv/bin/activate
python dynamic_search.py collect --run-dir runs/aqp2_backward
python dynamic_search.py step --run-dir runs/aqp2_backward
```

Repeat `submit_round.sh`, `collect`, and `step` for each newly prepared round.
Use `python dynamic_search.py status --run-dir ...` at any time. The run is
resumable because `state.json` is written atomically and raw structural scores
are cached by unordered pair, protocol ID, and metric in
`runtime/pair_cache.sqlite3`.

Operational notes:

- AlphaPulldown is loaded from the NIH module; no installation is attempted by
  the job scripts.
- Feature generation is capped at ten concurrent tasks, following Biowulf's
  warning about alignment pressure on shared filesystems.
- Each prediction array task requests one A100 GPU. Set
  `GBI_MAX_CONCURRENT_GPU_JOBS` before submission to change concurrency.
- Set `GBI_ALPHAFOLD_DATA_DIR` if the installed module requires an explicit
  database path.
- Set `GBI_DB_PRESET` to `full_dbs` (the default) or `reduced_dbs`. The selected
  preset is passed to AlphaPulldown, written to `feature_protocol.tsv`, and
  repeated in every compact pair manifest so the feature protocol remains
  auditable.
- Set `GBI_MAX_TEMPLATE_DATE` to the scientifically intended cutoff. The
  default is the Version 2 implementation date (`2026-09-11`) so reruns do not
  silently use later templates.
- `run_get_good_pae.sh` creates the tabular ipTM output consumed by `collect`.
  Failed analyses are listed in `collection_summary.json` and are not imputed.
- FASTA identifiers are restricted to letters, digits, period, underscore, and
  hyphen so they can be passed safely through the batch scripts.

### Recovering from an HH-suite full-database feature failure

Some proteins can make HH-suite's full-database search fail with a
`MergeMasterSlave` match-state error. Do not exclude the protein and do not let
the dependent GPU array proceed with an incomplete feature set. Cancel both
arrays, move the run-private feature directory to `/scratch`, and regenerate
the entire round with one consistent reduced-database protocol:

```bash
RUN=runs/RUN_ID
ROUND="$RUN/backward_search/round_000"
OLD_FEATURES="$RUN/backward_search/features"
BACKUP="/scratch/$USER/gbi_v2_feature_backups/RUN_ID_full_dbs_$(date +%Y%m%d_%H%M%S)"

scancel FEATURE_JOB_ID PREDICTION_JOB_ID
mkdir -p "$(dirname "$BACKUP")"
mv -- "$OLD_FEATURES" "$BACKUP"

GBI_DB_PRESET=reduced_dbs \
GBI_MAX_CONCURRENT_GPU_JOBS=4 \
bash biowulf/submit_round.sh "$ROUND"
```

Replace all placeholders with the run and job IDs printed by the GUI. Moving
the old directory is deliberate: `--skip_existing=True` must not silently mix
full- and reduced-database feature objects. The old features remain recoverable
in `/scratch` while the replacement jobs run.

## 5. Compact structural storage policy

Each unordered pair now receives exactly one prediction: one prediction from
`model_1_multimer_v3`. Both settings are frozen in each round's
`round_config.json`; `num_predictions_per_model=1` alone would still invoke all
five default AlphaFold-Multimer models.

GPU model generation and analysis occur under the job's `$LSCRATCH`. Only a
validated `predictions_with_good_interpae.csv`, `ranked_0.pdb` when available,
small ranking/timing metadata, and `prediction_manifest.tsv` return to `/data`.
Large result pickles and PAE JSON products are not copied back.

Features are stored at `runs/<run>/backward_search/features/`, never in a
cross-run shared folder. After a completed round is integrated and the next
frontier is selected, the engine deletes features for proteins it moved away
from and rejected candidates. It retains only new-frontier features; at search
completion it retains none. `feature_cleanup.tsv` records all removals. Do not
manually remove features while `gbi_features` or `gbi_pairs` jobs are running.

## 6. Import the existing ppi_screen without copying it

The end-to-end pipeline automatically detects a sibling
`/data/$USER/ppi_screen` by default. To validate and import it before starting a
GUI run, submit:

```bash
cd /data/$USER/graphical_bayesian_inference_v2
bash biowulf/submit_ppi_screen_import.sh
squeue -u "$USER"
```

This copies no MSA or model directory. It validates each finite `[0,1]` ipTM
against a nonempty `hits/<pair>/ranked_0.pdb`, retains low scores as negative
structural evidence, and writes only compact records into
`runtime/pair_cache.sqlite3`. Review
`runtime/legacy_ppi_screen_import/legacy_ppi_screen_import_summary.json` and the
compressed row audit afterward. Legacy results are tagged as a five-model
protocol and reused only through the explicit compatibility setting.

## 7. Optional HuRI binary interaction evidence

On a login node, download the published HuRI pair table:

```bash
cd /data/$USER/graphical_bayesian_inference_v2
bash biowulf/fetch_huri.sh
```

Then enable HuRI in the V2 GUI. Reported mapped pairs receive the chosen
positive BF. Unreported pairs remain neutral by default because the public
network is not a complete pair-level failed-test table. The GUI can apply a
weak negative BF only within an explicit node scope; this assumption is written
to every run's HuRI summary and evidence ledger.

Leave **Verified HuRI positives replace AlphaPulldown** enabled to use a
reported HuRI hit as the structural answer for that pair. The HuRI BF is counted
once and the pair is omitted from the AlphaPulldown arrays. Unreported or weak-
negative HuRI pairs still require AlphaPulldown when structurally applicable.

See `docs/V2_STORAGE_CACHE_HURI_2026-09-15.md` for the complete scientific and
storage rationale.

Primary references:

- NIH AlphaPulldown instructions: https://helixweb.nih.gov/apps/alphapulldown.html
- NIH Biowulf user guide: https://hpc.nih.gov/docs/userguide.html
- AlphaPulldown source/manuals: https://github.com/KosinskiLab/AlphaPulldown
