# Version 2 compact AlphaPulldown, legacy-screen reuse, and HuRI evidence

## Purpose

This revision reduces Biowulf storage use while preserving an auditable
structural-evidence record. It also reuses validated results from the earlier
`ppi_screen` and adds the published Human Reference Interactome (HuRI) as an
optional inexpensive edge-evidence source.

## Exactly one AlphaFold prediction per pair

The V2 AlphaPulldown job now fixes both controls that determine prediction
count:

- `num_predictions_per_model = 1`; and
- `model_names = model_1_multimer_v3`.

This means one model seed from one AlphaFold-Multimer model: **one total
prediction per unordered protein pair**. Setting only
`num_predictions_per_model=1` is insufficient because AlphaPulldown otherwise
runs its default list of five multimer models.

Each GPU task performs modeling and `run_get_good_pae.sh` analysis entirely in
the Slurm job's node-local `$LSCRATCH`. A task is valid only if it produces a
nonempty `predictions_with_good_interpae.csv` containing at least one finite
ipTM score in `[0,1]`. Only these compact artifacts return to `/data`:

- `predictions_with_good_interpae.csv`;
- `ranked_0.pdb`, when emitted;
- `ranking_debug.json` and `timings.json`, when emitted; and
- `prediction_manifest.tsv`, which records the pair, model, recycle count,
  prediction count, and Slurm job.

Large result pickles, PAE JSON files, and intermediate model products remain on
node-local scratch and disappear when the task ends. Failed tasks exit nonzero;
the collector does not impute a score.

## MSA/feature lifetime

AlphaPulldown monomer feature objects are written to a **run-private**
`backward_search/features/` directory. This isolation is essential: one run's
cleanup can never invalidate another run's active jobs or cache.

Cleanup occurs only after all scores for a round have been collected,
integrated, and the next frontier has been selected. The new frontier proteins'
feature objects are retained because they are needed in the next round. Feature
objects and any explicit MSA directories for the old current node and rejected
candidates are deleted. At search completion, all features generated for the
completed round are eligible for deletion. Source FASTA, compact pair results,
raw ipTM scores, Bayesian ledgers, and the SQLite pair cache are retained.

Every deletion is recorded in `feature_cleanup.tsv` with symbol, exact path,
byte count, completed round, timestamp, and reason. This deliberately trades
possible future MSA recomputation for substantially lower persistent storage.
No cleanup runs while the Slurm array is active.

## Reusing `/data/$USER/ppi_screen`

The legacy importer scans `ppi_screen/scores/*.csv`. A row is accepted only
when:

1. its name can be parsed as `<symbol>_and_<symbol>`;
2. ipTM is finite and lies in `[0,1]`; and
3. `ppi_screen/hits/<pair>/ranked_0.pdb` exists and is nonempty.

Canonical unordered duplicates are resolved by ranking confidence and then
ipTM. Importantly, **low ipTM scores are retained**. The importer does not keep
only favorable predictions, because that would erase structural evidence
against an interaction.

Accepted scores enter `runtime/pair_cache.sqlite3` under
`legacy_ppi_screen_5models_predictions1`. They are not relabeled as outputs of
the new one-model protocol. The active search checks its exact protocol first
and then only the explicitly configured compatible protocol IDs. The score
ledger records which protocol and source file supplied each reused value.

Run the importer directly:

```bash
python dynamic_search.py import-legacy-screen \
  --screen-dir /data/$USER/ppi_screen \
  --cache runtime/pair_cache.sqlite3 \
  --audit-dir runtime/legacy_ppi_screen_import
```

Or submit its small Slurm job:

```bash
bash biowulf/submit_ppi_screen_import.sh
```

The end-to-end pipeline also auto-imports the sibling `../ppi_screen` when it
exists. Only scores are transferred into SQLite; the large model and MSA trees
are not duplicated. The row audit and JSON summary explain every rejection.

## HuRI evidence and the meaning of an absent pair

HuRI (Luck et al., Nature 2020, DOI `10.1038/s41586-020-2188-x`) is a map of
experimentally detected human binary protein interactions. Its downloadable
TSV uses human Ensembl gene identifiers. V2 maps those identifiers to mouse
nodes using the existing audited mouse-human HPA orthology table.

A reported mapped pair receives the configurable `positive_bayes_factor`
(default 5.0) and stream weight (default 1.0). These defaults are transparent
heuristics, not an experimentally estimated likelihood ratio, and should be
calibrated or varied in sensitivity analyses.

HuRI does **not** provide a complete public pair-level table saying that every
unreported pair was tested successfully and found not to interact. It is also
an incomplete yeast-two-hybrid map. Therefore absence is neutral by default.
The GUI offers two optional weak-negative assumptions:

- `huri_interactor_set`: an unreported pair receives the configurable weak
  negative BF (default 0.9) only if both orthologs occur somewhere in HuRI's
  positive network; or
- `screened_gene_file`: the weak negative applies only if both orthologs occur
  in a user-supplied Ensembl screen/testability list.

Unmapped proteins, small molecules, and pairs outside the chosen scope always
remain neutral. A sparse positive-factor table, scope-node table, mapping
counts, chosen assumption, and caveat are written under the run's
`huri_evidence/` directory. The stepwise cheap-evidence ledger marks each row
as a reported record, a scoped non-report, or an out-of-scope neutral absence.

Download the published table on a login node with:

```bash
bash biowulf/fetch_huri.sh
```

Then enable HuRI in the V2 GUI. The expected default path is
`runtime/reference/huri/HuRI.tsv`.

## Git-based Biowulf updates

Git tracks code, configuration, tests, and documentation. `data/`, `results/`,
`runs/`, `runtime/`, virtual environments, model outputs, and caches remain
ignored. Consequently a pull updates code without uploading laboratory data or
overwriting Biowulf results.

After the repository has an `origin`, update an existing Biowulf checkout with:

```bash
cd /data/$USER/graphical_bayesian_inference_v2
bash biowulf/update_from_github.sh
```

The helper refuses to pull over tracked local edits and uses a fast-forward-only
merge. Existing ignored data and runtime files remain in place.

## Primary references

- AlphaPulldown manual: <https://github.com/KosinskiLab/AlphaPulldown/blob/main/manuals/example_1.md>
- NIH AlphaPulldown guidance: <https://hpc.nih.gov/apps/alphapulldown.html>
- HuRI paper: <https://www.nature.com/articles/s41586-020-2188-x>
- HuRI download portal: <https://www.interactome-atlas.org/download>
