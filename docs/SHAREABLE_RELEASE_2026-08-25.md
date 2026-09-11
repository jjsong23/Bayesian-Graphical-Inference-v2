# Shareable release — 2026-08-25

## Deliverables

The release builder creates three local deliverables under `deliverables/`:

1. `BGI_Plots_2026-08-25.zip`
2. `BGI_App_2026-08-25.zip`
3. `SHA256SUMS_2026-08-25.txt`

The directory is ignored by Git because these are generated binary releases.
The source code, methods, tests, and package builder remain version-controlled.

The short archive names and internal `BGI_Plots/` and `BGI_App/` roots are
intentional. The original nested plot layout made
`collecting_duct_audit_overview.png` resolve to exactly 260 characters when the
archive was extracted beside the ZIP in the project's `deliverables` folder,
which triggered Windows Explorer error `0x80010135: Path too long`. The builder
now rejects any layout whose hypothetical extracted path exceeds 240 characters.
At this release, the longest plot path is 188 characters and the longest app
path is 238 characters in that same extraction location.

## Plot archive

The plot archive contains all project-generated PNG figures found at release
time, excluding Git metadata, package-manager dependencies, and the deliverables
directory itself. PNGs use short, numbered filenames inside the archive so
Windows Explorer can extract them from a deep directory without exceeding its
legacy path limit. `PLOT_MANIFEST.tsv` maps every numbered PNG back to its full
project-relative source path and records byte size, dimensions, and SHA-256 hash.

The archive includes:

- all sensitivity-analysis figures;
- node and edge Tq-versus-Bayes-factor curves;
- node and edge Tq-versus-KL-information curves;
- all 20 individual evidence-score distributions with their cutoff/reference;
- the consolidated 12-stream node and 8-stream edge distribution panels;
- posterior, rank-progression, class-composition, and nonbaseline-count plots;
- edge-probability histogram/violin figures;
- stepwise adjacency heatmaps;
- localization and colocalization figures; and
- path, PPI-filtering, and AlphaFold planning graphics present in the project.

Six sensitivity figures that existed only as SVG were rendered to equivalent
PNG files before packaging. The original SVGs remain available in the analysis
directory for vector editing.

## Application archive

The application ZIP has a single top-level
`BGI_App/` directory and is intended to launch after
normal Python dependency installation:

```text
python -m pip install -r requirements.txt
python launch.py --check
python launch.py
```

`python launch.py` validates the release, starts the local workbench, and opens
`http://127.0.0.1:8765/`. `--foreground` keeps the server attached to the
terminal, `--no-browser` suppresses browser opening, and `--check` validates
without starting the server.

The package includes:

- `launch.py`, the GUI, workflow engine, browser client, analysis code, tests,
  and project documentation;
- the 9,170-candidate node factor catalog and the 891-node selected seed;
- all eight declared edge-stream catalogs, including the optional alternatives
  and derived closure definition;
- the raw score tables/matrices needed to vary each supported Tq or odds
  reference in the GUI;
- raw localization, STRING, HPA, OmniPath, STITCH, UniProt, and KinasePredictor
  sources required to characterize newly added nodes;
- the official KinasePredictor output matrices used by the local scorer; and
- the phosphoproteomic workbook used for optional temporal path validation.

The package deliberately excludes:

- the 3.35 GB working SQLite incremental-edge cache;
- historical `results/gui_runs/` outputs;
- previous generated analyses and duplicate exploratory workbooks;
- Git history, Python environments, `node_modules`, logs, and PID files.

The omitted SQLite database is a disposable performance cache. The extracted
application creates a fresh database and computes only previously unseen pairs
involving added nodes. The immutable 891-node seed factors remain available in
the packaged precomputed catalogs.

## Verification performed

The release process performed all of the following:

- opened every PNG with Pillow and recorded its dimensions;
- CRC-tested both ZIP archives with Python's `zipfile.testzip()`;
- extracted the application ZIP into a clean verification directory;
- extracted both archives with Windows PowerShell `Expand-Archive` into folders
  named after their ZIPs inside the deep project `deliverables` directory;
- verified that the previously failing collecting-duct PNG exists and matches
  the SHA-256 value in `PLOT_MANIFEST.tsv`;
- ran `python launch.py --check` successfully from the extracted copy;
- started the extracted application on an alternate local port;
- received `{"status": "ok"}` from `/api/health`; and
- received HTTP 200 for the main page.

The clean-install report contained 12 node evidence streams, 8 edge evidence
streams, 9,170 node candidates, 4 default-enabled node streams, and 6
default-enabled edge streams. Optional packages unavailable in the verification
runtime were correctly reported as optional; installing `requirements.txt`
enables all current features, including SciPy parameter calibration.

## Rebuild command

From the project root:

```text
python code/release/build_shareable_packages.py
```

The builder checks that every required runtime input exists, writes package
manifests, atomically replaces the two named release archives, verifies both
archives, and rewrites the companion checksum file.

## Redistribution warning

This technical package is suitable for controlled lab sharing. Before public
distribution, a scientist or data steward must review the licenses and citation
requirements of every bundled external dataset and the official KinasePredictor
files. A functioning local application bundle is not automatically a legally
approved public data release.
