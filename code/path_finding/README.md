# External-target graph extension

`build_target_adjacency_vector.py` adds one mouse protein that is not already
in the active node universe to the current undirected Bayesian graph. It
calculates one target-to-node probability for every current node, preserves the
existing square matrix exactly, and optionally appends the target as its last
row and column.

## Command-line use

From the project root:

```powershell
python code/path_finding/build_target_adjacency_vector.py Aqp2
```

The target may be a mouse gene symbol or UniProt accession. By default, files
are written to `results/path_finding/target_extensions/<symbol>/`.

Useful options:

```text
--output-dir PATH       Override the output directory
--top-n N               Retain N distinct eligible kinases per target site
--no-extended-matrix    Write the vector without the enlarged square matrix
```

## Python use

```python
from pathlib import Path
import sys

project = Path("graphical_bayesian_inference").resolve()
sys.path.insert(0, str(project / "code" / "path_finding"))
from build_target_adjacency_vector import build_target_adjacency_vector

result = build_target_adjacency_vector(
    "Aqp2",
    project_root=project,
)

vector = result.vector
matrix = result.extended_matrix
```

The explicit path insertion avoids a collision between this project's `code`
directory and Python's standard-library `code` module.

## Behavior

- The five streams are mpkCCD localization, observed/curated-site
  KinasePredictor, STRING v12, HPA v25.1 primary localization, and OmniPath
  core.
- Every target-node hypothesis begins at 0.5.
- Missing coverage in a stream is exactly neutral (`BF = 1`).
- Evidence streams are applied sequentially as odds multipliers.
- The graph is unsigned and undirected. OmniPath direction and sign remain in
  audit columns only.
- A target already in the active universe is rejected because its adjacency
  vector already exists in the current matrix.

See `docs/path_finding_target_extension.md` for the complete method and the
Aqp2 validation result.

## Ranked path finding

`find_ranked_paths.py` searches the current graph between a signaling start
node and either an existing or external target. An external target is first
added with the function above; a validated current extension is reused when
available.

```powershell
python code/path_finding/find_ranked_paths.py "PKA subunit a" Aqp2
```

The PKA alias in this example is explicitly mapped to mouse `Prkaca`. Gene
symbols are preferred for all other inputs. Defaults are 25 paths, at most six
hops, edges strictly above the neutral 0.5 baseline, and mechanistic
signal-relay intermediates only.

Eligible internal nodes must have at least one of `receptor`,
`receptor_regulator`, `kinase`, `phosphatase`, `phosphodiesterase`, `gtpase`,
`gtpase_regulator`, `signaling_regulation`, or `second_messenger`. Broad labels
such as `adaptor_scaffold`, `ligand`, `kinase_phosphatase_binding`, and generic
`signaling_process` do not qualify by themselves. The selected start and target
are exempt because they are endpoints rather than intermediates.

Controls include:

```text
--top-k N
--max-hops N
--minimum-edge-probability P
--rebuild-target-extension
--no-auto-extend
--allow-all-intermediates
--exclude-any-scaffold
```

Paths are loopless and undirected. They are ranked by descending product of
edge probabilities, implemented as a hop-limited Yen search with Dijkstra
subproblems and additive `-log(p)` edge costs. Results are stored under
`results/path_finding/ranked_paths/<start>_to_<target>_signal_relay/`.

`--exclude-any-scaffold` is stricter than the default role-based rule: it also
removes multi-role proteins that have a genuine relay class (for example,
kinase) but additionally carry an `adaptor_scaffold` label.
