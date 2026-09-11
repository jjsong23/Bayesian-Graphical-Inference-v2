#!/usr/bin/env python3
"""Build undirected secondary-messenger/protein edge evidence from STITCH v5.

STRING itself contains protein-protein associations only. Its sister database,
STITCH, uses compatible confidence scores for protein-chemical associations.
This script maps the project's curated second messengers to PubChem CIDs,
filters the mouse STITCH network, maps STITCH protein IDs through the audited
STRING mapping, and integrates each unique messenger-protein pair once.

STITCH reports scores from 0.150 to 0.999. The reporting floor is used as a
conservative reference probability. For score s and reference r=0.150:

    BF = max(1, odds(s) / odds(r))

Missing records are neutral (BF=1). CIDm/CIDs representations and multiple
PubChem records for one curated messenger are collapsed by maximum confidence,
so no pair is counted twice. Direction and effect are intentionally not modeled.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DEFAULT = Path(__file__).resolve().parents[2]
STITCH_VERSION = "5.0"
STITCH_REPORTING_FLOOR = 0.150
EPSILON = 1e-12
SOURCE_URL = (
    "http://stitch.embl.de/download/protein_chemical.links.v5.0/"
    "10090.protein_chemical.links.v5.0.tsv.gz"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_DEFAULT,
        help=f"Project directory (default: {PROJECT_DEFAULT})",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Rows per STITCH input chunk.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_logit(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(np.asarray(values, dtype=float), 1e-12, 1 - 1e-12)
    return np.log(clipped) - np.log1p(-clipped)


def stable_expit(values: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(values, -700.0, 700.0)))


def stitch_bayes_factor(
    score: np.ndarray | pd.Series,
    reference: float = STITCH_REPORTING_FLOOR,
) -> np.ndarray:
    if not 0 < reference < 1:
        raise ValueError("STITCH reference score must be between zero and one")
    probability = np.asarray(score, dtype=float)
    if (~np.isfinite(probability) | (probability <= 0) | (probability >= 1)).any():
        raise ValueError("STITCH scores must be finite probabilities between zero and one")
    log_bf = stable_logit(probability) - math.log(reference / (1.0 - reference))
    return np.maximum(1.0, np.exp(np.clip(log_bf, -700.0, 700.0)))


def joined(values: pd.Series) -> str:
    return ";".join(sorted({str(value) for value in values if str(value)}))


def load_messenger_mapping(path: Path) -> tuple[pd.DataFrame, dict[int, str]]:
    mapping = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {
        "messenger_node",
        "display_symbol",
        "preferred_name",
        "pubchem_cids",
        "mapping_method",
        "mapping_note",
    }
    if set(mapping.columns) != required:
        raise ValueError(f"Unexpected messenger mapping columns: {mapping.columns.tolist()}")
    if mapping["messenger_node"].duplicated().any():
        raise ValueError("Messenger mapping contains duplicate nodes")
    cid_to_messenger: dict[int, str] = {}
    for row in mapping.to_dict(orient="records"):
        for text in str(row["pubchem_cids"]).split(";"):
            cid = int(text)
            previous = cid_to_messenger.get(cid)
            if previous is not None and previous != row["messenger_node"]:
                raise ValueError(f"PubChem CID {cid} maps to multiple messenger nodes")
            cid_to_messenger[cid] = str(row["messenger_node"])
    return mapping, cid_to_messenger


def filter_stitch(
    path: Path,
    cid_to_messenger: dict[int, str],
    chunk_size: int,
) -> tuple[pd.DataFrame, int]:
    if chunk_size <= 0:
        raise ValueError("chunk size must be positive")
    retained: list[pd.DataFrame] = []
    raw_rows = 0
    expected = ["chemical", "protein", "combined_score"]
    for chunk in pd.read_csv(
        path,
        sep="\t",
        dtype={"chemical": str, "protein": str, "combined_score": int},
        chunksize=chunk_size,
    ):
        if chunk.columns.tolist() != expected:
            raise ValueError(f"Unexpected STITCH columns: {chunk.columns.tolist()}")
        raw_rows += len(chunk)
        pubchem = pd.to_numeric(chunk["chemical"].str[4:], errors="coerce")
        keep = pubchem.isin(cid_to_messenger)
        if not keep.any():
            continue
        selected = chunk.loc[keep].copy()
        selected["pubchem_cid"] = pubchem.loc[keep].astype(np.int64)
        selected["messenger_node"] = selected["pubchem_cid"].map(cid_to_messenger)
        selected["chemical_representation"] = selected["chemical"].str[:4]
        retained.append(selected)
    if not retained:
        return pd.DataFrame(), raw_rows
    result = pd.concat(retained, ignore_index=True)
    if not result["combined_score"].between(150, 999).all():
        raise ValueError("STITCH score outside the documented 150-999 range")
    return result, raw_rows


def collapse_pairs(filtered: pd.DataFrame) -> pd.DataFrame:
    if filtered.empty:
        return filtered
    grouped = filtered.groupby(
        ["messenger_node", "protein"], sort=False, as_index=False
    ).agg(
        pubchem_cids=("pubchem_cid", lambda values: joined(values.astype(str))),
        stitch_chemical_ids=("chemical", joined),
        chemical_representations=("chemical_representation", joined),
        supporting_stitch_rows=("chemical", "size"),
        stitch_combined_score=("combined_score", "max"),
    )
    grouped["stitch_score"] = grouped["stitch_combined_score"] / 1000.0
    return grouped


def write_matrix(path: Path, symbols: list[str], values: np.ndarray) -> None:
    frame = pd.DataFrame(values, index=symbols, columns=symbols)
    frame.index.name = "symbol"
    frame.to_csv(path, sep="\t", float_format="%.9g")


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    universe_path = project / "data/node_selection/node_universe_combined_nonzero.tsv"
    mapping_path = (
        project
        / "data/edge_characterization/stitch/v5.0/processed/"
        "secondary_messenger_pubchem_mapping.tsv"
    )
    raw_path = (
        project
        / "data/edge_characterization/stitch/v5.0/raw/"
        "10090.protein_chemical.links.v5.0.tsv.gz"
    )
    string_mapping_path = (
        project
        / "data/edge_characterization/string/v12.0/processed/"
        "node_to_string_mapping.tsv"
    )
    current_dir = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa_omnipath"
    )
    current_matrix_path = current_dir / "combined_adjacency_matrix.tsv"
    processed_dir = mapping_path.parent
    output_dir = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa_omnipath_stitch"
    )
    backend_dir = project / "results/backend_bayes_factor_catalogs/edge_factors_891"
    for path in (
        universe_path,
        mapping_path,
        raw_path,
        string_mapping_path,
        current_matrix_path,
    ):
        if not path.exists():
            raise FileNotFoundError(f"Required input is missing: {path}")
    processed_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    backend_dir.mkdir(parents=True, exist_ok=True)

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    symbols = universe["symbol"].astype(str).tolist()
    if len(symbols) != 891 or universe["symbol"].duplicated().any():
        raise ValueError("Expected a unique 891-node universe")
    messenger_symbols = set(
        universe.loc[universe["node_type"].eq("molecule"), "symbol"]
    )
    mapping, cid_to_messenger = load_messenger_mapping(mapping_path)
    if set(mapping["messenger_node"]) != messenger_symbols:
        raise ValueError("PubChem mapping does not match the 20 messenger nodes")

    filtered, raw_rows = filter_stitch(raw_path, cid_to_messenger, args.chunk_size)
    collapsed = collapse_pairs(filtered)
    if collapsed.empty:
        raise ValueError("No curated messenger records were found in STITCH")

    string_mapping = pd.read_csv(
        string_mapping_path, sep="\t", dtype=str
    ).fillna("")
    mapped = string_mapping.loc[string_mapping["mapping_status"].eq("mapped")]
    if mapped["string_id"].duplicated().any():
        raise ValueError("STRING protein mapping is not one-to-one")
    stitch_to_node = mapped.set_index("string_id")["symbol"]
    collapsed["protein_node"] = collapsed["protein"].map(stitch_to_node).fillna("")
    collapsed["protein_in_node_universe"] = collapsed["protein_node"].ne("")
    collapsed["stitch_reference_score"] = STITCH_REPORTING_FLOOR
    collapsed["stitch_bayes_factor"] = stitch_bayes_factor(
        collapsed["stitch_score"].to_numpy(float)
    )
    collapsed["stitch_log_bayes_factor"] = np.log(
        collapsed["stitch_bayes_factor"].to_numpy(float)
    )
    collapsed = collapsed.sort_values(
        ["messenger_node", "stitch_combined_score", "protein"],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    all_mouse_path = processed_dir / "stitch_secondary_messenger_all_mouse_edges.tsv.gz"
    collapsed.to_csv(all_mouse_path, sep="\t", index=False, compression="gzip")

    in_universe = collapsed.loc[collapsed["protein_in_node_universe"]].copy()
    order = {symbol: index for index, symbol in enumerate(symbols)}
    pair_nodes = [
        sorted((protein, messenger), key=order.__getitem__)
        for protein, messenger in zip(
            in_universe["protein_node"],
            in_universe["messenger_node"],
            strict=True,
        )
    ]
    in_universe["node_a"] = [pair[0] for pair in pair_nodes]
    in_universe["node_b"] = [pair[1] for pair in pair_nodes]
    in_universe["edge_id"] = in_universe["node_a"] + "|" + in_universe["node_b"]
    if in_universe["edge_id"].duplicated().any():
        raise ValueError("Collapsed STITCH evidence contains duplicate universe pairs")

    current = pd.read_csv(current_matrix_path, sep="\t", index_col=0).reindex(
        index=symbols, columns=symbols
    )
    if current.isna().any().any():
        raise ValueError("Current adjacency matrix does not align to the universe")
    current_values = current.to_numpy(float)
    if not np.allclose(current_values, current_values.T, atol=1e-12):
        raise ValueError("Current adjacency matrix is not symmetric")
    final_values = current_values.copy()
    in_universe["pre_stitch_probability"] = np.nan
    in_universe["final_probability"] = np.nan
    for row_index, row in in_universe.iterrows():
        i = order[str(row["node_a"])]
        j = order[str(row["node_b"])]
        prior = float(current_values[i, j])
        posterior = float(
            stable_expit(
                np.asarray(
                    [stable_logit(np.asarray([prior]))[0] + row["stitch_log_bayes_factor"]]
                )
            )[0]
        )
        final_values[i, j] = final_values[j, i] = posterior
        in_universe.at[row_index, "pre_stitch_probability"] = prior
        in_universe.at[row_index, "final_probability"] = posterior
    np.fill_diagonal(final_values, 0.0)

    evidence_columns = [
        "edge_id",
        "node_a",
        "node_b",
        "messenger_node",
        "protein_node",
        "protein",
        "pubchem_cids",
        "stitch_chemical_ids",
        "chemical_representations",
        "supporting_stitch_rows",
        "stitch_combined_score",
        "stitch_score",
        "stitch_reference_score",
        "stitch_bayes_factor",
        "stitch_log_bayes_factor",
        "pre_stitch_probability",
        "final_probability",
    ]
    evidence = in_universe[evidence_columns].sort_values(
        ["stitch_bayes_factor", "node_a", "node_b"],
        ascending=[False, True, True],
    )
    evidence_path = output_dir / "stitch_secondary_messenger_edges.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False, float_format="%.12g")
    factor_catalog = evidence.loc[
        evidence["stitch_bayes_factor"] > 1.0 + EPSILON
    ].copy()
    factor_path = backend_dir / "stitch_secondary_messenger_bf_gt1.tsv.gz"
    factor_catalog.to_csv(
        factor_path, sep="\t", index=False, compression="gzip", float_format="%.12g"
    )

    factor_matrix = np.ones_like(final_values)
    np.fill_diagonal(factor_matrix, 0.0)
    for row in evidence.to_dict(orient="records"):
        i = order[str(row["node_a"])]
        j = order[str(row["node_b"])]
        factor_matrix[i, j] = factor_matrix[j, i] = float(row["stitch_bayes_factor"])
    write_matrix(output_dir / "stitch_bayes_factor_matrix.tsv", symbols, factor_matrix)
    write_matrix(output_dir / "combined_adjacency_matrix.tsv", symbols, final_values)

    audit_rows: list[dict[str, object]] = []
    for row in mapping.to_dict(orient="records"):
        messenger = str(row["messenger_node"])
        subset = collapsed.loc[collapsed["messenger_node"].eq(messenger)]
        universe_subset = subset.loc[subset["protein_in_node_universe"]]
        audit_rows.append(
            {
                **row,
                "found_in_stitch": bool(len(subset)),
                "unique_mouse_stitch_proteins": int(subset["protein"].nunique()),
                "in_universe_protein_nodes": int(universe_subset["protein_node"].nunique()),
                "in_universe_edges_above_reference": int(
                    (universe_subset["stitch_bayes_factor"] > 1.0 + EPSILON).sum()
                ),
                "maximum_stitch_score": (
                    float(subset["stitch_score"].max()) if len(subset) else np.nan
                ),
            }
        )
    audit = pd.DataFrame(audit_rows)
    audit_path = processed_dir / "secondary_messenger_stitch_coverage.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)

    upper = final_values[np.triu_indices(len(symbols), 1)]
    summary = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": "STITCH v5 mouse protein-chemical links",
        "source_url": SOURCE_URL,
        "source_file": str(raw_path),
        "source_sha256": sha256_file(raw_path),
        "source_file_bytes": raw_path.stat().st_size,
        "raw_stitch_rows_scanned": int(raw_rows),
        "curated_messenger_nodes": int(len(mapping)),
        "messengers_found_by_exact_pubchem_cid": int(
            audit["found_in_stitch"].sum()
        ),
        "messengers_with_in_universe_edges": int(
            (audit["in_universe_protein_nodes"] > 0).sum()
        ),
        "filtered_messenger_stitch_rows": int(len(filtered)),
        "unique_messenger_mouse_protein_pairs": int(len(collapsed)),
        "in_universe_messenger_protein_pairs": int(len(evidence)),
        "in_universe_pairs_above_reference": int(len(factor_catalog)),
        "stitch_score_minimum": float(evidence["stitch_score"].min()),
        "stitch_score_median": float(evidence["stitch_score"].median()),
        "stitch_score_maximum": float(evidence["stitch_score"].max()),
        "reference_score": STITCH_REPORTING_FLOOR,
        "bayes_factor_formula": (
            "max(1, odds(STITCH combined score / 1000) / odds(0.150))"
        ),
        "duplicate_rule": (
            "Collapse CIDm/CIDs and alternate PubChem records to one unordered "
            "messenger-protein pair using maximum combined confidence."
        ),
        "missing_evidence_rule": "Missing STITCH records remain neutral at BF=1.",
        "direction_rule": "All retained chemical-protein edges are treated as undirected.",
        "protein_mapping": (
            "STITCH v5 mouse Ensembl protein ID matched to the project's audited "
            "STRING v12 node mapping."
        ),
        "combined_unique_pairs_above_0_5": int((upper > 0.5).sum()),
        "limitations": [
            "STITCH v5 is the latest available STITCH bulk release but dates to 2016.",
            "Combined scores can include experiments, curated databases, predictions, and text mining.",
            "Generic lipid messengers are mapped to one explicitly audited PubChem class record.",
            "An absent STITCH record is not treated as evidence of no interaction.",
        ],
        "outputs": {
            "all_mouse_filtered_edges": str(all_mouse_path),
            "coverage_audit": str(audit_path),
            "in_universe_edges": str(evidence_path),
            "factor_catalog": str(factor_path),
            "factor_matrix": str(output_dir / "stitch_bayes_factor_matrix.tsv"),
            "combined_adjacency_matrix": str(output_dir / "combined_adjacency_matrix.tsv"),
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    readme = f"""STITCH secondary-messenger edge integration
================================================

This directory adds undirected protein-chemical evidence for the project's 20
curated second messengers. STRING itself is protein-only; STITCH v{STITCH_VERSION}
is its protein-chemical sister database.

Mapping
-------
Messenger nodes were mapped to explicitly audited PubChem CIDs in
secondary_messenger_pubchem_mapping.tsv. Both STITCH merged (CIDm) and
stereospecific (CIDs) representations were accepted. Duplicate representations
were collapsed to one messenger-protein pair by maximum combined score. STITCH
mouse protein IDs were mapped through the existing STRING node mapping.

Evidence conversion
-------------------
STITCH stores confidence scores from 0.150 to 0.999. The reporting floor 0.150
is the neutral reference:

  BF = max(1, odds(score / 1000) / odds(0.150))

Missing pairs remain BF=1. Edges are undirected. The output retains every source
score and identifier needed for review. STITCH evidence is association evidence,
not a claim of direct binding or direction of signal flow.
"""
    (output_dir / "README.txt").write_text(readme, encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
