#!/usr/bin/env python3
"""Extend the current undirected graph with one external mouse protein target.

The existing square adjacency block is preserved exactly. The function
calculates one edge probability between the target and every existing node by
applying the same five evidence streams used in edge characterization:

1. mpkCCD differential-fraction localization;
2. observed/curated-site KinasePredictor;
3. STRING v12 functional associations;
4. HPA v25.1 primary localization; and
5. OmniPath core signaling interactions, collapsed to undirected evidence.

Missing coverage is neutral for every stream. The output vector is in the
canonical node-universe order and can optionally be appended as one row and
column to create a symmetric matrix with one additional node.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
CODE_DIR = SCRIPT_DIR.parent
EDGE_CODE_DIR = CODE_DIR / "edge_characterization"
for search_path in (CODE_DIR, EDGE_CODE_DIR):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from build_localization_adjacency import (  # noqa: E402
    FRACTIONS,
    collapse_source_profiles,
)
from build_phosphosite_database import (  # noqa: E402
    MOD_RES_RE,
    SITE_RE,
    complete_window,
    entry_rank,
    load_uniprot,
    normalized_sequence,
    split_gene_names,
)
from integrate_hpa_localization_edges import (  # noqa: E402
    MINIMUM_LIKELIHOOD as HPA_MINIMUM_LIKELIHOOD,
    PRIMARY_TIERS,
    build_profiles,
    extract_location_dictionary,
    load_hpa,
    load_one_to_one_orthology,
)
from integrate_observed_phosphosite_kinase_edges import (  # noqa: E402
    AA_INDEX,
    AA_ORDER,
    MINIMUM_SITE_LIKELIHOOD,
    TOP_N,
    load_matrix_set,
    load_metadata,
    map_models,
    score_site,
)
from integrate_string_edges import (  # noqa: E402
    CHANNEL_COLUMNS,
    STRING_PRIOR,
    choose_mapping,
    collect_alias_hits,
    read_string_info,
)


PROJECT_DEFAULT = Path(__file__).resolve().parents[2]
Q = 0.75
NEUTRAL_LIKELIHOOD = 0.5
NEUTRAL_BAYES_FACTOR = 1.0


@dataclass
class TargetAdjacencyResult:
    """In-memory result returned by :func:`build_target_adjacency_vector`."""

    target_symbol: str
    target_uniprot: str
    vector: pd.DataFrame
    extended_matrix: pd.DataFrame | None
    phosphosites: pd.DataFrame
    kinase_predictions: pd.DataFrame
    summary: dict[str, Any]
    output_directory: Path | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "target",
        help="Mouse gene symbol or UniProt accession outside the node universe.",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_DEFAULT,
        help=f"Project directory (default: {PROJECT_DEFAULT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional output directory; defaults to results/path_finding/target_extensions/<symbol>.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=TOP_N,
        help=f"Top distinct in-universe kinases retained per phosphosite (default: {TOP_N}).",
    )
    parser.add_argument(
        "--no-extended-matrix",
        action="store_true",
        help="Write only the adjacency vector, not the appended square matrix.",
    )
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "-"} else text


def split_semicolon(value: object) -> set[str]:
    return {item.strip() for item in clean(value).split(";") if item.strip()}


def parse_bool(value: object) -> bool:
    return clean(value).casefold() in {"1", "true", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._") or "target"


def stable_logit(values: np.ndarray | pd.Series) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    clipped = np.clip(array, 1e-15, 1.0 - 1e-15)
    return np.log(clipped / (1.0 - clipped))


def stable_expit(values: np.ndarray | pd.Series) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return 1.0 / (1.0 + np.exp(-np.clip(array, -700.0, 700.0)))


def apply_bayes_factor(
    prior: np.ndarray,
    bayes_factor: np.ndarray,
) -> np.ndarray:
    """Apply an edge-wise factor while preserving BF=1 entries exactly."""

    prior = np.asarray(prior, dtype=float)
    bayes_factor = np.asarray(bayes_factor, dtype=float)
    if prior.shape != bayes_factor.shape:
        raise ValueError("prior and bayes_factor shapes differ")
    if (
        not np.isfinite(prior).all()
        or not np.isfinite(bayes_factor).all()
        or (prior < 0).any()
        or (prior > 1).any()
        or (bayes_factor <= 0).any()
    ):
        raise ValueError("invalid prior or Bayes-factor values")
    result = prior.copy()
    changed = ~np.isclose(
        bayes_factor,
        NEUTRAL_BAYES_FACTOR,
        atol=1e-15,
        rtol=0.0,
    )
    result[changed] = stable_expit(
        stable_logit(prior[changed]) + np.log(bayes_factor[changed])
    )
    return result


def complement_minimum(values: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    thresholds = np.asarray(thresholds, dtype=float)
    factors = np.full(values.shape, NEUTRAL_LIKELIHOOD, dtype=float)
    valid = (
        np.isfinite(values)
        & np.isfinite(thresholds)
        & (thresholds > 0)
    )
    if valid.any():
        z = values[valid] / thresholds[valid]
        factors[valid] = np.maximum(
            NEUTRAL_LIKELIHOOD,
            1.0 - np.exp(-0.5 * np.square(z)),
        )
    return factors


def resolve_target(
    identifier: str,
    uniprot_path: Path,
) -> tuple[str, dict[str, str], str]:
    by_accession, by_primary, by_synonym = load_uniprot(uniprot_path)
    key = clean(identifier).upper()
    if not key:
        raise ValueError("target identifier is empty")

    if key in by_accession:
        entry = by_accession[key]
        method = "exact_uniprot_accession"
    elif by_primary.get(key):
        entry = max(
            by_primary[key],
            key=lambda item: entry_rank(item, key),
        )
        method = "exact_primary_gene_symbol"
    elif by_synonym.get(key):
        entry = max(
            by_synonym[key],
            key=lambda item: entry_rank(item, key),
        )
        method = "exact_gene_synonym"
    else:
        raise ValueError(
            f"Target {identifier!r} was not found in the cached mouse "
            "UniProt reference proteome."
        )

    primary_names = split_gene_names(entry["Gene Names (primary)"])
    if not primary_names:
        raise ValueError(
            f"Resolved UniProt entry {entry['Entry']} lacks a primary gene symbol"
        )
    return primary_names[0], entry, method


def target_localization_evidence(
    project: Path,
    target_symbol: str,
    universe: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    profiles_path = (
        project
        / "data/edge_characterization/localization/processed"
        / "node_localization_profiles.tsv"
    )
    source_path = (
        project
        / "data/edge_characterization/localization/raw"
        / "Proteomics_of_subcellular_fractions.xlsx"
    )
    existing = pd.read_csv(profiles_path, sep="\t", dtype=str).fillna("")
    if existing["symbol"].tolist() != universe["symbol"].tolist():
        raise ValueError("mpkCCD profile order does not match node universe")

    raw = pd.read_excel(
        source_path,
        sheet_name="Web Page w Links",
        header=3,
        usecols="A,C,E:I,Q",
    )
    collapsed = collapse_source_profiles(raw)
    target_key = target_symbol.upper()
    matched = collapsed.loc[target_key] if target_key in collapsed.index else None
    target_profile = (
        np.asarray([float(matched[item]) for item in FRACTIONS], dtype=float)
        if matched is not None
        else np.full(len(FRACTIONS), np.nan)
    )
    target_observed = bool(
        matched is not None
        and np.isfinite(target_profile).all()
        and target_profile.sum() > 0
    )

    n = len(universe)
    observed_existing = existing["localization_observed"].map(parse_bool).to_numpy()
    existing_profiles = existing[FRACTIONS].apply(
        pd.to_numeric,
        errors="coerce",
    ).to_numpy(dtype=float)
    node_tq = pd.to_numeric(
        existing["localization_Tq"],
        errors="coerce",
    ).to_numpy(dtype=float)
    dot_product = np.full(n, np.nan, dtype=float)
    factor_target = np.full(n, NEUTRAL_LIKELIHOOD, dtype=float)
    factor_node = np.full(n, NEUTRAL_LIKELIHOOD, dtype=float)
    likelihood = np.full(n, NEUTRAL_LIKELIHOOD, dtype=float)
    target_tq = math.nan

    if target_observed and observed_existing.any():
        dot_product[observed_existing] = (
            existing_profiles[observed_existing] @ target_profile
        )
        target_tq = float(np.quantile(dot_product[observed_existing], Q))
        if not math.isfinite(target_tq) or target_tq <= 0:
            raise ValueError(
                f"Target mpkCCD localization Tq is not positive: {target_tq}"
            )
        target_thresholds = np.full(n, target_tq, dtype=float)
        factor_target[observed_existing] = complement_minimum(
            dot_product[observed_existing],
            target_thresholds[observed_existing],
        )
        factor_node[observed_existing] = complement_minimum(
            dot_product[observed_existing],
            node_tq[observed_existing],
        )
        likelihood[observed_existing] = (
            factor_target[observed_existing]
            + factor_node[observed_existing]
        ) / 2.0

    evidence = pd.DataFrame(
        {
            "mpkccd_localization_observed_for_target": target_observed,
            "mpkccd_localization_observed_for_node": observed_existing,
            "mpkccd_localization_observed_for_both": (
                target_observed & observed_existing
            ),
            "mpkccd_dot_product": dot_product,
            "mpkccd_target_Tq": target_tq,
            "mpkccd_node_Tq": node_tq,
            "mpkccd_factor_using_target_Tq": factor_target,
            "mpkccd_factor_using_node_Tq": factor_node,
            "mpkccd_likelihood": likelihood,
            "mpkccd_bayes_factor": likelihood / NEUTRAL_LIKELIHOOD,
            "mpkccd_log_bayes_factor": np.log(
                likelihood / NEUTRAL_LIKELIHOOD
            ),
        }
    )
    metadata = {
        "observed": target_observed,
        "source_gene_symbol": (
            clean(matched["source_gene_symbol"])
            if matched is not None
            else ""
        ),
        "selected_gi_number": (
            clean(matched["selected_gi_number"])
            if matched is not None
            else ""
        ),
        "source_entry_count": (
            int(matched["source_entry_count"])
            if matched is not None
            else 0
        ),
        "distinct_profile_count": (
            int(matched["distinct_profile_count"])
            if matched is not None
            else 0
        ),
        "fractions": {
            name: (
                float(matched[name]) if matched is not None else None
            )
            for name in FRACTIONS
        },
        "Tq": target_tq if math.isfinite(target_tq) else None,
    }
    return evidence, metadata


def target_phosphosites(
    project: Path,
    target_symbol: str,
    target_entry: dict[str, str],
    audit_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    sequence = normalized_sequence(target_entry["Sequence"])
    target_uniprot = clean(target_entry["Entry"])
    records: dict[tuple[str, int], dict[str, object]] = {}

    def ensure(residue: str, position: int) -> dict[str, object]:
        key = (residue, position)
        if key not in records:
            canonical = complete_window(sequence, position)
            records[key] = {
                "target_symbol": target_symbol,
                "target_uniprot": target_uniprot,
                "site": f"{residue}{position}",
                "position": position,
                "residue": residue,
                "sources": set(),
                "modifications": set(),
                "evidence": set(),
                "canonical_window": canonical,
                "dataset_windows": set(),
                "dataset_raw_log2_change": set(),
                "dataset_site_p_value": set(),
                "dataset_site_key": set(),
            }
        return records[key]

    for match in MOD_RES_RE.finditer(clean(target_entry["Modified residue"])):
        start_text, end_text, note, evidence = match.groups()
        if "phospho" not in note.casefold():
            continue
        position = int(start_text)
        if end_text and int(end_text) != position:
            continue
        if not 1 <= position <= len(sequence):
            continue
        residue = sequence[position - 1]
        record = ensure(residue, position)
        record["sources"].add("uniprot_annotated")
        record["modifications"].add(note)
        if evidence:
            record["evidence"].add(evidence)

    audit_path = project / "results/phosphoprotein_evidence/phosphosite_audit.tsv"
    if audit_path.exists():
        audit = (
            audit_table.fillna("")
            if audit_table is not None
            else pd.read_csv(audit_path, sep="\t", dtype=str).fillna("")
        )
        target_rows = audit.loc[
            audit["gene_symbol"].str.casefold().eq(target_symbol.casefold())
            | audit["uniprot"].str.casefold().eq(target_uniprot.casefold())
        ]
        for row in target_rows.to_dict(orient="records"):
            match = SITE_RE.match(clean(row.get("site", "")))
            if match is None:
                continue
            residue = match.group(1).upper()
            position = int(match.group(2))
            record = ensure(residue, position)
            record["sources"].add("dataset_observed")
            record["modifications"].add(
                {
                    "S": "Phosphoserine",
                    "T": "Phosphothreonine",
                    "Y": "Phosphotyrosine",
                }.get(residue, "Phosphorylated residue")
            )
            record["evidence"].add("PKA-KO phosphoproteomics")
            dataset_window = normalized_sequence(
                row.get("centralized_sequence", "")
            )
            if len(dataset_window) == 13:
                record["dataset_windows"].add(dataset_window)
            for field, destination in (
                ("raw_log2_change", "dataset_raw_log2_change"),
                ("site_p_value", "dataset_site_p_value"),
                ("site_key", "dataset_site_key"),
            ):
                value = clean(row.get(field, ""))
                if value:
                    record[destination].add(value)

    rows: list[dict[str, object]] = []
    for (_, _), record in sorted(
        records.items(),
        key=lambda item: (item[0][1], item[0][0]),
    ):
        residue = str(record["residue"])
        canonical = str(record["canonical_window"])
        dataset_windows = sorted(record["dataset_windows"])
        selected_window = canonical if len(canonical) == 13 else ""
        window_source = "canonical_uniprot" if selected_window else ""
        if not selected_window:
            for candidate in dataset_windows:
                if candidate[6] == residue:
                    selected_window = candidate
                    window_source = "dataset_supplied"
                    break
        unsupported = sorted(set(selected_window).difference(AA_INDEX))
        scorable = bool(
            len(selected_window) == 13
            and selected_window[6] == residue
            and residue in {"S", "T", "Y"}
            and not unsupported
        )

        def join_set(name: str) -> str:
            return ";".join(sorted(str(value) for value in record[name]))

        rows.append(
            {
                "target_symbol": target_symbol,
                "target_uniprot": target_uniprot,
                "site": record["site"],
                "position": record["position"],
                "residue": residue,
                "sources": join_set("sources"),
                "modification": join_set("modifications"),
                "evidence": join_set("evidence"),
                "canonical_sequence_13mer": canonical,
                "dataset_sequence_13mers": ";".join(dataset_windows),
                "centralized_sequence_13mer": selected_window,
                "window_source": window_source,
                "kinasepredictor_scorable": scorable,
                "scoring_exclusion_reason": (
                    ""
                    if scorable
                    else (
                        "no_complete_supported_13mer"
                        if not selected_window
                        else "unsupported_or_mismatched_sequence"
                    )
                ),
                "dataset_raw_log2_change": join_set(
                    "dataset_raw_log2_change"
                ),
                "dataset_site_p_value": join_set(
                    "dataset_site_p_value"
                ),
                "dataset_site_key": join_set("dataset_site_key"),
            }
        )
    return pd.DataFrame(rows)


def target_kinase_evidence(
    project: Path,
    universe: pd.DataFrame,
    phosphosites: pd.DataFrame,
    top_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    n = len(universe)
    symbol_to_index = {
        symbol: index
        for index, symbol in enumerate(universe["symbol"].astype(str))
    }
    universe_ci = {
        symbol.upper(): symbol for symbol in universe["symbol"].astype(str)
    }
    kinase_nodes = {
        str(row["symbol"])
        for row in universe.to_dict(orient="records")
        if "kinase" in clean(row.get("classes", "")).split(";")
    }
    matrix_root = (
        project
        / "data/kinase_predictor/v0.8/official_package/output_matrices"
    )
    metadata_path = (
        project
        / "data/edge_characterization/kinase_predictor"
        / "kinase_logos_metadata.html"
    )
    st_labels, st_matrices = load_matrix_set(
        matrix_root / "Ser-Thr_output_matrices"
    )
    tyr_labels, tyr_matrices = load_matrix_set(
        matrix_root / "Tyr_output_matrices"
    )
    metadata = load_metadata(metadata_path)
    st_mapping = map_models(
        "serine_threonine",
        st_labels,
        metadata,
        universe_ci,
        kinase_nodes,
    )
    tyr_mapping = map_models(
        "tyrosine",
        tyr_labels,
        metadata,
        universe_ci,
        kinase_nodes,
    )

    log_bf = np.zeros(n, dtype=float)
    hit_count = np.zeros(n, dtype=int)
    site_sets: list[set[str]] = [set() for _ in range(n)]
    label_sets: list[set[str]] = [set() for _ in range(n)]
    prediction_rows: list[dict[str, object]] = []
    scored_sites = 0
    fallback_sites = 0

    for site in phosphosites.to_dict(orient="records"):
        if not bool(site.get("kinasepredictor_scorable", False)):
            continue
        sequence = clean(site["centralized_sequence_13mer"]).upper()
        residue = clean(site["residue"]).upper()
        if len(sequence) != 13 or sequence[6] != residue:
            continue
        if residue in {"S", "T"}:
            labels, matrices, mapping = (
                st_labels,
                st_matrices,
                st_mapping,
            )
            matrix_set = "serine_threonine"
        elif residue == "Y":
            labels, matrices, mapping = (
                tyr_labels,
                tyr_matrices,
                tyr_mapping,
            )
            matrix_set = "tyrosine"
        else:
            continue
        result = score_site(sequence, labels, matrices, mapping, top_n)
        scored_sites += 1
        if result["threshold_method"] == "positive_score_q75_fallback":
            fallback_sites += 1
        for prediction in result["predictions"]:
            kinase = clean(prediction["mapped_mouse_node"])
            if kinase not in symbol_to_index:
                continue
            index = symbol_to_index[kinase]
            site_bf = float(prediction["site_bayes_factor"])
            log_bf[index] += math.log(site_bf)
            hit_count[index] += 1
            site_sets[index].add(clean(site["site"]))
            label_sets[index].add(clean(prediction["predictor_label"]))
            prediction_rows.append(
                {
                    "target_symbol": site["target_symbol"],
                    "target_site": site["site"],
                    "position": site["position"],
                    "residue": residue,
                    "observed_sources": site["sources"],
                    "centralized_sequence_13mer": sequence,
                    "matrix_set": matrix_set,
                    "kinase_node": kinase,
                    "predictor_label": prediction["predictor_label"],
                    "predictor_family": prediction["predictor_family"],
                    "in_universe_rank": prediction["in_universe_rank"],
                    "global_model_rank": prediction["global_model_rank"],
                    "raw_score": prediction["raw_score"],
                    "raw_site_Tq_q75": result["raw_threshold"],
                    "site_Tq_q75": result["threshold"],
                    "Tq_method": result["threshold_method"],
                    "site_edge_likelihood": prediction[
                        "site_edge_likelihood"
                    ],
                    "site_bayes_factor": site_bf,
                }
            )

    bayes_factor = np.exp(np.clip(log_bf, -700.0, 700.0))
    evidence = pd.DataFrame(
        {
            "kinasepredictor_evidence_observed": hit_count > 0,
            "kinasepredictor_prediction_hit_count": hit_count,
            "kinasepredictor_unique_target_site_count": [
                len(values) for values in site_sets
            ],
            "kinasepredictor_supporting_sites": [
                ";".join(sorted(values)) for values in site_sets
            ],
            "kinasepredictor_predictor_labels": [
                ";".join(sorted(values)) for values in label_sets
            ],
            "kinasepredictor_bayes_factor": bayes_factor,
            "kinasepredictor_log_bayes_factor": log_bf,
        }
    )
    predictions = pd.DataFrame(prediction_rows)
    summary = {
        "observed_or_annotated_phosphosites": int(len(phosphosites)),
        "scorable_phosphosites": int(
            phosphosites.get(
                "kinasepredictor_scorable",
                pd.Series(dtype=bool),
            ).sum()
        ),
        "sites_scored": scored_sites,
        "sites_using_positive_q75_fallback": fallback_sites,
        "predictions_emitted": int(len(predictions)),
        "kinase_nodes_with_ranked_hits": int((hit_count > 0).sum()),
        "kinase_nodes_with_non_neutral_factor": int(
            (bayes_factor > 1.0 + 1e-12).sum()
        ),
        "top_n_per_site": top_n,
    }
    return evidence, predictions, summary


def target_string_evidence(
    project: Path,
    universe: pd.DataFrame,
    target_symbol: str,
    target_uniprot: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    n = len(universe)
    mapping_path = (
        project
        / "data/edge_characterization/string/v12.0/processed"
        / "node_to_string_mapping.tsv"
    )
    raw_dir = project / "data/edge_characterization/string/v12.0/raw"
    info_path = raw_dir / "10090.protein.info.v12.0.txt.gz"
    alias_path = raw_dir / "10090.protein.aliases.v12.0.txt.gz"
    links_path = raw_dir / "10090.protein.links.detailed.v12.0.txt.gz"

    existing_mapping = pd.read_csv(
        mapping_path,
        sep="\t",
        dtype=str,
    ).fillna("")
    existing_id_to_symbol = {
        row["string_id"]: row["symbol"]
        for row in existing_mapping.to_dict(orient="records")
        if row["mapping_status"] == "mapped" and row["string_id"]
    }
    info = read_string_info(info_path)
    alias_hits = collect_alias_hits(
        alias_path,
        {target_symbol, target_uniprot},
    )
    target_mapping = choose_mapping(
        target_symbol,
        [target_uniprot],
        info,
        alias_hits,
    )
    target_string_id = clean(target_mapping["string_id"])

    records: dict[str, dict[str, object]] = {}
    raw_rows_with_target = 0
    expected_header = [
        "protein1",
        "protein2",
        "neighborhood",
        "fusion",
        "cooccurence",
        "coexpression",
        "experimental",
        "database",
        "textmining",
        "combined_score",
    ]
    if target_string_id:
        with gzip.open(links_path, "rt", encoding="utf-8") as handle:
            header = handle.readline().strip().split()
            if header != expected_header:
                raise ValueError(f"Unexpected STRING links header: {header}")
            for line in handle:
                fields = line.split()
                if len(fields) != 10:
                    continue
                protein_a, protein_b = fields[:2]
                if protein_a == target_string_id:
                    partner_id = protein_b
                elif protein_b == target_string_id:
                    partner_id = protein_a
                else:
                    continue
                raw_rows_with_target += 1
                partner = existing_id_to_symbol.get(partner_id)
                if not partner:
                    continue
                scores = [int(value) / 1000.0 for value in fields[2:]]
                record = {
                    "node_symbol": partner,
                    "target_string_id": target_string_id,
                    "node_string_id": partner_id,
                    "neighborhood_score": scores[0],
                    "fusion_score": scores[1],
                    "cooccurrence_score": scores[2],
                    "coexpression_score": scores[3],
                    "experimental_score": scores[4],
                    "database_score": scores[5],
                    "textmining_score": scores[6],
                    "string_combined_score": scores[7],
                }
                previous = records.get(partner)
                if (
                    previous is None
                    or record["string_combined_score"]
                    > previous["string_combined_score"]
                ):
                    records[partner] = record

    symbol_to_index = {
        symbol: index
        for index, symbol in enumerate(universe["symbol"].astype(str))
    }
    observed = np.zeros(n, dtype=bool)
    score = np.full(n, np.nan, dtype=float)
    log_bf = np.zeros(n, dtype=float)
    target_id = np.full(n, "", dtype=object)
    node_id = np.full(n, "", dtype=object)
    channels = {
        column: np.full(n, np.nan, dtype=float)
        for column in CHANNEL_COLUMNS
    }
    prior_log_odds = math.log(STRING_PRIOR / (1.0 - STRING_PRIOR))
    for node, record in records.items():
        index = symbol_to_index[node]
        observed[index] = True
        combined = float(record["string_combined_score"])
        score[index] = combined
        log_bf[index] = math.log(combined / (1.0 - combined)) - prior_log_odds
        target_id[index] = record["target_string_id"]
        node_id[index] = record["node_string_id"]
        for column in CHANNEL_COLUMNS:
            channels[column][index] = float(record[column])
    bayes_factor = np.exp(np.clip(log_bf, -700.0, 700.0))

    evidence_data: dict[str, object] = {
        "string_evidence_observed": observed,
        "string_target_id": target_id,
        "string_node_id": node_id,
        "string_combined_score": score,
        "string_bayes_factor": bayes_factor,
        "string_log_bayes_factor": log_bf,
    }
    evidence_data.update({f"string_{key}": value for key, value in channels.items()})
    evidence = pd.DataFrame(evidence_data)
    summary = {
        "target_mapping": target_mapping,
        "target_string_id": target_string_id,
        "raw_network_rows_incident_to_target": raw_rows_with_target,
        "in_universe_unique_associations": int(observed.sum()),
        "non_neutral_associations": int((bayes_factor > 1.0 + 1e-12).sum()),
        "mapping_rule": (
            "Exact target UniProt plus preferred-name mapping, using the "
            "same audited STRING mapper as the current edge workflow."
        ),
    }
    return evidence, summary


def target_hpa_evidence(
    project: Path,
    universe: pd.DataFrame,
    target_symbol: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw_dir = project / "data/edge_characterization/localization/hpa/v25.1/raw"
    hpa_path = raw_dir / "subcellular_location.tsv.zip"
    orthology_path = raw_dir / "HOM_ProteinCoding.rpt"
    node_profiles_path = (
        project
        / "data/edge_characterization/localization/hpa/v25.1/processed"
        / "node_hpa_localization_profiles.tsv"
    )
    hpa = load_hpa(hpa_path)
    _, orthology_lookup, ambiguous_mouse = load_one_to_one_orthology(
        orthology_path
    )
    location_dictionary = extract_location_dictionary(hpa)
    locations = sorted(
        location_dictionary.loc[
            location_dictionary["included_in_primary_profile"],
            "location",
        ].tolist()
    )
    target_universe = pd.DataFrame(
        [
            {
                "symbol": target_symbol,
                "display_symbol": target_symbol,
                "node_type": "protein",
            }
        ]
    )
    target_audit, target_binary, _ = build_profiles(
        target_universe,
        hpa,
        orthology_lookup,
        ambiguous_mouse,
        locations,
    )
    target_row = target_audit.iloc[0].to_dict()
    target_observed = bool(target_binary[0].sum() > 0)

    existing = pd.read_csv(
        node_profiles_path,
        sep="\t",
        dtype=str,
    ).fillna("")
    if existing["symbol"].tolist() != universe["symbol"].tolist():
        raise ValueError("HPA profile order does not match node universe")
    profile_columns = []
    for location in locations:
        safe = re.sub(r"[^A-Za-z0-9]+", "_", location).strip("_").lower()
        profile_columns.append(f"location__{safe}")
    existing_binary = existing[profile_columns].apply(
        pd.to_numeric,
        errors="coerce",
    ).fillna(0).to_numpy(dtype=float)
    existing_observed = existing_binary.sum(axis=1) > 0
    node_tq = pd.to_numeric(
        existing["hpa_Tq_q75"],
        errors="coerce",
    ).to_numpy(dtype=float)
    n = len(universe)
    similarity = np.full(n, np.nan, dtype=float)
    factor_target = np.full(n, HPA_MINIMUM_LIKELIHOOD, dtype=float)
    factor_node = np.full(n, HPA_MINIMUM_LIKELIHOOD, dtype=float)
    likelihood = np.full(n, HPA_MINIMUM_LIKELIHOOD, dtype=float)
    target_raw_tq = math.nan
    target_tq = math.nan
    threshold_method = "profile_not_observed"
    positive_background_size = 0

    if target_observed and existing_observed.any():
        target_unit = target_binary[0].astype(float)
        target_unit /= np.linalg.norm(target_unit)
        existing_unit = existing_binary[existing_observed]
        existing_unit /= np.linalg.norm(
            existing_unit,
            axis=1,
            keepdims=True,
        )
        similarity[existing_observed] = existing_unit @ target_unit
        background = similarity[existing_observed]
        target_raw_tq = float(np.quantile(background, Q))
        positive = background[background > 0]
        positive_background_size = int(len(positive))
        if math.isfinite(target_raw_tq) and target_raw_tq > 0:
            target_tq = target_raw_tq
            threshold_method = "all_observed_nodes_q75"
        elif len(positive):
            target_tq = float(np.quantile(positive, Q))
            threshold_method = "positive_overlap_q75_fallback"
        else:
            threshold_method = "no_positive_overlap_neutral"
        if math.isfinite(target_tq) and target_tq > 0:
            factor_target[existing_observed] = complement_minimum(
                similarity[existing_observed],
                np.full(existing_observed.sum(), target_tq),
            )
            factor_node[existing_observed] = complement_minimum(
                similarity[existing_observed],
                node_tq[existing_observed],
            )
            likelihood[existing_observed] = (
                factor_target[existing_observed]
                + factor_node[existing_observed]
            ) / 2.0

    evidence = pd.DataFrame(
        {
            "hpa_profile_observed_for_target": target_observed,
            "hpa_profile_observed_for_node": existing_observed,
            "hpa_profiles_observed_for_both": (
                target_observed & existing_observed
            ),
            "hpa_cosine_similarity": similarity,
            "hpa_target_Tq": target_tq,
            "hpa_target_Tq_method": threshold_method,
            "hpa_node_Tq": node_tq,
            "hpa_factor_using_target_Tq": factor_target,
            "hpa_factor_using_node_Tq": factor_node,
            "hpa_likelihood": likelihood,
            "hpa_bayes_factor": likelihood / NEUTRAL_LIKELIHOOD,
            "hpa_log_bayes_factor": np.log(
                likelihood / NEUTRAL_LIKELIHOOD
            ),
        }
    )
    metadata = {
        **target_row,
        "hpa_raw_Tq_q75": (
            target_raw_tq if math.isfinite(target_raw_tq) else None
        ),
        "hpa_Tq_q75": target_tq if math.isfinite(target_tq) else None,
        "hpa_Tq_method": threshold_method,
        "hpa_positive_overlap_background_size": positive_background_size,
        "included_tiers": list(PRIMARY_TIERS),
    }
    return evidence, metadata


def target_omnipath_evidence(
    project: Path,
    universe: pd.DataFrame,
    target_symbol: str,
    target_uniprot: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    raw_path = (
        project
        / "data/edge_characterization/omnipath/2026-07-30/raw"
        / "omnipath_mouse_core_post_translational.tsv"
    )
    summary_path = (
        project
        / "results/edge_characterization"
        / "localization_kinase_predictor_string_hpa_omnipath"
        / "analysis_summary.json"
    )
    protein_index_path = (
        project
        / "data/edge_characterization/kinase_predictor"
        / "phosphosite_database/protein_index.tsv"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    t_q = float(summary["Tq_curation_effort"])
    raw = pd.read_csv(raw_path, sep="\t", dtype=str).fillna("")
    raw["curation_effort"] = pd.to_numeric(
        raw["curation_effort"],
        errors="raise",
    ).astype(float)
    protein_index = pd.read_csv(
        protein_index_path,
        sep="\t",
        dtype=str,
    ).fillna("")
    universe_ci = {
        symbol.upper(): symbol
        for symbol in universe.loc[
            universe["node_type"].eq("protein"),
            "symbol",
        ].astype(str)
    }
    accession_to_symbol = {
        row["mapped_uniprot"]: row["symbol"]
        for row in protein_index.to_dict(orient="records")
        if row["mapped_uniprot"]
    }
    target_key = target_symbol.upper()

    def map_endpoint(api_symbol: str, accession: str) -> str:
        symbol_key = clean(api_symbol).upper()
        if symbol_key == target_key or clean(accession) == target_uniprot:
            return target_symbol
        if symbol_key in universe_ci:
            return universe_ci[symbol_key]
        return accession_to_symbol.get(clean(accession), "")

    raw["mapped_source"] = [
        map_endpoint(symbol, accession)
        for symbol, accession in zip(
            raw["source_genesymbol"],
            raw["source"],
            strict=True,
        )
    ]
    raw["mapped_target"] = [
        map_endpoint(symbol, accession)
        for symbol, accession in zip(
            raw["target_genesymbol"],
            raw["target"],
            strict=True,
        )
    ]
    incident = raw.loc[
        (
            raw["mapped_source"].eq(target_symbol)
            & raw["mapped_target"].isin(universe_ci.values())
        )
        | (
            raw["mapped_target"].eq(target_symbol)
            & raw["mapped_source"].isin(universe_ci.values())
        )
    ].copy()
    incident["node_symbol"] = np.where(
        incident["mapped_source"].eq(target_symbol),
        incident["mapped_target"],
        incident["mapped_source"],
    )
    for column in (
        "is_stimulation",
        "is_inhibition",
        "consensus_direction",
        "consensus_stimulation",
        "consensus_inhibition",
    ):
        incident[f"{column}_bool"] = incident[column].map(parse_bool)

    rows: list[dict[str, object]] = []
    for node, group in incident.groupby("node_symbol", sort=False):
        resources: set[str] = set()
        references: set[str] = set()
        for record in group.to_dict(orient="records"):
            resources.update(split_semicolon(record["sources"]))
            references.update(split_semicolon(record["references"]))
        directions = sorted(
            {
                f"{source}->{target}"
                for source, target in zip(
                    group["mapped_source"],
                    group["mapped_target"],
                    strict=True,
                )
            }
        )
        effort = float(group["curation_effort"].max())
        support_component = 1.0 - math.exp(
            -0.5 * (effort / t_q) ** 2
        )
        likelihood = (
            NEUTRAL_LIKELIHOOD
            + (1.0 - NEUTRAL_LIKELIHOOD) * support_component
        )
        rows.append(
            {
                "node_symbol": node,
                "omnipath_evidence_observed": True,
                "omnipath_directed_record_count": int(len(group)),
                "omnipath_directional_interactions": ";".join(directions),
                "omnipath_bidirectional": len(directions) > 1,
                "omnipath_any_stimulation_annotation": bool(
                    group["is_stimulation_bool"].any()
                ),
                "omnipath_any_inhibition_annotation": bool(
                    group["is_inhibition_bool"].any()
                ),
                "omnipath_any_unsigned_annotation": bool(
                    (
                        ~group["is_stimulation_bool"]
                        & ~group["is_inhibition_bool"]
                    ).any()
                ),
                "omnipath_any_consensus_direction": bool(
                    group["consensus_direction_bool"].any()
                ),
                "omnipath_curation_effort": effort,
                "omnipath_Tq": t_q,
                "omnipath_resources": ";".join(sorted(resources)),
                "omnipath_references": ";".join(sorted(references)),
                "omnipath_likelihood": likelihood,
                "omnipath_bayes_factor": (
                    likelihood / NEUTRAL_LIKELIHOOD
                ),
                "omnipath_log_bayes_factor": math.log(
                    likelihood / NEUTRAL_LIKELIHOOD
                ),
            }
        )

    columns = [
        "omnipath_evidence_observed",
        "omnipath_directed_record_count",
        "omnipath_directional_interactions",
        "omnipath_bidirectional",
        "omnipath_any_stimulation_annotation",
        "omnipath_any_inhibition_annotation",
        "omnipath_any_unsigned_annotation",
        "omnipath_any_consensus_direction",
        "omnipath_curation_effort",
        "omnipath_Tq",
        "omnipath_resources",
        "omnipath_references",
        "omnipath_likelihood",
        "omnipath_bayes_factor",
        "omnipath_log_bayes_factor",
    ]
    evidence = pd.DataFrame({"node_symbol": universe["symbol"]})
    supported = pd.DataFrame(rows)
    if supported.empty:
        supported = pd.DataFrame(columns=["node_symbol", *columns])
    evidence = evidence.merge(supported, on="node_symbol", how="left")
    defaults: dict[str, object] = {
        "omnipath_evidence_observed": False,
        "omnipath_directed_record_count": 0,
        "omnipath_directional_interactions": "",
        "omnipath_bidirectional": False,
        "omnipath_any_stimulation_annotation": False,
        "omnipath_any_inhibition_annotation": False,
        "omnipath_any_unsigned_annotation": False,
        "omnipath_any_consensus_direction": False,
        "omnipath_curation_effort": np.nan,
        "omnipath_Tq": t_q,
        "omnipath_resources": "",
        "omnipath_references": "",
        "omnipath_likelihood": NEUTRAL_LIKELIHOOD,
        "omnipath_bayes_factor": NEUTRAL_BAYES_FACTOR,
        "omnipath_log_bayes_factor": 0.0,
    }
    for column, default in defaults.items():
        evidence[column] = evidence[column].fillna(default)
    evidence = evidence.drop(columns="node_symbol")
    stream_summary = {
        "directed_rows_incident_to_target_and_universe": int(len(incident)),
        "unique_undirected_target_node_pairs": int(len(rows)),
        "Tq_curation_effort": t_q,
        "direction_and_sign_rule": (
            "Direction and stimulation/inhibition are retained for audit "
            "but ignored by the undirected model."
        ),
    }
    return evidence, stream_summary


def build_target_adjacency_vector(
    target: str,
    *,
    project_root: Path | str = PROJECT_DEFAULT,
    output_dir: Path | str | None = None,
    write_outputs: bool = True,
    write_extended_matrix: bool = True,
    top_n: int = TOP_N,
) -> TargetAdjacencyResult:
    """Calculate target-to-universe edge probabilities and optionally persist.

    Parameters
    ----------
    target
        Mouse gene symbol or UniProt accession. The resolved protein must not
        already be present in the active node universe.
    project_root
        Root of ``graphical_bayesian_inference``.
    output_dir
        Optional output directory. By default, results are written under
        ``results/path_finding/target_extensions/<resolved_symbol>``.
    write_outputs
        If false, return in-memory results without writing files.
    write_extended_matrix
        If true, append the target vector to the current square matrix.
    top_n
        Number of distinct in-universe kinase genes retained per target site.

    Returns
    -------
    TargetAdjacencyResult
        The canonical-order vector, optional extended matrix, phosphosite and
        kinase-prediction audits, metadata, and output location.
    """

    project = Path(project_root).resolve()
    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    current_matrix_path = (
        project
        / "results/edge_characterization"
        / "localization_kinase_predictor_string_hpa_omnipath"
        / "combined_adjacency_matrix.tsv"
    )
    uniprot_path = (
        project
        / "data/edge_characterization/kinase_predictor"
        / "phosphosite_database/raw/uniprot_mouse_reference_proteome.tsv.gz"
    )
    for path in (universe_path, current_matrix_path, uniprot_path):
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    if universe["symbol"].duplicated().any():
        raise ValueError("Node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    if not symbols:
        raise ValueError("The active node universe is empty")
    target_symbol, target_entry, target_mapping_method = resolve_target(
        target,
        uniprot_path,
    )
    universe_ci = {symbol.casefold(): symbol for symbol in symbols}
    if target_symbol.casefold() in universe_ci:
        existing = universe_ci[target_symbol.casefold()]
        raise ValueError(
            f"Target {target_symbol} is already node {existing!r} in the "
            "active universe; use the existing matrix column instead."
        )
    target_uniprot = clean(target_entry["Entry"])

    localization, localization_meta = target_localization_evidence(
        project,
        target_symbol,
        universe,
    )
    phosphosites = target_phosphosites(
        project,
        target_symbol,
        target_entry,
    )
    kinase, kinase_predictions, kinase_summary = target_kinase_evidence(
        project,
        universe,
        phosphosites,
        top_n,
    )
    string, string_summary = target_string_evidence(
        project,
        universe,
        target_symbol,
        target_uniprot,
    )
    hpa, hpa_meta = target_hpa_evidence(
        project,
        universe,
        target_symbol,
    )
    omnipath, omnipath_summary = target_omnipath_evidence(
        project,
        universe,
        target_symbol,
        target_uniprot,
    )

    n = len(universe)
    prior = np.full(n, 0.5, dtype=float)
    posterior_localization = apply_bayes_factor(
        prior,
        localization["mpkccd_bayes_factor"].to_numpy(dtype=float),
    )
    posterior_kinase = apply_bayes_factor(
        posterior_localization,
        kinase["kinasepredictor_bayes_factor"].to_numpy(dtype=float),
    )
    posterior_string = apply_bayes_factor(
        posterior_kinase,
        string["string_bayes_factor"].to_numpy(dtype=float),
    )
    posterior_hpa = apply_bayes_factor(
        posterior_string,
        hpa["hpa_bayes_factor"].to_numpy(dtype=float),
    )
    final_posterior = apply_bayes_factor(
        posterior_hpa,
        omnipath["omnipath_bayes_factor"].to_numpy(dtype=float),
    )

    vector = universe[
        ["symbol", "display_symbol", "name", "classes", "node_type"]
    ].copy()
    vector.insert(0, "universe_index", np.arange(n))
    vector.insert(1, "target_symbol", target_symbol)
    vector.insert(2, "target_uniprot", target_uniprot)
    vector["prior_edge_probability"] = prior
    vector = pd.concat(
        [
            vector.reset_index(drop=True),
            localization.reset_index(drop=True),
        ],
        axis=1,
    )
    vector["posterior_after_mpkccd_localization"] = posterior_localization
    vector = pd.concat(
        [vector, kinase.reset_index(drop=True)],
        axis=1,
    )
    vector["posterior_after_kinasepredictor"] = posterior_kinase
    vector = pd.concat(
        [vector, string.reset_index(drop=True)],
        axis=1,
    )
    vector["posterior_after_string"] = posterior_string
    vector = pd.concat(
        [vector, hpa.reset_index(drop=True)],
        axis=1,
    )
    vector["posterior_after_hpa"] = posterior_hpa
    vector = pd.concat(
        [vector, omnipath.reset_index(drop=True)],
        axis=1,
    )
    vector["final_edge_probability"] = final_posterior
    vector["any_non_neutral_edge_evidence"] = (
        localization["mpkccd_bayes_factor"].to_numpy(dtype=float)
        > 1.0 + 1e-12
    ) | (
        kinase["kinasepredictor_bayes_factor"].to_numpy(dtype=float)
        > 1.0 + 1e-12
    ) | (
        string["string_bayes_factor"].to_numpy(dtype=float)
        > 1.0 + 1e-12
    ) | (
        hpa["hpa_bayes_factor"].to_numpy(dtype=float)
        > 1.0 + 1e-12
    ) | (
        omnipath["omnipath_bayes_factor"].to_numpy(dtype=float)
        > 1.0 + 1e-12
    )
    vector["final_rank"] = (
        pd.Series(final_posterior)
        .rank(method="min", ascending=False)
        .astype(int)
        .to_numpy()
    )

    current = pd.read_csv(
        current_matrix_path,
        sep="\t",
        index_col=0,
    )
    if (
        current.index.astype(str).tolist() != symbols
        or current.columns.astype(str).tolist() != symbols
    ):
        raise ValueError("Current adjacency matrix order differs from universe")
    current_values = current.to_numpy(dtype=float)
    if not np.array_equal(current_values, current_values.T):
        raise ValueError("Current adjacency matrix is not symmetric")
    extended: pd.DataFrame | None = None
    if write_extended_matrix:
        extended_symbols = [*symbols, target_symbol]
        extended_values = np.zeros((n + 1, n + 1), dtype=float)
        extended_values[:n, :n] = current_values
        extended_values[:n, n] = final_posterior
        extended_values[n, :n] = final_posterior
        extended = pd.DataFrame(
            extended_values,
            index=extended_symbols,
            columns=extended_symbols,
        )
        extended.index.name = "symbol"

    baseline = np.isclose(
        final_posterior,
        0.5,
        atol=1e-12,
        rtol=0.0,
    )
    validation = {
        "vector_length_matches_universe": len(vector) == n,
        "vector_preserves_universe_order": (
            vector["symbol"].tolist() == symbols
        ),
        "final_probabilities_bounded_zero_one": bool(
            ((final_posterior >= 0.0) & (final_posterior <= 1.0)).all()
        ),
        "neutral_all_stream_pairs_equal_exact_baseline": bool(
            np.array_equal(
                final_posterior[
                    ~vector["any_non_neutral_edge_evidence"].to_numpy()
                ],
                np.full(
                    int(
                        (~vector["any_non_neutral_edge_evidence"]).sum()
                    ),
                    0.5,
                ),
            )
        ),
        "existing_matrix_block_unchanged": bool(
            extended is None
            or np.array_equal(
                extended.to_numpy(dtype=float)[:n, :n],
                current_values,
            )
        ),
        "extended_matrix_symmetric": bool(
            extended is None
            or np.array_equal(
                extended.to_numpy(dtype=float),
                extended.to_numpy(dtype=float).T,
            )
        ),
        "extended_matrix_zero_diagonal": bool(
            extended is None
            or np.array_equal(
                np.diag(extended.to_numpy(dtype=float)),
                np.zeros(n + 1),
            )
        ),
        "target_not_in_original_universe": target_symbol not in symbols,
    }
    if not all(validation.values()):
        raise AssertionError(f"Target-extension validation failed: {validation}")

    stream_counts = {
        "mpkccd_bf_gt1": int(
            (vector["mpkccd_bayes_factor"] > 1.0 + 1e-12).sum()
        ),
        "kinasepredictor_bf_gt1": int(
            (
                vector["kinasepredictor_bayes_factor"]
                > 1.0 + 1e-12
            ).sum()
        ),
        "string_bf_gt1": int(
            (vector["string_bayes_factor"] > 1.0 + 1e-12).sum()
        ),
        "hpa_bf_gt1": int(
            (vector["hpa_bayes_factor"] > 1.0 + 1e-12).sum()
        ),
        "omnipath_bf_gt1": int(
            (vector["omnipath_bayes_factor"] > 1.0 + 1e-12).sum()
        ),
    }
    summary: dict[str, Any] = {
        "target_input": target,
        "target_symbol": target_symbol,
        "target_uniprot": target_uniprot,
        "target_protein_name": clean(target_entry["Protein names"]),
        "target_sequence_length": len(normalized_sequence(target_entry["Sequence"])),
        "target_mapping_method": target_mapping_method,
        "original_node_count": n,
        "extended_node_count": n + 1 if write_extended_matrix else None,
        "target_to_universe_pair_count": n,
        "nonbaseline_target_edges": int((~baseline).sum()),
        "baseline_target_edges": int(baseline.sum()),
        "maximum_target_edge_probability": float(final_posterior.max()),
        "median_nonbaseline_target_edge_probability": (
            float(np.median(final_posterior[~baseline]))
            if (~baseline).any()
            else None
        ),
        "stream_non_neutral_counts": stream_counts,
        "mpkccd_target_metadata": localization_meta,
        "kinasepredictor_target_summary": kinase_summary,
        "string_target_summary": string_summary,
        "hpa_target_metadata": hpa_meta,
        "omnipath_target_summary": omnipath_summary,
        "evidence_order": [
            "mpkCCD localization",
            "observed-site KinasePredictor",
            "STRING v12",
            "HPA v25.1 primary localization",
            "OmniPath core undirected",
        ],
        "missing_evidence_rule": (
            "A missing or floor-valued stream receives BF=1 and leaves the "
            "target-node edge probability exactly unchanged."
        ),
        "matrix_extension_rule": (
            f"The current {n}x{n} matrix is copied without modification. The "
            "target vector is appended as one symmetric row and column; the "
            "new target diagonal is zero."
        ),
        "validation": validation,
        "input_sha256": {
            universe_path.name: sha256_file(universe_path),
            current_matrix_path.name: sha256_file(current_matrix_path),
            uniprot_path.name: sha256_file(uniprot_path),
        },
    }

    resolved_output: Path | None = None
    if write_outputs:
        resolved_output = (
            Path(output_dir).resolve()
            if output_dir is not None
            else (
                project
                / "results/path_finding/target_extensions"
                / safe_name(target_symbol)
            )
        )
        resolved_output.mkdir(parents=True, exist_ok=True)
        vector_path = resolved_output / "target_adjacency_vector.tsv"
        ranked_path = resolved_output / "target_adjacency_vector_ranked.tsv"
        vector.to_csv(
            vector_path,
            sep="\t",
            index=False,
            float_format="%.9g",
        )
        vector.sort_values(
            ["final_edge_probability", "symbol"],
            ascending=[False, True],
            kind="stable",
        ).to_csv(
            ranked_path,
            sep="\t",
            index=False,
            float_format="%.9g",
        )
        phosphosite_path = resolved_output / "target_phosphosites.tsv"
        phosphosites.to_csv(
            phosphosite_path,
            sep="\t",
            index=False,
        )
        prediction_path = resolved_output / "target_kinase_predictions.tsv"
        kinase_predictions.to_csv(
            prediction_path,
            sep="\t",
            index=False,
            float_format="%.9g",
        )
        if extended is not None:
            extended.to_csv(
                resolved_output / "extended_adjacency_matrix.tsv",
                sep="\t",
                float_format="%.9g",
            )
        summary_path = resolved_output / "analysis_summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2, default=str) + "\n",
            encoding="utf-8",
        )
        readme = f"""# Target adjacency extension: {target_symbol}

This folder adds mouse protein `{target_symbol}` (`{target_uniprot}`) to the
current undirected {n}-node graph without changing the existing matrix block.

## Output semantics

- `target_adjacency_vector.tsv` contains one row per original node in canonical
  universe order. `final_edge_probability` is the appended adjacency value.
- `target_adjacency_vector_ranked.tsv` contains the same complete records sorted
  by final probability for inspection.
- `extended_adjacency_matrix.tsv` is the symmetric {n + 1} by {n + 1} matrix with the
  target appended last. The target diagonal is zero.
- `target_phosphosites.tsv` records the UniProt-annotated and PKA-KO-observed
  phosphosites used for target-specific KinasePredictor scoring.
- `target_kinase_predictions.tsv` records the top {top_n} distinct eligible
  kinase genes per scorable site.
- `analysis_summary.json` records coverage, counts, rules, hashes, and QC.

## Evidence order

1. mpkCCD differential-fraction localization;
2. observed/annotated-site KinasePredictor;
3. STRING v12 functional association;
4. HPA v25.1 primary localization; and
5. OmniPath core interaction evidence collapsed to an unsigned, undirected pair.

Each target-node hypothesis starts at 0.5. Missing evidence is neutral. Existing
node-specific thresholds are reused; target-specific localization thresholds are
calculated against all covered nodes in the fixed {n}-node universe. The
existing {n} by {n} posterior matrix is never recomputed.

## Path-finding handoff

Path finding can use `extended_adjacency_matrix.tsv` directly. Because larger
values mean more probable edges, convert probabilities to a nonnegative path
cost before a shortest-path algorithm, for example `cost = -log(probability)`.
Applying an edge-probability cutoff or a maximum hop count is recommended to
prevent long paths composed mostly of 0.5 baseline edges.
"""
        (resolved_output / "README.md").write_text(readme, encoding="utf-8")
        output_files = [
            vector_path,
            ranked_path,
            phosphosite_path,
            prediction_path,
            summary_path,
            resolved_output / "README.md",
        ]
        if extended is not None:
            output_files.append(resolved_output / "extended_adjacency_matrix.tsv")
        summary["output_sha256"] = {
            path.name: sha256_file(path) for path in output_files
        }
        summary_path.write_text(
            json.dumps(summary, indent=2, default=str) + "\n",
            encoding="utf-8",
        )

    return TargetAdjacencyResult(
        target_symbol=target_symbol,
        target_uniprot=target_uniprot,
        vector=vector,
        extended_matrix=extended,
        phosphosites=phosphosites,
        kinase_predictions=kinase_predictions,
        summary=summary,
        output_directory=resolved_output,
    )


def main() -> int:
    args = parse_args()
    result = build_target_adjacency_vector(
        args.target,
        project_root=args.project_root,
        output_dir=args.output_dir,
        write_outputs=True,
        write_extended_matrix=not args.no_extended_matrix,
        top_n=args.top_n,
    )
    print(json.dumps(result.summary, indent=2, default=str))
    if result.output_directory is not None:
        print(f"Output directory: {result.output_directory}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
