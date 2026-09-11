#!/usr/bin/env python3
"""Infer undirected kinase-protein edges from observed phosphosites.

For every observed/annotated phosphosite in the current node universe:

1. Score its centered 13-mer with the official KinasePredictor v0.8 matrices.
2. Calculate the site's T_q as the 75th percentile of all applicable matrix
   scores, following the Deshpande et al. Kinase-Target Match protocol. If that
   percentile is nonpositive, use the positive-score q75 as a recorded fallback.
3. Collapse multiple predictor models that map to the same mouse kinase gene.
4. Retain the top 10 kinase genes that are protein-kinase nodes in the current
   node universe.
5. Treat every retained site prediction as a candidate evidence event for the
   unordered kinase-protein pair; weak scores at the 0.5 floor are neutral.
6. Sequentially update the existing localization edge posterior.

No prediction is interpreted as evidence against an edge. Pairs without a
retained prediction therefore remain exactly at their localization posterior.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd


TOP_N = 10
Q = 0.75
MINIMUM_SITE_LIKELIHOOD = 0.5
NO_EDGE_LIKELIHOOD = 0.5
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWYBJ")
AA_INDEX = {amino_acid: index for index, amino_acid in enumerate(AA_ORDER)}
KINASE_LOGOS_URL = "https://esbl.nhlbi.nih.gov/Databases/Kinase_Logos/"
KINASE_PREDICTOR_URL = (
    "https://esbl.nhlbi.nih.gov/Databases/Kinase_Logos/KinasePredictor.html"
)

# Matrix filenames that differ from the official metadata-table labels.
EXPLICIT_MOUSE_ALIASES = {
    "AMPKa1": "Prkaa1",
    "AMPKa2": "Prkaa2",
    "CDC2": "Cdk1",
    "CDC7": "Cdc7",
    "CDK2": "Cdk2",
    "CDK3": "Cdk3",
    "CDK4": "Cdk4",
    "CDK5": "Cdk5",
    "CDK6": "Cdk6",
    "CDK7": "Cdk7",
    "CDK8": "Cdk8",
    "CDK9": "Cdk9",
    "EGFR[T790M-L858R]": "Egfr",
}

MAPPING_COLUMNS = [
    "matrix_set",
    "predictor_label",
    "predictor_family",
    "human_gene_symbol",
    "metadata_mouse_gene_symbol",
    "mapped_mouse_node",
    "mapping_method",
    "in_node_universe",
    "node_has_kinase_class",
    "eligible_for_prediction",
]

SITE_AUDIT_COLUMNS = [
    "target_protein",
    "site",
    "position",
    "residue",
    "observed_sources",
    "centralized_sequence_13mer",
    "input_kinasepredictor_scorable",
    "scoring_status",
    "matrix_set",
    "raw_site_Tq_q75",
    "site_Tq_q75",
    "Tq_method",
    "applicable_matrix_count",
    "eligible_kinase_gene_count",
    "predictions_emitted",
    "self_predictions_excluded_from_edges",
]

PREDICTION_COLUMNS = [
    "target_protein",
    "target_site",
    "position",
    "residue",
    "observed_sources",
    "centralized_sequence_13mer",
    "kinase_node",
    "predictor_label",
    "predictor_family",
    "in_universe_rank",
    "global_model_rank",
    "raw_score",
    "raw_site_Tq_q75",
    "site_Tq_q75",
    "Tq_method",
    "z",
    "site_edge_likelihood",
    "site_bayes_factor",
    "used_for_undirected_edge",
    "undirected_node_a",
    "undirected_node_b",
]

SUPPORTED_EDGE_COLUMNS = [
    "node_a",
    "node_b",
    "predicted_kinase_nodes",
    "phosphosite_protein_nodes",
    "prediction_hit_count",
    "unique_observed_site_count",
    "supporting_site_examples",
    "observed_sources",
    "best_in_universe_rank",
    "best_global_model_rank",
    "best_raw_score",
    "maximum_site_edge_likelihood",
    "kinase_log_bayes_factor",
    "kinase_bayes_factor",
    "localization_data_observed",
    "localization_posterior",
    "kinase_only_posterior",
    "combined_posterior",
    "posterior_increase_over_localization",
    "kinase_kinase_edge",
]


def parse_args() -> argparse.Namespace:
    default_project = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=default_project,
        help=f"Project directory (default: {default_project})",
    )
    parser.add_argument("--top-n", type=int, default=TOP_N)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"nan", "none", "-"} else text


def symbol_candidates(value: object) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    candidates = []
    for token in re.split(r"[,;/\s]+", text):
        token = token.strip()
        if not token:
            continue
        candidates.append(token)
        without_variant = re.sub(r"\[.*$", "", token)
        if without_variant and without_variant != token:
            candidates.append(without_variant)
    return list(dict.fromkeys(candidates))


def load_matrix_set(folder: Path) -> tuple[np.ndarray, np.ndarray]:
    labels: list[str] = []
    matrices: list[np.ndarray] = []
    for path in sorted(folder.glob("*.csv"), key=lambda item: item.stem.casefold()):
        frame = pd.read_csv(path, index_col=0).reindex(AA_ORDER)
        if frame.shape != (22, 13) or frame.isna().any().any():
            raise ValueError(f"Unexpected KinasePredictor matrix: {path} {frame.shape}")
        labels.append(path.stem)
        matrices.append(frame.to_numpy(dtype=float))
    if not matrices:
        raise FileNotFoundError(f"No KinasePredictor matrices found in {folder}")
    return np.asarray(labels, dtype=object), np.stack(matrices)


def load_metadata(path: Path) -> pd.DataFrame:
    raw = path.read_bytes()
    text = None
    for encoding in ("utf-8", "windows-1252", "latin-1"):
        try:
            text = raw.decode(encoding)
            break
        except UnicodeDecodeError:
            continue
    if text is None:
        raise UnicodeError(f"Could not decode {path}")
    table = max(pd.read_html(StringIO(text)), key=len).copy()
    if table.shape[1] < 9:
        raise ValueError(f"Unexpected Kinase Logos metadata table shape: {table.shape}")
    table = table.iloc[:, :9]
    table.columns = [
        "predictor_family",
        "predictor_label",
        "sugiyama_id",
        "human_gene_symbol",
        "human_uniprot",
        "description",
        "anti_logo",
        "mouse_gene_symbol",
        "mouse_uniprot",
    ]
    for column in table.columns:
        table[column] = table[column].map(clean_text)
    return table


def map_models(
    matrix_set: str,
    labels: np.ndarray,
    metadata: pd.DataFrame,
    universe_ci: dict[str, str],
    kinase_nodes: set[str],
) -> pd.DataFrame:
    metadata_by_label: dict[str, dict[str, str]] = {}
    for row in metadata.to_dict(orient="records"):
        metadata_by_label.setdefault(str(row["predictor_label"]), row)

    rows: list[dict[str, object]] = []
    for matrix_index, label_value in enumerate(labels):
        label = str(label_value)
        meta = metadata_by_label.get(label, {})
        mapped_node = ""
        method = "unmapped"

        alias = EXPLICIT_MOUSE_ALIASES.get(label)
        if alias and alias.upper() in universe_ci:
            mapped_node = universe_ci[alias.upper()]
            method = "explicit matrix-label alias"

        if not mapped_node:
            for source_name, value in (
                ("official mouse symbol", meta.get("mouse_gene_symbol", "")),
                ("human symbol case-normalized", meta.get("human_gene_symbol", "")),
                ("predictor label case-normalized", label),
            ):
                for candidate in symbol_candidates(value):
                    match = universe_ci.get(candidate.upper())
                    if match:
                        mapped_node = match
                        method = source_name
                        break
                if mapped_node:
                    break

        in_universe = bool(mapped_node)
        is_kinase = mapped_node in kinase_nodes
        rows.append(
            {
                "matrix_index": matrix_index,
                "matrix_set": matrix_set,
                "predictor_label": label,
                "predictor_family": clean_text(meta.get("predictor_family", "")),
                "human_gene_symbol": clean_text(meta.get("human_gene_symbol", "")),
                "metadata_mouse_gene_symbol": clean_text(
                    meta.get("mouse_gene_symbol", "")
                ),
                "mapped_mouse_node": mapped_node,
                "mapping_method": method,
                "in_node_universe": in_universe,
                "node_has_kinase_class": is_kinase,
                "eligible_for_prediction": in_universe and is_kinase,
            }
        )
    return pd.DataFrame(rows)


def write_matrix(path: Path, symbols: list[str], matrix: np.ndarray) -> None:
    frame = pd.DataFrame(matrix, index=symbols, columns=symbols)
    frame.index.name = "symbol"
    frame.to_csv(path, sep="\t", float_format="%.9g")


def stable_expit(values: np.ndarray | float) -> np.ndarray | float:
    array = np.asarray(values, dtype=float)
    result = 1.0 / (1.0 + np.exp(-np.clip(array, -700.0, 700.0)))
    if np.isscalar(values):
        return float(result)
    return result


def score_site(
    sequence: str,
    labels: np.ndarray,
    matrices: np.ndarray,
    model_mapping: pd.DataFrame,
    top_n: int,
) -> dict[str, object]:
    residue_indices = np.asarray([AA_INDEX[amino_acid] for amino_acid in sequence])
    scores = matrices[:, residue_indices, np.arange(13)].sum(axis=1)
    raw_threshold = float(np.quantile(scores, Q))
    threshold = raw_threshold
    threshold_method = "raw_score_q75"
    if not math.isfinite(threshold):
        raise ValueError(f"Nonfinite site T_q for {sequence}: {threshold}")
    if threshold <= 0:
        positive_scores = scores[scores > 0]
        if positive_scores.size:
            threshold = float(np.quantile(positive_scores, Q))
            threshold_method = "positive_score_q75_fallback"
        else:
            # All scores imply no motif compatibility. A finite scale is used
            # only to keep the calculation defined; every hit remains neutral.
            threshold = 1.0
            threshold_method = "no_positive_scores_neutral"

    global_order = np.argsort(-scores, kind="stable")
    global_ranks = np.empty(len(scores), dtype=int)
    global_ranks[global_order] = np.arange(1, len(scores) + 1)

    eligible = model_mapping.loc[model_mapping["eligible_for_prediction"]].copy()
    eligible["raw_score"] = scores[eligible["matrix_index"].to_numpy(dtype=int)]
    eligible["global_model_rank"] = global_ranks[
        eligible["matrix_index"].to_numpy(dtype=int)
    ]
    eligible = eligible.sort_values(
        ["mapped_mouse_node", "raw_score", "global_model_rank", "predictor_label"],
        ascending=[True, False, True, True],
        kind="stable",
    )
    eligible = eligible.drop_duplicates("mapped_mouse_node", keep="first")
    eligible = eligible.sort_values(
        ["raw_score", "mapped_mouse_node"],
        ascending=[False, True],
        kind="stable",
    ).head(top_n)
    eligible["in_universe_rank"] = np.arange(1, len(eligible) + 1)

    predictions: list[dict[str, object]] = []
    for row in eligible.to_dict(orient="records"):
        raw_score = float(row["raw_score"])
        z = max(raw_score, 0.0) / threshold
        likelihood = max(
            MINIMUM_SITE_LIKELIHOOD,
            1.0 - math.exp(-0.5 * z * z),
        )
        predictions.append(
            {
                **row,
                "raw_score": raw_score,
                "z": z,
                "site_edge_likelihood": likelihood,
                "site_bayes_factor": likelihood / NO_EDGE_LIKELIHOOD,
            }
        )
    return {
        "raw_threshold": raw_threshold,
        "threshold": threshold,
        "threshold_method": threshold_method,
        "applicable_matrix_count": len(labels),
        "eligible_kinase_gene_count": int(
            model_mapping.loc[
                model_mapping["eligible_for_prediction"], "mapped_mouse_node"
            ].nunique()
        ),
        "predictions": predictions,
    }


def main() -> int:
    args = parse_args()
    if args.top_n <= 0:
        raise ValueError("--top-n must be positive")

    project = args.project_root.resolve()
    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    observed_path = (
        project
        / "data/edge_characterization/kinase_predictor/phosphosite_database"
        / "observed_phosphosites.tsv"
    )
    metadata_html_path = (
        project
        / "data/edge_characterization/kinase_predictor"
        / "kinase_logos_metadata.html"
    )
    metadata_tsv_path = (
        project
        / "data/edge_characterization/kinase_predictor"
        / "kinase_logos_metadata.tsv"
    )
    matrix_root = (
        project
        / "data/kinase_predictor/v0.8/official_package/output_matrices"
    )
    localization_matrix_path = (
        project
        / "results/edge_characterization/localization"
        / "localization_adjacency_matrix.tsv"
    )
    localization_profiles_path = (
        project
        / "data/edge_characterization/localization/processed"
        / "node_localization_profiles.tsv"
    )
    output_dir = (
        project
        / "results/edge_characterization/localization_kinase_predictor"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    required_paths = [
        universe_path,
        observed_path,
        metadata_html_path,
        localization_matrix_path,
        localization_profiles_path,
    ]
    for path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    universe = pd.read_csv(universe_path, sep="\t", dtype=str)
    if universe["symbol"].duplicated().any():
        raise ValueError("Node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}
    universe_ci = {symbol.upper(): symbol for symbol in symbols}
    if len(universe_ci) != len(symbols):
        raise ValueError("Node universe contains case-insensitive symbol duplicates")

    class_tokens = universe["classes"].fillna("").map(
        lambda value: set(str(value).split(";"))
    )
    kinase_mask = (
        universe["node_type"].eq("protein")
        & class_tokens.map(lambda values: "kinase" in values)
    )
    kinase_nodes = set(universe.loc[kinase_mask, "symbol"])

    metadata = load_metadata(metadata_html_path)
    metadata.to_csv(metadata_tsv_path, sep="\t", index=False)

    st_labels, st_matrices = load_matrix_set(
        matrix_root / "Ser-Thr_output_matrices"
    )
    tyr_labels, tyr_matrices = load_matrix_set(
        matrix_root / "Tyr_output_matrices"
    )
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
    model_mapping = pd.concat([st_mapping, tyr_mapping], ignore_index=True)
    model_mapping[MAPPING_COLUMNS].to_csv(
        output_dir / "kinase_predictor_node_mapping.tsv",
        sep="\t",
        index=False,
    )

    localization = pd.read_csv(
        localization_matrix_path,
        sep="\t",
        index_col=0,
    ).reindex(index=symbols, columns=symbols)
    if localization.isna().any().any():
        raise ValueError("Localization matrix could not be aligned to the node universe")
    localization_matrix = localization.to_numpy(dtype=float)
    if not np.allclose(localization_matrix, localization_matrix.T):
        raise ValueError("Input localization adjacency matrix is not symmetric")
    if not np.allclose(np.diag(localization_matrix), 0.0):
        raise ValueError("Input localization adjacency diagonal is not zero")

    localization_profiles = pd.read_csv(
        localization_profiles_path,
        sep="\t",
        dtype={"symbol": str},
    ).set_index("symbol")
    localization_observed_by_node = (
        localization_profiles["localization_observed"]
        .astype(str)
        .str.lower()
        .eq("true")
        .reindex(symbols)
        .fillna(False)
        .to_numpy(dtype=bool)
    )

    observed_sites = pd.read_csv(observed_path, sep="\t", dtype=str).fillna("")
    protein_nodes = set(
        universe.loc[universe["node_type"].eq("protein"), "symbol"]
    )
    if not set(observed_sites["symbol"]).issubset(protein_nodes):
        raise ValueError("Observed-site table contains a target outside protein nodes")

    site_audit_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    score_cache: dict[tuple[str, str], dict[str, object]] = {}

    for site in observed_sites.to_dict(orient="records"):
        target = str(site["symbol"])
        target_site = str(site["site"])
        position = int(site["position"])
        residue = str(site["residue"])
        sequence = re.sub(
            r"[^A-Za-z]",
            "",
            str(site["centralized_sequence_13mer"]),
        ).upper()
        input_scorable = str(site["kinasepredictor_scorable"]).lower() == "true"
        status = "scored"
        matrix_set = ""
        result: dict[str, object] | None = None

        if not input_scorable:
            status = "not_scored_input_flag"
        elif len(sequence) != 13:
            status = f"not_scored_sequence_length_{len(sequence)}"
        elif sequence[6] != residue:
            status = f"not_scored_center_mismatch_{sequence[6]}"
        elif residue not in {"S", "T", "Y"}:
            status = f"not_scored_center_residue_{residue}"
        else:
            unsupported = sorted(set(sequence).difference(AA_INDEX))
            if unsupported:
                status = "not_scored_unsupported_residues_" + ",".join(unsupported)
            else:
                if residue == "Y":
                    matrix_set = "tyrosine"
                    labels, matrices, mapping = (
                        tyr_labels,
                        tyr_matrices,
                        tyr_mapping,
                    )
                else:
                    matrix_set = "serine_threonine"
                    labels, matrices, mapping = (
                        st_labels,
                        st_matrices,
                        st_mapping,
                    )
                cache_key = (matrix_set, sequence)
                if cache_key not in score_cache:
                    score_cache[cache_key] = score_site(
                        sequence,
                        labels,
                        matrices,
                        mapping,
                        args.top_n,
                    )
                result = score_cache[cache_key]

        emitted = 0
        self_excluded = 0
        if result is not None:
            for prediction in result["predictions"]:
                kinase_node = str(prediction["mapped_mouse_node"])
                used = kinase_node != target
                if not used:
                    self_excluded += 1
                else:
                    emitted += 1
                if symbol_to_index[kinase_node] <= symbol_to_index[target]:
                    node_a, node_b = kinase_node, target
                else:
                    node_a, node_b = target, kinase_node
                prediction_rows.append(
                    {
                        "target_protein": target,
                        "target_site": target_site,
                        "position": position,
                        "residue": residue,
                        "observed_sources": site["sources"],
                        "centralized_sequence_13mer": sequence,
                        "kinase_node": kinase_node,
                        "predictor_label": prediction["predictor_label"],
                        "predictor_family": prediction["predictor_family"],
                        "in_universe_rank": prediction["in_universe_rank"],
                        "global_model_rank": prediction["global_model_rank"],
                        "raw_score": prediction["raw_score"],
                        "raw_site_Tq_q75": result["raw_threshold"],
                        "site_Tq_q75": result["threshold"],
                        "Tq_method": result["threshold_method"],
                        "z": prediction["z"],
                        "site_edge_likelihood": prediction[
                            "site_edge_likelihood"
                        ],
                        "site_bayes_factor": prediction["site_bayes_factor"],
                        "used_for_undirected_edge": bool_text(used),
                        "undirected_node_a": node_a if used else "",
                        "undirected_node_b": node_b if used else "",
                    }
                )

        site_audit_rows.append(
            {
                "target_protein": target,
                "site": target_site,
                "position": position,
                "residue": residue,
                "observed_sources": site["sources"],
                "centralized_sequence_13mer": sequence,
                "input_kinasepredictor_scorable": bool_text(input_scorable),
                "scoring_status": status,
                "matrix_set": matrix_set,
                "raw_site_Tq_q75": result["raw_threshold"] if result else "",
                "site_Tq_q75": result["threshold"] if result else "",
                "Tq_method": result["threshold_method"] if result else "",
                "applicable_matrix_count": (
                    result["applicable_matrix_count"] if result else 0
                ),
                "eligible_kinase_gene_count": (
                    result["eligible_kinase_gene_count"] if result else 0
                ),
                "predictions_emitted": emitted,
                "self_predictions_excluded_from_edges": self_excluded,
            }
        )

    site_audit = pd.DataFrame(site_audit_rows, columns=SITE_AUDIT_COLUMNS)
    predictions = pd.DataFrame(prediction_rows, columns=PREDICTION_COLUMNS)
    site_audit.to_csv(
        output_dir / "observed_phosphosite_scoring_audit.tsv",
        sep="\t",
        index=False,
        float_format="%.9g",
    )
    with gzip.open(
        output_dir / "observed_phosphosite_top10_predictions.tsv.gz",
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        predictions.to_csv(
            handle,
            sep="\t",
            index=False,
            float_format="%.9g",
        )

    used_predictions = predictions.loc[
        predictions["used_for_undirected_edge"].eq("true")
    ].copy()
    grouped_predictions: dict[tuple[str, str], list[dict[str, object]]] = (
        defaultdict(list)
    )
    for row in used_predictions.to_dict(orient="records"):
        grouped_predictions[
            (str(row["undirected_node_a"]), str(row["undirected_node_b"]))
        ].append(row)

    n = len(symbols)
    kinase_log_bf_matrix = np.zeros((n, n), dtype=float)
    hit_count_matrix = np.zeros((n, n), dtype=int)
    supported_edge_rows: list[dict[str, object]] = []

    for (node_a, node_b), rows in grouped_predictions.items():
        i = symbol_to_index[node_a]
        j = symbol_to_index[node_b]
        log_bf = float(
            sum(math.log(float(row["site_bayes_factor"])) for row in rows)
        )
        kinase_log_bf_matrix[i, j] = kinase_log_bf_matrix[j, i] = log_bf
        hit_count_matrix[i, j] = hit_count_matrix[j, i] = len(rows)
        localization_probability = float(localization_matrix[i, j])
        localization_log_odds = math.log(
            localization_probability / (1.0 - localization_probability)
        )
        kinase_only_probability = float(stable_expit(log_bf))
        combined_probability = float(
            stable_expit(localization_log_odds + log_bf)
        )
        source_values = sorted(
            {
                source
                for row in rows
                for source in str(row["observed_sources"]).split(";")
                if source
            }
        )
        site_keys = sorted(
            {
                f"{row['target_protein']}:{row['target_site']}"
                for row in rows
            }
        )
        kinase_values = sorted({str(row["kinase_node"]) for row in rows})
        target_values = sorted({str(row["target_protein"]) for row in rows})
        supported_edge_rows.append(
            {
                "node_a": node_a,
                "node_b": node_b,
                "predicted_kinase_nodes": ";".join(kinase_values),
                "phosphosite_protein_nodes": ";".join(target_values),
                "prediction_hit_count": len(rows),
                "unique_observed_site_count": len(site_keys),
                "supporting_site_examples": ";".join(site_keys[:20]),
                "observed_sources": ";".join(source_values),
                "best_in_universe_rank": min(
                    int(row["in_universe_rank"]) for row in rows
                ),
                "best_global_model_rank": min(
                    int(row["global_model_rank"]) for row in rows
                ),
                "best_raw_score": max(float(row["raw_score"]) for row in rows),
                "maximum_site_edge_likelihood": max(
                    float(row["site_edge_likelihood"]) for row in rows
                ),
                "kinase_log_bayes_factor": log_bf,
                "kinase_bayes_factor": math.exp(min(log_bf, 700.0)),
                "localization_data_observed": bool_text(
                    bool(
                        localization_observed_by_node[i]
                        and localization_observed_by_node[j]
                    )
                ),
                "localization_posterior": localization_probability,
                "kinase_only_posterior": kinase_only_probability,
                "combined_posterior": combined_probability,
                "posterior_increase_over_localization": (
                    combined_probability - localization_probability
                ),
                "kinase_kinase_edge": bool_text(
                    node_a in kinase_nodes and node_b in kinase_nodes
                ),
            }
        )

    supported_edges = pd.DataFrame(
        supported_edge_rows,
        columns=SUPPORTED_EDGE_COLUMNS,
    ).sort_values(
        ["combined_posterior", "prediction_hit_count", "node_a", "node_b"],
        ascending=[False, False, True, True],
        kind="stable",
    )

    off_diagonal_localization = localization_matrix.copy()
    np.fill_diagonal(off_diagonal_localization, 0.5)
    clipped_localization = np.clip(
        off_diagonal_localization,
        1e-15,
        1.0 - 1e-15,
    )
    localization_log_odds_matrix = np.log(
        clipped_localization / (1.0 - clipped_localization)
    )
    kinase_only_matrix = stable_expit(kinase_log_bf_matrix)
    combined_matrix = stable_expit(
        localization_log_odds_matrix + kinase_log_bf_matrix
    )
    np.fill_diagonal(kinase_only_matrix, 0.0)
    np.fill_diagonal(combined_matrix, 0.0)

    if not np.allclose(kinase_log_bf_matrix, kinase_log_bf_matrix.T):
        raise AssertionError("Kinase evidence matrix is not symmetric")
    if not np.array_equal(hit_count_matrix, hit_count_matrix.T):
        raise AssertionError("Prediction hit-count matrix is not symmetric")
    if not np.allclose(combined_matrix, combined_matrix.T):
        raise AssertionError("Combined adjacency matrix is not symmetric")
    if not np.allclose(np.diag(combined_matrix), 0.0):
        raise AssertionError("Combined adjacency diagonal is not zero")

    write_matrix(
        output_dir / "kinase_predictor_log_bayes_factor_matrix.tsv",
        symbols,
        kinase_log_bf_matrix,
    )
    write_matrix(
        output_dir / "kinase_predictor_hit_count_matrix.tsv",
        symbols,
        hit_count_matrix,
    )
    write_matrix(
        output_dir / "kinase_predictor_edge_probability_matrix.tsv",
        symbols,
        kinase_only_matrix,
    )
    write_matrix(
        output_dir / "combined_adjacency_matrix.tsv",
        symbols,
        combined_matrix,
    )
    supported_edges.to_csv(
        output_dir / "kinase_supported_undirected_edges.tsv",
        sep="\t",
        index=False,
        float_format="%.9g",
    )

    upper_i, upper_j = np.triu_indices(n, k=1)
    pair_localization_observed = (
        localization_observed_by_node[upper_i]
        & localization_observed_by_node[upper_j]
    )
    all_edges = pd.DataFrame(
        {
            "node_a": np.asarray(symbols)[upper_i],
            "node_b": np.asarray(symbols)[upper_j],
            "localization_data_observed": pair_localization_observed,
            "localization_posterior": localization_matrix[upper_i, upper_j],
            "kinase_prediction_evidence": (
                hit_count_matrix[upper_i, upper_j] > 0
            ),
            "kinase_prediction_hit_count": hit_count_matrix[upper_i, upper_j],
            "kinase_log_bayes_factor": kinase_log_bf_matrix[upper_i, upper_j],
            "kinase_only_posterior": kinase_only_matrix[upper_i, upper_j],
            "combined_posterior": combined_matrix[upper_i, upper_j],
            "posterior_increase_over_localization": (
                combined_matrix[upper_i, upper_j]
                - localization_matrix[upper_i, upper_j]
            ),
        }
    )
    with gzip.open(
        output_dir / "combined_undirected_edges.tsv.gz",
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        all_edges.to_csv(
            handle,
            sep="\t",
            index=False,
            float_format="%.9g",
        )
    all_edges.sort_values(
        [
            "combined_posterior",
            "kinase_prediction_hit_count",
            "node_a",
            "node_b",
        ],
        ascending=[False, False, True, True],
        kind="stable",
    ).head(1000).to_csv(
        output_dir / "top_combined_undirected_edges.tsv",
        sep="\t",
        index=False,
        float_format="%.9g",
    )

    scored_audit = site_audit["scoring_status"].eq("scored")
    supported_pair_mask = hit_count_matrix[upper_i, upper_j] > 0
    localization_upper = localization_matrix[upper_i, upper_j]
    combined_upper = combined_matrix[upper_i, upper_j]
    mapping_eligible = model_mapping["eligible_for_prediction"]
    summary = {
        "node_count": n,
        "protein_node_count": int(universe["node_type"].eq("protein").sum()),
        "molecule_node_count": int(universe["node_type"].eq("molecule").sum()),
        "nodes_with_kinase_class": len(kinase_nodes),
        "kinase_nodes_with_at_least_one_predictor_model": int(
            model_mapping.loc[
                mapping_eligible, "mapped_mouse_node"
            ].nunique()
        ),
        "kinase_predictor_matrix_count": len(model_mapping),
        "serine_threonine_matrix_count": len(st_labels),
        "tyrosine_matrix_count": len(tyr_labels),
        "eligible_in_universe_predictor_model_count": int(mapping_eligible.sum()),
        "observed_or_annotated_phosphosite_count": len(observed_sites),
        "observed_phosphosites_scored": int(scored_audit.sum()),
        "observed_phosphosites_not_scored": int((~scored_audit).sum()),
        "scored_serine_threonine_sites": int(
            (
                scored_audit
                & site_audit["matrix_set"].eq("serine_threonine")
            ).sum()
        ),
        "scored_tyrosine_sites": int(
            (scored_audit & site_audit["matrix_set"].eq("tyrosine")).sum()
        ),
        "sites_using_positive_score_Tq_fallback": int(
            site_audit["Tq_method"].eq("positive_score_q75_fallback").sum()
        ),
        "sites_with_no_positive_scores": int(
            site_audit["Tq_method"].eq("no_positive_scores_neutral").sum()
        ),
        "top_n_in_universe_kinase_genes_per_site": args.top_n,
        "site_predictions_recorded": len(predictions),
        "site_predictions_with_non_neutral_likelihood": int(
            (
                pd.to_numeric(predictions["site_edge_likelihood"])
                > MINIMUM_SITE_LIKELIHOOD + 1e-12
            ).sum()
        ),
        "site_predictions_at_neutral_likelihood_floor": int(
            (
                pd.to_numeric(predictions["site_edge_likelihood"])
                <= MINIMUM_SITE_LIKELIHOOD + 1e-12
            ).sum()
        ),
        "self_predictions_excluded_from_edges": int(
            predictions["used_for_undirected_edge"].eq("false").sum()
        ),
        "site_predictions_used_as_edge_evidence": len(used_predictions),
        "unique_undirected_edges_with_kinase_prediction_evidence": len(
            supported_edges
        ),
        "unique_undirected_edges_with_non_neutral_kinase_evidence": int(
            (
                pd.to_numeric(supported_edges["kinase_log_bayes_factor"])
                > 1e-12
            ).sum()
        ),
        "unique_undirected_edges_with_only_neutral_ranked_hits": int(
            (
                pd.to_numeric(supported_edges["kinase_log_bayes_factor"])
                <= 1e-12
            ).sum()
        ),
        "kinase_kinase_edges_with_prediction_evidence": int(
            supported_edges["kinase_kinase_edge"].eq("true").sum()
        ),
        "kinase_nonkinase_edges_with_prediction_evidence": int(
            supported_edges["kinase_kinase_edge"].eq("false").sum()
        ),
        "possible_undirected_pairs": len(all_edges),
        "pairs_strengthened_beyond_localization": int(
            (combined_upper > localization_upper + 1e-12).sum()
        ),
        "pairs_unchanged_from_localization": int(
            np.isclose(combined_upper, localization_upper, atol=1e-12).sum()
        ),
        "localization_only_maximum_probability": float(localization_upper.max()),
        "combined_minimum_off_diagonal_probability": float(combined_upper.min()),
        "combined_maximum_off_diagonal_probability": float(combined_upper.max()),
        "combined_median_supported_edge_probability": float(
            np.median(combined_upper[supported_pair_mask])
        ),
        "matrix_symmetric": bool(np.allclose(combined_matrix, combined_matrix.T)),
        "diagonal_all_zero": bool(
            np.allclose(np.diag(combined_matrix), 0.0)
        ),
        "edge_directionality": (
            "Undirected. Predictions retain kinase and phosphosite-protein "
            "roles for audit, but all evidence is accumulated on the unordered "
            "{node_a,node_b} pair."
        ),
        "prediction_scope": (
            "Top 10 distinct mapped kinase genes among kinase-class nodes in "
            f"the {n}-node universe; duplicate predictor labels mapping to one "
            "gene are collapsed by maximum score."
        ),
        "phosphosite_scope": (
            "Only observed_phosphosites.tsv; sequence-only S/T/Y candidates "
            "are excluded."
        ),
        "Tq_background_definition": (
            "For each observed phosphosite, the complete applicable official "
            "KinasePredictor matrix-score distribution (237 S/T or 98 Y "
            "models), before restriction to the node universe. The raw q75 is "
            "used when positive; if it is nonpositive, q75 of the positive "
            "scores is used and flagged in the site audit."
        ),
        "site_likelihood_formula": (
            "max(0.5, 1-exp(-0.5*(max(score,0)/Tq)^2)); "
            "Tq=raw q75 or recorded positive-score q75 fallback"
        ),
        "integration_formula": (
            "logit(P_combined)=logit(P_localization)+"
            "sum_site_hits(log(site_likelihood/0.5))"
        ),
        "absence_of_prediction_rule": (
            "No negative evidence is assigned; unsupported pairs retain the "
            "localization posterior exactly."
        ),
        "localization_matrix": str(localization_matrix_path),
        "observed_phosphosite_table": str(observed_path),
        "kinase_predictor_version": "0.8",
        "kinase_predictor_url": KINASE_PREDICTOR_URL,
        "kinase_metadata_url": KINASE_LOGOS_URL,
        "input_sha256": {
            "localization_adjacency_matrix.tsv": sha256_file(
                localization_matrix_path
            ),
            "observed_phosphosites.tsv": sha256_file(observed_path),
            "kinase_logos_metadata.html": sha256_file(metadata_html_path),
        },
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    readme = f"""OBSERVED-PHOSPHOSITE KINASE EDGE INTEGRATION

This directory contains an UNDIRECTED {n} x {n} adjacency model integrating
subcellular-localization evidence with KinasePredictor motif predictions.

SCOPE
- Only phosphosites in observed_phosphosites.tsv are scored.
- Sequence-only S/T/Y candidates are excluded.
- Only predictor models that map to kinase-class protein nodes in the current
  {n}-node universe can generate predictions.
- For every scored site, the top {args.top_n} distinct eligible kinase genes are
  retained. Multiple mutant/isoform predictor labels mapping to the same gene
  are collapsed to the strongest model for that site.
- Self-predictions are recorded in the site-prediction audit but do not create
  diagonal edges.

UNDIRECTED INTERPRETATION
The site-level audit retains "kinase_node" and "target_protein" so the motif
assignment can be checked. For the graph, however, each hit contributes to the
unordered {{kinase, protein}} pair. If both nodes are kinases, evidence in either
biochemical direction accumulates on the same symmetric edge.

BAYESIAN INTEGRATION
For each site, T_q is the 75th percentile of all applicable official
KinasePredictor scores (237 S/T models or 98 tyrosine models). If this raw
percentile is nonpositive, T_q is conservatively recalculated as the 75th
percentile of positive scores and that fallback is marked in the site audit.
Each retained prediction receives:

    L_site = max(0.5, 1 - exp(-0.5 * (max(score, 0) / T_q)^2))

The no-edge reference likelihood is 0.5. Each retained site assignment is a
separate candidate evidence event:

    logit(P_combined) = logit(P_localization)
                          + sum(log(L_site / 0.5))

No prediction is not evidence against an edge. An unsupported pair therefore
retains its localization probability exactly. A weak top-ranked assignment
whose L_site remains at the 0.5 floor is recorded for audit but is neutral.

FILES
- combined_adjacency_matrix.tsv: final symmetric {n} x {n} edge probabilities.
- combined_undirected_edges.tsv.gz: all {len(all_edges):,} unique node pairs.
- kinase_supported_undirected_edges.tsv: only pairs with motif evidence.
- observed_phosphosite_top10_predictions.tsv.gz: site-level prediction audit.
- observed_phosphosite_scoring_audit.tsv: every observed site's scoring status.
- kinase_predictor_node_mapping.tsv: predictor-model to node mapping.
- kinase_predictor_log_bayes_factor_matrix.tsv: accumulated motif evidence.
- kinase_predictor_hit_count_matrix.tsv: supporting prediction counts.
- kinase_predictor_edge_probability_matrix.tsv: motif-only edge posterior from
  a neutral 0.5 prior.
- top_combined_undirected_edges.tsv: 1,000 highest final probabilities.
- analysis_summary.json: parameters, provenance, validation, and counts.

INTERPRETATION LIMIT
KinasePredictor scores sequence-motif compatibility. These edges are predicted
physical/functional interaction hypotheses, not proof of direct phosphorylation
and not a statement about direction of information flow.

Official KinasePredictor:
{KINASE_PREDICTOR_URL}
"""
    (output_dir / "README.txt").write_text(readme, encoding="utf-8")

    summary["output_sha256"] = {
        filename: sha256_file(output_dir / filename)
        for filename in (
            "combined_adjacency_matrix.tsv",
            "kinase_supported_undirected_edges.tsv",
            "observed_phosphosite_top10_predictions.tsv.gz",
            "observed_phosphosite_scoring_audit.tsv",
            "kinase_predictor_node_mapping.tsv",
            "README.txt",
        )
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
