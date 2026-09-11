#!/usr/bin/env python3
"""Build an undirected protein-colocalization graph for AlphaFold screening.

Three localization evidence streams are scored independently:

1. mpkCCD five-fraction proteomics (the project's existing likelihood matrix),
2. Human Protein Atlas localization through mouse-human orthology (existing
   primary HPA likelihood matrix), and
3. mouse COMPARTMENTS knowledge-channel annotations.

COMPARTMENTS annotations with confidence scores 4 or 5 are collapsed to the
maximum score per GO cellular-component term for each mouse gene. The resulting
confidence-weighted profiles are L2-normalized and compared by cosine
similarity. Per-node q75 backgrounds, the positive-overlap fallback, the
complement-of-minimum likelihood kernel, and reciprocal symmetrization match
the existing HPA workflow.

For every protein pair:

    logit(P_colocalized) =
        logit(0.5)
        + log(BF_mpkCCD)
        + log(BF_HPA)
        + log(BF_COMPARTMENTS)

All graph edges are undirected. The 20 small-molecule nodes are retained in a
full-universe audit matrix but are excluded from AlphaFold candidate pairs.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


PROJECT_DEFAULT = Path(__file__).resolve().parents[2]
Q = 0.75
MINIMUM_LIKELIHOOD = 0.5
NO_COLOCALIZATION_LIKELIHOOD = 0.5
COMPARTMENTS_MIN_SCORE = 4
COMPARTMENTS_SENSITIVITY_MIN_SCORE = 3
EPSILON = 1e-12
COMPARTMENTS_COLUMNS = [
    "ensembl_protein_id",
    "gene_symbol",
    "go_id",
    "location_name",
    "evidence_source",
    "go_evidence_code",
    "confidence_score",
]
COMPARTMENTS_PAGE_URL = "https://compartments.jensenlab.org/Downloads"
COMPARTMENTS_DATA_URL = (
    "https://download.jensenlab.org/"
    "mouse_compartment_knowledge_filtered.tsv"
)
COMPARTMENTS_PAPER_URL = "https://pubmed.ncbi.nlm.nih.gov/24573882/"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_DEFAULT,
        help=f"Project directory (default: {PROJECT_DEFAULT})",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_matrix(
    path: Path,
    symbols: list[str],
    matrix: np.ndarray,
    *,
    float_format: str = "%.9g",
) -> None:
    frame = pd.DataFrame(matrix, index=symbols, columns=symbols)
    frame.index.name = "symbol"
    frame.to_csv(path, sep="\t", float_format=float_format)


def read_aligned_matrix(path: Path, symbols: list[str]) -> np.ndarray:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    if frame.index.tolist() != symbols or frame.columns.tolist() != symbols:
        raise ValueError(f"Node order mismatch in {path}")
    return frame.to_numpy(dtype=float)


def load_compartments(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=COMPARTMENTS_COLUMNS,
        dtype={
            "ensembl_protein_id": str,
            "gene_symbol": str,
            "go_id": str,
            "location_name": str,
            "evidence_source": str,
            "go_evidence_code": str,
        },
    )
    if frame.shape[1] != len(COMPARTMENTS_COLUMNS):
        raise ValueError(f"Unexpected COMPARTMENTS shape: {frame.shape}")
    frame["confidence_score"] = pd.to_numeric(
        frame["confidence_score"], errors="raise"
    ).astype(int)
    if not frame["confidence_score"].between(1, 5).all():
        raise ValueError("COMPARTMENTS confidence scores must be between 1 and 5")
    if frame[COMPARTMENTS_COLUMNS[:-1]].isna().any().any():
        raise ValueError("COMPARTMENTS contains missing identifier/text fields")
    return frame


def collapse_node_locations(
    rows: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse repeated sources and isoforms to one maximum score per GO term."""
    collapsed_rows: list[dict[str, object]] = []
    for go_id, group in rows.groupby("go_id", sort=True):
        maximum = int(group["confidence_score"].max())
        strongest = group.loc[group["confidence_score"].eq(maximum)].sort_values(
            ["location_name", "evidence_source", "go_evidence_code"],
            kind="stable",
        )
        collapsed_rows.append(
            {
                "go_id": str(go_id),
                "location_name": str(strongest.iloc[0]["location_name"]),
                "maximum_confidence_score": maximum,
                "evidence_sources": ";".join(
                    sorted(set(group["evidence_source"].astype(str)))
                ),
                "go_evidence_codes": ";".join(
                    sorted(set(group["go_evidence_code"].astype(str)))
                ),
                "ensembl_protein_ids": ";".join(
                    sorted(set(group["ensembl_protein_id"].astype(str)))
                ),
                "raw_evidence_rows": int(len(group)),
            }
        )
    return pd.DataFrame(collapsed_rows)


def build_compartments_profiles(
    universe: pd.DataFrame,
    compartments: pd.DataFrame,
    string_mapping: pd.DataFrame,
    *,
    minimum_score: int,
    include_audit: bool,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Map COMPARTMENTS to the node universe and construct weighted profiles."""
    eligible = compartments.loc[
        compartments["confidence_score"].ge(minimum_score)
    ].copy()
    by_symbol = {
        str(symbol): group.copy()
        for symbol, group in eligible.groupby("gene_symbol", sort=False)
    }
    by_ensembl = {
        str(identifier): group.copy()
        for identifier, group in eligible.groupby(
            "ensembl_protein_id", sort=False
        )
    }
    string_by_symbol = string_mapping.set_index("symbol", drop=False)

    mapped: dict[str, pd.DataFrame] = {}
    mapping_rows: list[dict[str, object]] = []
    collapsed_audit: list[pd.DataFrame] = []

    for row in universe.to_dict(orient="records"):
        symbol = str(row["symbol"])
        node_type = str(row["node_type"])
        selected_string_id = ""
        selected_ensembl = ""
        mapping_status = "excluded_nonprotein_node"
        matched = pd.DataFrame(columns=eligible.columns)

        if symbol in string_by_symbol.index:
            selected_string_id = str(
                string_by_symbol.loc[symbol].get("string_id", "")
            )
            selected_ensembl = selected_string_id.removeprefix("10090.")

        if node_type == "protein":
            if symbol in by_symbol:
                matched = by_symbol[symbol]
                protein_ids = matched["ensembl_protein_id"].nunique()
                mapping_status = (
                    "mapped_exact_symbol_single_protein_id"
                    if protein_ids == 1
                    else "mapped_exact_symbol_multiple_protein_ids"
                )
            elif selected_ensembl in by_ensembl:
                matched = by_ensembl[selected_ensembl]
                mapping_status = "mapped_selected_ensembl_protein_id"
            else:
                mapping_status = "unmapped_no_high_confidence_annotation"

        collapsed = (
            collapse_node_locations(matched)
            if not matched.empty
            else pd.DataFrame(
                columns=[
                    "go_id",
                    "location_name",
                    "maximum_confidence_score",
                    "evidence_sources",
                    "go_evidence_codes",
                    "ensembl_protein_ids",
                    "raw_evidence_rows",
                ]
            )
        )
        if not collapsed.empty:
            mapped[symbol] = collapsed
            if include_audit:
                node_audit = collapsed.copy()
                node_audit.insert(0, "symbol", symbol)
                collapsed_audit.append(node_audit)

        mapping_rows.append(
            {
                "symbol": symbol,
                "display_symbol": row.get("display_symbol", symbol),
                "node_type": node_type,
                "mapping_status": mapping_status,
                "selected_string_id": selected_string_id,
                "selected_ensembl_protein_id": selected_ensembl,
                "matched_ensembl_protein_ids": ";".join(
                    sorted(set(matched["ensembl_protein_id"].astype(str)))
                )
                if not matched.empty
                else "",
                "raw_matched_rows": int(len(matched)),
                "profile_location_count": int(len(collapsed)),
                "profile_go_ids": ";".join(collapsed["go_id"].astype(str))
                if not collapsed.empty
                else "",
                "profile_location_names": ";".join(
                    collapsed["location_name"].astype(str)
                )
                if not collapsed.empty
                else "",
                "minimum_compartments_score": minimum_score,
            }
        )

    go_ids = sorted(
        {
            str(go_id)
            for collapsed in mapped.values()
            for go_id in collapsed["go_id"]
        }
    )
    go_to_index = {go_id: index for index, go_id in enumerate(go_ids)}
    profiles = np.zeros((len(universe), len(go_ids)), dtype=float)
    symbol_to_index = {
        str(symbol): index
        for index, symbol in enumerate(universe["symbol"].astype(str))
    }
    for symbol, collapsed in mapped.items():
        row_index = symbol_to_index[symbol]
        for annotation in collapsed.to_dict(orient="records"):
            profiles[row_index, go_to_index[str(annotation["go_id"])]] = (
                float(annotation["maximum_confidence_score"]) / 5.0
            )

    mapping = pd.DataFrame(mapping_rows)
    evidence_audit = (
        pd.concat(collapsed_audit, ignore_index=True)
        if collapsed_audit
        else pd.DataFrame()
    )
    location_rows: list[dict[str, object]] = []
    for go_id in go_ids:
        column = profiles[:, go_to_index[go_id]]
        source_rows = eligible.loc[eligible["go_id"].eq(go_id)]
        location_rows.append(
            {
                "go_id": go_id,
                "location_name": str(
                    source_rows.sort_values(
                        ["confidence_score", "location_name"],
                        ascending=[False, True],
                        kind="stable",
                    ).iloc[0]["location_name"]
                ),
                "profiled_node_count": int((column > 0).sum()),
                "maximum_score_in_universe": float(column.max() * 5.0),
                "minimum_included_score": minimum_score,
            }
        )
    location_dictionary = pd.DataFrame(location_rows)
    return profiles, mapping, evidence_audit, location_dictionary


def score_profiles(weighted_profiles: np.ndarray) -> dict[str, object]:
    """Apply the HPA-style per-node q75 complement-of-minimum procedure."""
    observed_mask = weighted_profiles.sum(axis=1) > 0
    observed_indices = np.flatnonzero(observed_mask)
    observed = weighted_profiles[observed_mask].astype(float)
    norms = np.linalg.norm(observed, axis=1, keepdims=True)
    unit_profiles = observed / norms
    similarity_observed = unit_profiles @ unit_profiles.T
    similarity_observed = np.clip(similarity_observed, 0.0, 1.0)

    threshold_background = similarity_observed.copy()
    np.fill_diagonal(threshold_background, np.nan)
    raw_thresholds = np.nanquantile(threshold_background, Q, axis=1)
    thresholds = raw_thresholds.copy()
    methods: list[str] = []
    positive_background_sizes: list[int] = []

    for index, raw_threshold in enumerate(raw_thresholds):
        positive = threshold_background[index][
            np.isfinite(threshold_background[index])
            & (threshold_background[index] > 0)
        ]
        positive_background_sizes.append(int(len(positive)))
        if np.isfinite(raw_threshold) and raw_threshold > 0:
            methods.append("all_observed_nodes_q75")
        elif len(positive):
            thresholds[index] = float(np.quantile(positive, Q))
            methods.append("positive_overlap_q75_fallback")
        else:
            thresholds[index] = np.nan
            methods.append("no_positive_overlap_neutral")

    directed = np.full_like(similarity_observed, MINIMUM_LIKELIHOOD)
    for index, threshold in enumerate(thresholds):
        if np.isfinite(threshold) and threshold > 0:
            z = similarity_observed[index] / threshold
            directed[index] = np.maximum(
                MINIMUM_LIKELIHOOD,
                1.0 - np.exp(-0.5 * np.square(z)),
            )
    np.fill_diagonal(directed, MINIMUM_LIKELIHOOD)
    symmetric = (directed + directed.T) / 2.0

    n = len(weighted_profiles)
    similarity = np.zeros((n, n), dtype=float)
    likelihood = np.full((n, n), MINIMUM_LIKELIHOOD, dtype=float)
    factor_a = np.full((n, n), MINIMUM_LIKELIHOOD, dtype=float)
    similarity[np.ix_(observed_indices, observed_indices)] = (
        similarity_observed
    )
    likelihood[np.ix_(observed_indices, observed_indices)] = symmetric
    factor_a[np.ix_(observed_indices, observed_indices)] = directed
    np.fill_diagonal(similarity, 0.0)
    np.fill_diagonal(likelihood, 0.0)
    np.fill_diagonal(factor_a, 0.0)

    threshold_full = np.full(n, np.nan, dtype=float)
    raw_threshold_full = np.full(n, np.nan, dtype=float)
    method_full = np.full(n, "profile_not_observed", dtype=object)
    positive_size_full = np.zeros(n, dtype=int)
    threshold_full[observed_indices] = thresholds
    raw_threshold_full[observed_indices] = raw_thresholds
    method_full[observed_indices] = methods
    positive_size_full[observed_indices] = positive_background_sizes
    return {
        "observed_mask": observed_mask,
        "similarity": similarity,
        "likelihood": likelihood,
        "directed_factor": factor_a,
        "threshold": threshold_full,
        "raw_threshold": raw_threshold_full,
        "threshold_method": method_full,
        "positive_background_size": positive_size_full,
    }


def posterior_from_log_bf(log_bf: np.ndarray) -> np.ndarray:
    """Stable logistic transform for a 0.5 prior (zero prior log-odds)."""
    result = np.empty_like(log_bf, dtype=float)
    positive = log_bf >= 0
    result[positive] = 1.0 / (1.0 + np.exp(-log_bf[positive]))
    exponent = np.exp(log_bf[~positive])
    result[~positive] = exponent / (1.0 + exponent)
    return result


def write_gzip_table(path: Path, frame: pd.DataFrame) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        frame.to_csv(
            handle,
            sep="\t",
            index=False,
            float_format="%.9g",
        )


def plot_summary(
    output_path: Path,
    matrix: np.ndarray,
    posterior_upper: np.ndarray,
    supported: np.ndarray,
    excluded: np.ndarray,
) -> None:
    """Draw a dependency-light summary with Pillow.

    Keeping plotting out of the numerical stack makes the workflow runnable
    with the project's bundled pandas/numpy/Pillow environment.
    """

    width, height = 1800, 610
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()

    # Panel 1: matrix heatmap.  A compact purple-teal-yellow ramp provides
    # sufficient contrast around the neutral posterior of 0.5.
    heatmap = matrix.copy()
    diagonal = np.eye(len(heatmap), dtype=bool)
    vmax = max(0.5001, float(np.max(heatmap[~diagonal])))
    scaled = np.clip((heatmap - 0.5) / (vmax - 0.5), 0.0, 1.0)
    red = (48 + 207 * scaled**1.8).astype(np.uint8)
    green = (24 + 205 * scaled).astype(np.uint8)
    blue = (105 - 65 * scaled + 15 * scaled**2).astype(np.uint8)
    rgb = np.stack([red, green, blue], axis=2)
    rgb[diagonal] = 255
    heat_image = Image.fromarray(rgb, mode="RGB").resize(
        (500, 500), resample=Image.Resampling.NEAREST
    )
    canvas.paste(heat_image, (45, 65))
    draw.rectangle((45, 65, 545, 565), outline="#333333", width=1)
    draw.text((45, 35), "Integrated protein colocalization", fill="black", font=font)
    draw.text((225, 575), "Protein node index", fill="black", font=font)
    draw.text((45, 575), f"posterior: 0.5 to {vmax:.3f}", fill="#555555", font=font)

    # Panel 2: posterior histogram.
    hist_left, hist_top, hist_right, hist_bottom = 630, 80, 1210, 520
    hist_counts, hist_edges = np.histogram(
        posterior_upper, bins=80, range=(0.5, 1.0)
    )
    hist_max = max(1, int(hist_counts.max()))
    bar_width = (hist_right - hist_left) / len(hist_counts)
    for index, count in enumerate(hist_counts):
        x0 = hist_left + index * bar_width
        x1 = hist_left + (index + 1) * bar_width
        y0 = hist_bottom - (hist_bottom - hist_top) * count / hist_max
        draw.rectangle((x0, y0, x1, hist_bottom), fill="#0F6B78")
    draw.line(
        (hist_left, hist_top, hist_left, hist_bottom), fill="#333333", width=1
    )
    draw.line(
        (hist_left, hist_bottom, hist_right, hist_bottom),
        fill="#333333",
        width=1,
    )
    draw.text((630, 35), "All undirected protein pairs", fill="black", font=font)
    draw.text((630, 535), "0.5", fill="black", font=font)
    draw.text((1190, 535), "1.0", fill="black", font=font)
    draw.text((835, 560), "Colocalization posterior", fill="black", font=font)
    draw.text((630, 60), f"maximum bin count: {hist_max:,}", fill="#555555", font=font)

    # Panel 3: AlphaFold pair-screening counts.
    counts = [
        int(len(posterior_upper)),
        int(supported.sum()),
        int(excluded.sum()),
        int((~excluded).sum()),
    ]
    labels = ["All pairs", "BF-supported", "Exclude", "Retain"]
    colors = ["#A6A6A6", "#2E75B6", "#C00000", "#70AD47"]
    plot_left, plot_top, plot_right, plot_bottom = 1290, 100, 1760, 520
    count_max = max(counts)
    slot = (plot_right - plot_left) / len(counts)
    for index, (label, count, color) in enumerate(
        zip(labels, counts, colors, strict=True)
    ):
        x0 = plot_left + index * slot + 15
        x1 = plot_left + (index + 1) * slot - 15
        y0 = plot_bottom - (plot_bottom - plot_top) * count / count_max
        draw.rectangle((x0, y0, x1, plot_bottom), fill=color)
        draw.text((x0, max(plot_top - 15, y0 - 18)), f"{count:,}", fill="black", font=font)
        draw.text((x0, 535), label, fill="black", font=font)
    draw.line(
        (plot_left, plot_bottom, plot_right, plot_bottom),
        fill="#333333",
        width=1,
    )
    draw.text((1290, 35), "AlphaFold screening", fill="black", font=font)
    draw.text(
        (1290, 60),
        "Recommended filter requires explicit disjoint localization",
        fill="#555555",
        font=font,
    )

    canvas.save(output_path, format="PNG", optimize=True)


def main() -> int:
    args = parse_args()
    project = args.project_root.resolve()
    universe_path = (
        project / "data/node_selection/node_universe_combined_nonzero.tsv"
    )
    string_mapping_path = (
        project
        / "data/edge_characterization/string/v12.0/processed/"
        "node_to_string_mapping.tsv"
    )
    compartments_path = (
        project
        / "data/colocalization/compartments/raw/"
        "mouse_compartment_knowledge_filtered.tsv"
    )
    compartments_full_path = (
        project
        / "data/colocalization/compartments/raw/"
        "mouse_compartment_knowledge_full.tsv"
    )
    mpk_likelihood_path = (
        project
        / "results/edge_characterization/localization/"
        "localization_evidence_matrix.tsv"
    )
    mpk_profile_path = (
        project
        / "data/edge_characterization/localization/processed/"
        "node_localization_profiles.tsv"
    )
    hpa_likelihood_path = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa/"
        "hpa_localization_likelihood_matrix.tsv"
    )
    hpa_similarity_path = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa/"
        "hpa_localization_similarity_matrix.tsv"
    )
    hpa_observed_pair_path = (
        project
        / "results/edge_characterization/"
        "localization_kinase_predictor_string_hpa/"
        "hpa_localization_observed_pair_matrix.tsv"
    )
    processed_dir = (
        project / "data/colocalization/compartments/processed"
    )
    output_dir = project / "results/colocalization"
    processed_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    required = [
        universe_path,
        string_mapping_path,
        compartments_path,
        compartments_full_path,
        mpk_likelihood_path,
        mpk_profile_path,
        hpa_likelihood_path,
        hpa_similarity_path,
        hpa_observed_pair_path,
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Required input not found: {path}")

    universe = pd.read_csv(universe_path, sep="\t", dtype=str).fillna("")
    if universe["symbol"].duplicated().any():
        raise ValueError("Node universe contains duplicate symbols")
    symbols = universe["symbol"].astype(str).tolist()
    if len(symbols) != 868:
        raise ValueError(f"Expected 868 nodes, found {len(symbols)}")
    node_type = universe["node_type"].astype(str).to_numpy()
    protein_indices = np.flatnonzero(node_type == "protein")
    protein_symbols = np.asarray(symbols)[protein_indices].tolist()
    if len(protein_symbols) != 848:
        raise ValueError(
            f"Expected 848 protein nodes, found {len(protein_symbols)}"
        )

    compartments = load_compartments(compartments_path)
    string_mapping = pd.read_csv(
        string_mapping_path, sep="\t", dtype=str
    ).fillna("")
    profiles, mapping, evidence_audit, location_dictionary = (
        build_compartments_profiles(
            universe,
            compartments,
            string_mapping,
            minimum_score=COMPARTMENTS_MIN_SCORE,
            include_audit=True,
        )
    )
    primary_score = score_profiles(profiles)
    sensitivity_profiles, sensitivity_mapping, _, _ = (
        build_compartments_profiles(
            universe,
            compartments,
            string_mapping,
            minimum_score=COMPARTMENTS_SENSITIVITY_MIN_SCORE,
            include_audit=False,
        )
    )
    sensitivity_score = score_profiles(sensitivity_profiles)

    mapping["raw_T_q75"] = primary_score["raw_threshold"]
    mapping["T_q75"] = primary_score["threshold"]
    mapping["T_q_method"] = primary_score["threshold_method"]
    mapping["positive_overlap_background_size"] = primary_score[
        "positive_background_size"
    ]
    mapping["profile_observed"] = primary_score["observed_mask"]
    mapping["score3plus_profile_observed"] = sensitivity_score[
        "observed_mask"
    ]
    mapping_path = processed_dir / "node_compartments_profiles.tsv"
    evidence_path = processed_dir / "mapped_node_location_evidence.tsv"
    location_path = processed_dir / "compartments_location_dictionary.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)
    evidence_audit.to_csv(evidence_path, sep="\t", index=False)
    location_dictionary.to_csv(location_path, sep="\t", index=False)

    mpk_likelihood = read_aligned_matrix(
        mpk_likelihood_path, symbols
    )
    hpa_likelihood = read_aligned_matrix(
        hpa_likelihood_path, symbols
    )
    hpa_similarity = read_aligned_matrix(
        hpa_similarity_path, symbols
    )
    hpa_observed_pair = read_aligned_matrix(
        hpa_observed_pair_path, symbols
    ).astype(bool)
    compartments_similarity = np.asarray(
        primary_score["similarity"], dtype=float
    )
    compartments_likelihood = np.asarray(
        primary_score["likelihood"], dtype=float
    )
    compartments_observed = np.asarray(
        primary_score["observed_mask"], dtype=bool
    )
    compartments_observed_pair = np.outer(
        compartments_observed, compartments_observed
    )
    np.fill_diagonal(compartments_observed_pair, False)

    mpk_profiles = pd.read_csv(
        mpk_profile_path, sep="\t", dtype=str
    ).fillna("")
    if mpk_profiles["symbol"].astype(str).tolist() != symbols:
        raise ValueError("mpkCCD profile order does not match node universe")
    mpk_observed = (
        mpk_profiles["localization_observed"]
        .astype(str)
        .str.casefold()
        .eq("true")
        .to_numpy()
    )
    mpk_observed_pair = np.outer(mpk_observed, mpk_observed)
    np.fill_diagonal(mpk_observed_pair, False)

    n = len(symbols)
    upper_i, upper_j = np.triu_indices(n, k=1)
    mpk_upper = mpk_likelihood[upper_i, upper_j]
    hpa_upper = hpa_likelihood[upper_i, upper_j]
    compartments_upper = compartments_likelihood[upper_i, upper_j]
    log_bf_mpk = np.log(mpk_upper / NO_COLOCALIZATION_LIKELIHOOD)
    log_bf_hpa = np.log(hpa_upper / NO_COLOCALIZATION_LIKELIHOOD)
    log_bf_compartments = np.log(
        compartments_upper / NO_COLOCALIZATION_LIKELIHOOD
    )
    combined_log_bf_upper = (
        log_bf_mpk + log_bf_hpa + log_bf_compartments
    )
    posterior_upper = posterior_from_log_bf(combined_log_bf_upper)

    full_matrix = np.zeros((n, n), dtype=float)
    full_matrix[upper_i, upper_j] = posterior_upper
    full_matrix[upper_j, upper_i] = posterior_upper
    combined_log_bf_matrix = np.zeros((n, n), dtype=float)
    combined_log_bf_matrix[upper_i, upper_j] = combined_log_bf_upper
    combined_log_bf_matrix[upper_j, upper_i] = combined_log_bf_upper

    protein_matrix = full_matrix[np.ix_(protein_indices, protein_indices)]
    protein_log_bf = combined_log_bf_matrix[
        np.ix_(protein_indices, protein_indices)
    ]
    write_matrix(
        output_dir / "colocalization_adjacency_matrix.tsv",
        protein_symbols,
        protein_matrix,
    )
    write_matrix(
        output_dir / "colocalization_log_bayes_factor_matrix.tsv",
        protein_symbols,
        protein_log_bf,
    )
    write_matrix(
        output_dir / "colocalization_adjacency_matrix_all_nodes.tsv",
        symbols,
        full_matrix,
    )
    write_matrix(
        output_dir / "compartments_similarity_matrix.tsv",
        symbols,
        compartments_similarity,
    )
    write_matrix(
        output_dir / "compartments_likelihood_matrix.tsv",
        symbols,
        compartments_likelihood,
    )
    write_matrix(
        output_dir / "compartments_log_bayes_factor_matrix.tsv",
        symbols,
        np.where(
            np.eye(n, dtype=bool),
            0.0,
            np.log(
                np.maximum(
                    compartments_likelihood, MINIMUM_LIKELIHOOD
                )
                / NO_COLOCALIZATION_LIKELIHOOD
            ),
        ),
    )
    write_matrix(
        output_dir / "compartments_observed_pair_matrix.tsv",
        symbols,
        compartments_observed_pair.astype(np.uint8),
        float_format="%d",
    )
    write_matrix(
        output_dir / "compartments_score3plus_likelihood_matrix.tsv",
        symbols,
        np.asarray(sensitivity_score["likelihood"], dtype=float),
    )

    protein_upper_i, protein_upper_j = np.triu_indices(
        len(protein_symbols), k=1
    )
    full_protein_i = protein_indices[protein_upper_i]
    full_protein_j = protein_indices[protein_upper_j]
    protein_posterior_upper = protein_matrix[
        protein_upper_i, protein_upper_j
    ]
    protein_log_bf_upper = protein_log_bf[
        protein_upper_i, protein_upper_j
    ]
    protein_mpk_likelihood = mpk_likelihood[
        full_protein_i, full_protein_j
    ]
    protein_hpa_likelihood = hpa_likelihood[
        full_protein_i, full_protein_j
    ]
    protein_comp_likelihood = compartments_likelihood[
        full_protein_i, full_protein_j
    ]
    protein_hpa_similarity = hpa_similarity[
        full_protein_i, full_protein_j
    ]
    protein_comp_similarity = compartments_similarity[
        full_protein_i, full_protein_j
    ]
    protein_mpk_covered = mpk_observed_pair[
        full_protein_i, full_protein_j
    ]
    protein_hpa_covered = hpa_observed_pair[
        full_protein_i, full_protein_j
    ]
    protein_comp_covered = compartments_observed_pair[
        full_protein_i, full_protein_j
    ]

    supported = protein_posterior_upper > 0.5 + EPSILON
    hpa_disjoint = protein_hpa_covered & (
        protein_hpa_similarity <= EPSILON
    )
    compartments_disjoint = protein_comp_covered & (
        protein_comp_similarity <= EPSILON
    )
    direct_disjoint_count = (
        hpa_disjoint.astype(int) + compartments_disjoint.astype(int)
    )
    any_pair_coverage = (
        protein_mpk_covered
        | protein_hpa_covered
        | protein_comp_covered
    )
    explicitly_shared = (
        (protein_hpa_covered & (protein_hpa_similarity > EPSILON))
        | (
            protein_comp_covered
            & (protein_comp_similarity > EPSILON)
        )
    )
    recommended_exclude = (
        ~supported
        & ~explicitly_shared
        & (hpa_disjoint | compartments_disjoint)
    )
    recommended_retain = ~recommended_exclude
    multi_source_disjoint = hpa_disjoint & compartments_disjoint

    protein_symbol_array = np.asarray(protein_symbols)
    edge_frame = pd.DataFrame(
        {
            "edge_id": [
                f"{protein_symbols[i]}|{protein_symbols[j]}"
                for i, j in zip(
                    protein_upper_i,
                    protein_upper_j,
                    strict=True,
                )
            ],
            "node_a": protein_symbol_array[protein_upper_i],
            "node_b": protein_symbol_array[protein_upper_j],
            "mpkccd_covered": protein_mpk_covered,
            "mpkccd_likelihood": protein_mpk_likelihood,
            "mpkccd_bayes_factor": (
                protein_mpk_likelihood
                / NO_COLOCALIZATION_LIKELIHOOD
            ),
            "hpa_covered": protein_hpa_covered,
            "hpa_cosine_similarity": np.where(
                protein_hpa_covered,
                protein_hpa_similarity,
                np.nan,
            ),
            "hpa_likelihood": protein_hpa_likelihood,
            "hpa_bayes_factor": (
                protein_hpa_likelihood
                / NO_COLOCALIZATION_LIKELIHOOD
            ),
            "compartments_covered": protein_comp_covered,
            "compartments_cosine_similarity": np.where(
                protein_comp_covered,
                protein_comp_similarity,
                np.nan,
            ),
            "compartments_likelihood": protein_comp_likelihood,
            "compartments_bayes_factor": (
                protein_comp_likelihood
                / NO_COLOCALIZATION_LIKELIHOOD
            ),
            "combined_log_bayes_factor": protein_log_bf_upper,
            "colocalization_posterior": protein_posterior_upper,
            "supported_colocalization_edge": supported,
            "any_pair_coverage": any_pair_coverage,
            "explicitly_shared_discrete_compartment": explicitly_shared,
            "hpa_profiles_disjoint": hpa_disjoint,
            "compartments_profiles_disjoint": compartments_disjoint,
            "direct_disjoint_source_count": direct_disjoint_count,
            "multi_source_disjoint": multi_source_disjoint,
            "recommended_exclude_from_alphafold": recommended_exclude,
            "recommended_retain_for_alphafold": recommended_retain,
        }
    )
    all_edges_path = output_dir / "colocalization_protein_pairs.tsv.gz"
    supported_path = output_dir / "supported_colocalization_edges.tsv.gz"
    candidate_path = (
        output_dir / "alphafold_candidate_pairs_recommended.tsv.gz"
    )
    excluded_path = (
        output_dir / "alphafold_pairs_excluded_by_localization.tsv.gz"
    )
    write_gzip_table(all_edges_path, edge_frame)
    write_gzip_table(supported_path, edge_frame.loc[supported])
    write_gzip_table(candidate_path, edge_frame.loc[recommended_retain])
    write_gzip_table(excluded_path, edge_frame.loc[recommended_exclude])
    edge_frame.sort_values(
        [
            "colocalization_posterior",
            "combined_log_bayes_factor",
            "node_a",
            "node_b",
        ],
        ascending=[False, False, True, True],
        kind="stable",
    ).head(2000).to_csv(
        output_dir / "top_colocalization_edges.tsv",
        sep="\t",
        index=False,
        float_format="%.9g",
    )

    plot_path = output_dir / "colocalization_graph_summary.png"
    plot_summary(
        plot_path,
        protein_matrix,
        protein_posterior_upper,
        supported,
        recommended_exclude,
    )

    sequential = np.full_like(protein_posterior_upper, 0.5)
    for likelihood in (
        protein_mpk_likelihood,
        protein_hpa_likelihood,
        protein_comp_likelihood,
    ):
        numerator = sequential * likelihood
        sequential = numerator / (
            numerator
            + (1.0 - sequential) * NO_COLOCALIZATION_LIKELIHOOD
        )
    reconstruction_error = float(
        np.max(np.abs(sequential - protein_posterior_upper))
    )
    protein_diagonal = np.diag(protein_matrix)
    primary_observed = np.asarray(
        primary_score["observed_mask"], dtype=bool
    )
    sensitivity_observed = np.asarray(
        sensitivity_score["observed_mask"], dtype=bool
    )
    primary_likelihood_upper = compartments_likelihood[
        full_protein_i, full_protein_j
    ]
    sensitivity_likelihood = np.asarray(
        sensitivity_score["likelihood"], dtype=float
    )
    sensitivity_likelihood_upper = sensitivity_likelihood[
        full_protein_i, full_protein_j
    ]
    comparable = (
        primary_observed[full_protein_i]
        & primary_observed[full_protein_j]
    )
    sensitivity_correlation = (
        float(
            np.corrcoef(
                primary_likelihood_upper[comparable],
                sensitivity_likelihood_upper[comparable],
            )[0, 1]
        )
        if comparable.sum() > 1
        else float("nan")
    )

    validation = {
        "full_matrix_shape": list(full_matrix.shape),
        "protein_matrix_shape": list(protein_matrix.shape),
        "full_matrix_symmetric": bool(
            np.array_equal(full_matrix, full_matrix.T)
        ),
        "protein_matrix_symmetric": bool(
            np.array_equal(protein_matrix, protein_matrix.T)
        ),
        "diagonal_all_zero": bool(
            np.array_equal(
                protein_diagonal,
                np.zeros(len(protein_symbols)),
            )
        ),
        "all_off_diagonal_probabilities_in_half_to_one": bool(
            np.all(
                (protein_posterior_upper >= 0.5)
                & (protein_posterior_upper <= 1.0)
            )
        ),
        "sequential_update_matches_log_bf_sum": bool(
            reconstruction_error <= 1e-12
        ),
        "candidate_and_excluded_partition_all_pairs": bool(
            recommended_retain.sum() + recommended_exclude.sum()
            == len(edge_frame)
        ),
        "excluded_pairs_have_no_bf_support": bool(
            not np.any(supported & recommended_exclude)
        ),
        "excluded_pairs_have_explicit_disjoint_evidence": bool(
            np.all(
                (hpa_disjoint | compartments_disjoint)[
                    recommended_exclude
                ]
            )
        ),
        "edge_ids_unique": bool(not edge_frame["edge_id"].duplicated().any()),
    }
    if not all(validation.values()):
        raise AssertionError(f"Validation failed: {validation}")

    mapping_counts = {
        str(key): int(value)
        for key, value in mapping["mapping_status"].value_counts().items()
    }
    tq_counts = {
        str(key): int(value)
        for key, value in mapping["T_q_method"].value_counts().items()
    }
    stream_supported = {
        "mpkccd": int(
            (protein_mpk_likelihood > MINIMUM_LIKELIHOOD + EPSILON).sum()
        ),
        "hpa": int(
            (protein_hpa_likelihood > MINIMUM_LIKELIHOOD + EPSILON).sum()
        ),
        "compartments": int(
            (
                protein_comp_likelihood
                > MINIMUM_LIKELIHOOD + EPSILON
            ).sum()
        ),
    }
    summary = {
        "purpose": (
            "Undirected colocalization graph and AlphaFold pair-screening "
            "table for the current signaling-node universe."
        ),
        "node_count_all": len(symbols),
        "protein_node_count": len(protein_symbols),
        "excluded_nonprotein_node_count": int(
            (node_type != "protein").sum()
        ),
        "possible_undirected_protein_pairs": int(len(edge_frame)),
        "supported_colocalization_edges": int(supported.sum()),
        "neutral_colocalization_pairs": int((~supported).sum()),
        "recommended_alphafold_pairs_retained": int(
            recommended_retain.sum()
        ),
        "recommended_alphafold_pairs_excluded": int(
            recommended_exclude.sum()
        ),
        "recommended_fraction_retained": float(
            recommended_retain.mean()
        ),
        "recommended_fraction_excluded": float(
            recommended_exclude.mean()
        ),
        "supported_only_fraction_retained": float(supported.mean()),
        "pairs_with_any_joint_source_coverage": int(
            any_pair_coverage.sum()
        ),
        "pairs_with_explicit_shared_discrete_location": int(
            explicitly_shared.sum()
        ),
        "pairs_hpa_profiles_disjoint": int(hpa_disjoint.sum()),
        "pairs_compartments_profiles_disjoint": int(
            compartments_disjoint.sum()
        ),
        "pairs_disjoint_in_both_discrete_sources": int(
            multi_source_disjoint.sum()
        ),
        "stream_non_neutral_edge_counts": stream_supported,
        "q": Q,
        "minimum_likelihood": MINIMUM_LIKELIHOOD,
        "edge_prior": 0.5,
        "no_colocalization_likelihood": (
            NO_COLOCALIZATION_LIKELIHOOD
        ),
        "combination_rule": (
            "logit(P_colocalized)=logit(0.5)+log(BF_mpkCCD)+"
            "log(BF_HPA)+log(BF_COMPARTMENTS)"
        ),
        "factor_rule": (
            "BF_source=undirected complement-of-minimum likelihood/0.5"
        ),
        "compartments_source_channel": "mouse knowledge channel, filtered",
        "compartments_minimum_confidence_score": (
            COMPARTMENTS_MIN_SCORE
        ),
        "compartments_raw_rows": int(len(compartments)),
        "compartments_profile_go_term_count": int(profiles.shape[1]),
        "protein_nodes_with_compartments_profile": int(
            primary_observed[protein_indices].sum()
        ),
        "compartments_mapping_status_counts": mapping_counts,
        "compartments_Tq_method_counts": tq_counts,
        "compartments_profile_definition": (
            "For each gene, maximum confidence score per GO cellular-"
            "component term across duplicate evidence and Ensembl protein "
            "isoforms; scores 4-5 included as weights score/5, followed by "
            "L2 normalization."
        ),
        "compartments_background_definition": (
            "For each profiled node, cosine similarities to every other "
            "profiled node in the current 868-node universe, self excluded; "
            "q75 of positive overlaps is used when the all-node q75 is zero."
        ),
        "compartments_directed_likelihood_formula": (
            "max(0.5, 1-exp(-0.5*(cosine_similarity/Tq_node)^2))"
        ),
        "symmetrization": (
            "Arithmetic mean of reciprocal endpoint-specific likelihoods."
        ),
        "recommended_filter_definition": (
            "Exclude a protein pair only when the integrated posterior is "
            "exactly neutral, neither HPA nor COMPARTMENTS records a shared "
            "discrete location, and at least one of those two sources has "
            "profiles for both endpoints with zero overlap. Unknown pairs "
            "and conflicts with positive evidence are retained."
        ),
        "sensitivity_score3plus": {
            "minimum_confidence_score": (
                COMPARTMENTS_SENSITIVITY_MIN_SCORE
            ),
            "profiled_protein_nodes": int(
                sensitivity_observed[protein_indices].sum()
            ),
            "likelihood_correlation_on_primary_covered_pairs": (
                sensitivity_correlation
            ),
        },
        "dependence_caveat": (
            "COMPARTMENTS knowledge annotations draw on GO, MGI, UniProt, "
            "Reactome and related sources. HPA reliability also references "
            "external annotations, so strict conditional independence is "
            "not guaranteed; the product-BF graph is a screening model."
        ),
        "interpretation_caveat": (
            "The complement-of-minimum model supplies positive or neutral "
            "evidence only. A 0.5 posterior is not proof of separation. The "
            "recommended exclusion list therefore additionally requires "
            "explicit zero overlap between observed discrete profiles."
        ),
        "validation": validation,
        "reconstruction_maximum_absolute_error": reconstruction_error,
        "source_urls": {
            "compartments_downloads": COMPARTMENTS_PAGE_URL,
            "compartments_data": COMPARTMENTS_DATA_URL,
            "compartments_paper": COMPARTMENTS_PAPER_URL,
            "mpkccd_localization": (
                "https://esbl.nhlbi.nih.gov/Databases/mpkFractions/"
                "proteomic_fractions_log_files/"
                "Proteomics_of_subcellular_fractions.xlsx"
            ),
            "hpa_subcellular": (
                "https://www.proteinatlas.org/humanproteome/subcellular"
            ),
        },
        "input_sha256": {
            path.name: sha256_file(path)
            for path in required
        },
    }
    summary_path = output_dir / "analysis_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    readme = f"""COLOCALIZATION GRAPH FOR ALPHAFOLD PAIR SCREENING

PURPOSE
This directory contains an UNDIRECTED protein-colocalization graph for the
current signaling-node universe. It is separate from the main signaling-edge
model. The graph is intended to reduce structure-prediction jobs by screening
protein pairs that have explicit, nonoverlapping localization profiles and no
positive colocalization evidence.

NODES
- Full project universe: {len(symbols):,} nodes.
- Protein nodes eligible for AlphaFold: {len(protein_symbols):,}.
- Small-molecule nodes excluded from AlphaFold: {int((node_type != "protein").sum()):,}.
- Possible unordered protein pairs: {len(edge_frame):,}.

EVIDENCE STREAMS
1. Mouse mpkCCD differential-centrifugation proteomics.
2. Human Protein Atlas primary localization mapped through stringent MGI
   one-to-one mouse-human orthology.
3. Mouse COMPARTMENTS knowledge-channel annotations scored 4 or 5 stars.

COMPARTMENTS PROFILE
The filtered knowledge file is used because the full ontology-propagated file
contains generic roots such as cellular_component that cannot discriminate
colocalization. Repeated evidence and Ensembl protein isoforms are collapsed to
the maximum confidence score for each gene/GO cellular-component term. The
profile weight is confidence/5 and vectors are L2-normalized.

For every profiled node, T_q is the 75th percentile of cosine similarities to
all other profiled nodes, excluding itself. If q75 is zero, the q75 of positive
overlaps is used and flagged.

    L_i = max(0.5, 1 - exp(-0.5 * (similarity / T_q_i)^2))
    L_undirected = (L_a + L_b) / 2
    BF_source = L_undirected / 0.5

INTEGRATION
Starting from P(colocalized)=0.5:

    logit(P_final) = log(BF_mpkCCD) + log(BF_HPA)
                     + log(BF_COMPARTMENTS)

The diagonal is zero by graph convention. A supported graph edge has
P_final > 0.5.

ALPHAFOLD FILTER
The recommended list is conservative. A pair is excluded only when:
- its integrated posterior remains exactly 0.5;
- neither HPA nor COMPARTMENTS records a shared discrete location; and
- at least one of HPA or COMPARTMENTS has profiles for both proteins with
  zero overlap.

Pairs with insufficient coverage or conflicting positive evidence are retained.
The supported-only edge list is also provided as a more aggressive alternative.

RESULT
- Supported colocalization edges: {int(supported.sum()):,}.
- Recommended AlphaFold pairs retained: {int(recommended_retain.sum()):,}.
- Recommended pairs excluded: {int(recommended_exclude.sum()):,}.
- Recommended reduction: {recommended_exclude.mean():.1%}.

FILES
- colocalization_adjacency_matrix.tsv: 848 x 848 protein posterior matrix.
- colocalization_log_bayes_factor_matrix.tsv: additive combined log-BF matrix.
- colocalization_adjacency_matrix_all_nodes.tsv: 868 x 868 audit matrix.
- colocalization_protein_pairs.tsv.gz: all unordered protein pairs and audits.
- supported_colocalization_edges.tsv.gz: graph edges with posterior > 0.5.
- alphafold_candidate_pairs_recommended.tsv.gz: conservative retained jobs.
- alphafold_pairs_excluded_by_localization.tsv.gz: recommended exclusions.
- top_colocalization_edges.tsv: top 2,000 graph edges.
- compartments_*_matrix.tsv: COMPARTMENTS components and coverage.
- colocalization_graph_summary.png: heatmap and pair-reduction summary.
- analysis_summary.json: methods, counts, hashes, validation and caveats.

PROCESSED COMPARTMENTS AUDITS
- {mapping_path}
- {evidence_path}
- {location_path}

IMPORTANT LIMITATION
These factors provide positive or neutral evidence only. A neutral posterior is
not proof that two proteins cannot meet. HPA and COMPARTMENTS are also not
strictly independent because both relate to curated external annotations. Use
the recommended excluded list as a computational screening heuristic, not as a
biological assertion that interaction is impossible.

OFFICIAL SOURCES
{COMPARTMENTS_PAGE_URL}
{COMPARTMENTS_DATA_URL}
{COMPARTMENTS_PAPER_URL}
https://www.proteinatlas.org/humanproteome/subcellular
"""
    readme_path = output_dir / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")

    summary["output_sha256"] = {
        path.name: sha256_file(path)
        for path in (
            mapping_path,
            evidence_path,
            location_path,
            output_dir / "colocalization_adjacency_matrix.tsv",
            output_dir / "colocalization_log_bayes_factor_matrix.tsv",
            all_edges_path,
            supported_path,
            candidate_path,
            excluded_path,
            plot_path,
            readme_path,
        )
    }
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
