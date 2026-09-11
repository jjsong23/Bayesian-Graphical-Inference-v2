#!/usr/bin/env python3
"""Render raw-score distributions and default cutoffs for every evidence stream.

The figures describe the quantities entering the Bayesian scoring kernels, not
the final posterior probabilities.  Eligible nondetections are represented as
zero for streams that support continuous negative evidence.  KinasePredictor
is the one intrinsically record-level source: a display-only maximum raw motif
score is assigned to each eligible kinase--protein pair, while the inference
engine continues to multiply the individual phosphosite prediction factors.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
GUI_DIR = PROJECT_ROOT / "gui"
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

import workflow_engine as workflow  # noqa: E402


DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "results/sensitivity_analysis/evidence_information_2026-08-24"
    / "score_distributions"
)
PRECLOSURE_MATRIX = (
    PROJECT_ROOT
    / "results/edge_characterization/"
    "localization_kinase_predictor_string_hpa_omnipath_stitch/"
    "combined_adjacency_matrix.tsv"
)

BACKGROUND = "#ffffff"
FOREGROUND = "#1f2937"
MUTED = "#6b7280"
GRID = "#d1d5db"
FRAME = "#9ca3af"
BAR = "#2563eb"
CUTOFF = "#dc2626"
CUTOFF_BAND = "#fecaca"
ZERO = "#7c3aed"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    filename = "arialbd.ttf" if bold else "arial.ttf"
    return ImageFont.truetype(str(Path("C:/Windows/Fonts") / filename), size=size)


TITLE_FONT = font(34, bold=True)
SUBTITLE_FONT = font(17)
PANEL_FONT = font(20, bold=True)
AXIS_FONT = font(14)
TICK_FONT = font(13)
NOTE_FONT = font(13)
SMALL_BOLD = font(13, bold=True)


@dataclass
class StreamDistribution:
    stage: str
    stream_id: str
    label: str
    scores: np.ndarray
    cutoffs: np.ndarray
    score_unit: str
    cutoff_label: str
    population_note: str
    inference_note: str


def _finite_nonnegative(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float).ravel()
    return np.maximum(array[np.isfinite(array)], 0.0)


def _compact(value: float) -> str:
    if not math.isfinite(value):
        return "NA"
    absolute = abs(value)
    if absolute >= 10_000 or (0 < absolute < 0.001):
        return f"{value:.2e}"
    if absolute >= 100:
        return f"{value:.0f}"
    if absolute >= 10:
        return f"{value:.1f}"
    return f"{value:.3g}"


def _safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()


def _node_distributions(
    project: Path,
    registry: dict,
) -> list[StreamDistribution]:
    factors = workflow.load_node_factor_catalog(project, registry)
    distributions: list[StreamDistribution] = []
    for definition in registry["node_streams"]:
        values, thresholds = workflow._node_stream_raw_values_and_thresholds(
            project, factors, definition
        )
        eligible, observed = workflow.node_stream_observation_masks(factors, definition)
        if not eligible.any():
            eligible = np.isfinite(values)
        zero_inclusive = np.where(np.isfinite(values), np.maximum(values, 0.0), 0.0)
        scores = zero_inclusive[eligible]
        cutoffs = np.asarray(thresholds, dtype=float)[eligible]
        cutoffs = cutoffs[np.isfinite(cutoffs) & (cutoffs > 0)]
        observed_count = int(np.count_nonzero(observed & eligible))
        handler = definition["normalization"]["handler"]
        score_unit = "raw abundance"
        if "transcript" in definition["id"] or "rna_" in definition["id"]:
            score_unit = "raw RNA abundance"
        if handler in {"kinase_activity", "phosphoprotein_response", "site_level_phosphoprotein_response"}:
            score_unit = "absolute log2 fold change"
        population_note = (
            f"{len(scores):,} eligible candidate nodes; {observed_count:,} observed; "
            f"{int(np.count_nonzero(scores == 0)):,} zero/nondetected"
        )
        distributions.append(
            StreamDistribution(
                stage="node",
                stream_id=definition["id"],
                label=definition["label"],
                scores=_finite_nonnegative(scores),
                cutoffs=_finite_nonnegative(cutoffs),
                score_unit=score_unit,
                cutoff_label="Tq",
                population_note=population_note,
                inference_note=(
                    "Candidate-level inputs; eligible missing values are plotted at zero, "
                    "matching continuous-negative scoring."
                ),
            )
        )
    return distributions


def _upper_eligible(eligible: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    left, right = np.triu_indices_from(eligible, 1)
    keep = eligible[left, right]
    return left[keep], right[keep]


def _matrix_distribution(
    *,
    definition: dict,
    symbols: list[str],
    eligible: np.ndarray,
    raw_matrix: np.ndarray,
    cutoffs: np.ndarray,
    score_unit: str,
    cutoff_label: str,
    note: str,
) -> StreamDistribution:
    left, right = _upper_eligible(eligible)
    scores = _finite_nonnegative(raw_matrix[left, right])
    return StreamDistribution(
        stage="edge",
        stream_id=definition["id"],
        label=definition["label"],
        scores=scores,
        cutoffs=_finite_nonnegative(cutoffs),
        score_unit=score_unit,
        cutoff_label=cutoff_label,
        population_note=(
            f"{len(scores):,} eligible unique pairs; "
            f"{int(np.count_nonzero(scores == 0)):,} zero/no-record"
        ),
        inference_note=note,
    )


def _edge_distributions(
    project: Path,
    registry: dict,
) -> list[StreamDistribution]:
    universe = pd.read_csv(project / workflow.UNIVERSE_RELATIVE, sep="\t", dtype=str).fillna("")
    symbols = universe["symbol"].astype(str).tolist()
    definitions = {item["id"]: item for item in registry["edge_streams"]}
    distributions: list[StreamDistribution] = []

    # mpkCCD localization: dot products, with endpoint-specific Tq values.
    definition = definitions["mpkccd_localization"]
    eligible = workflow.edge_stream_eligibility_matrix(
        project, definition, symbols, graph_metadata=universe
    )
    profiles = pd.read_csv(project / workflow.MPKCCD_PROFILES_RELATIVE, sep="\t").set_index("symbol").reindex(symbols)
    raw = profiles[["1K", "4K", "17K", "200Kp", "200Ks"]].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
    thresholds = pd.to_numeric(profiles["localization_Tq"], errors="coerce").to_numpy(float)
    active_nodes = np.flatnonzero(eligible.any(axis=1))
    distributions.append(
        _matrix_distribution(
            definition=definition,
            symbols=symbols,
            eligible=eligible,
            raw_matrix=raw @ raw.T,
            cutoffs=thresholds[active_nodes],
            score_unit="fraction-profile dot product",
            cutoff_label="endpoint Tq",
            note="The red band summarizes node-specific Tq values used before reciprocal symmetrization.",
        )
    )

    # KinasePredictor: show one display statistic per eligible pair while retaining
    # the per-site threshold population.  This does not alter inference.
    definition = definitions["kinase_predictor"]
    eligible = workflow.edge_stream_eligibility_matrix(
        project, definition, symbols, graph_metadata=universe
    )
    left, right = _upper_eligible(eligible)
    pair_index = {
        tuple(sorted((symbols[i], symbols[j]), key=lambda value: (value.casefold(), value))): position
        for position, (i, j) in enumerate(zip(left, right, strict=True))
    }
    pair_scores = np.zeros(len(left), dtype=float)
    predictions = pd.read_csv(project / workflow.KINASE_PREDICTIONS_RELATIVE, sep="\t")
    used = predictions["used_for_undirected_edge"]
    if used.dtype != bool:
        used = used.astype(str).str.casefold().eq("true")
    predictions = predictions.loc[used].copy()
    prediction_thresholds = pd.to_numeric(predictions["site_Tq_q75"], errors="coerce").to_numpy(float)
    for row in predictions[["undirected_node_a", "undirected_node_b", "raw_score"]].itertuples(index=False):
        key = tuple(sorted((str(row.undirected_node_a), str(row.undirected_node_b)), key=lambda value: (value.casefold(), value)))
        position = pair_index.get(key)
        if position is not None:
            pair_scores[position] = max(pair_scores[position], float(row.raw_score))
    distributions.append(
        StreamDistribution(
            stage="edge",
            stream_id=definition["id"],
            label=definition["label"],
            scores=_finite_nonnegative(pair_scores),
            cutoffs=_finite_nonnegative(prediction_thresholds),
            score_unit="maximum raw motif score per pair (display only)",
            cutoff_label="phosphosite Tq",
            population_note=(
                f"{len(pair_scores):,} eligible unique pairs; "
                f"{int(np.count_nonzero(pair_scores == 0)):,} without a retained prediction"
            ),
            inference_note=(
                "Inference multiplies every retained site-level factor; pairwise maxima are used only "
                "to make this distribution interpretable."
            ),
        )
    )

    # STRING association matrix.
    definition = definitions["string_v12"]
    eligible = workflow.edge_stream_eligibility_matrix(
        project, definition, symbols, graph_metadata=universe
    )
    matrix = workflow._read_aligned_matrix(project / workflow.STRING_SCORE_MATRIX_RELATIVE, symbols)
    distributions.append(
        _matrix_distribution(
            definition=definition,
            symbols=symbols,
            eligible=eligible,
            raw_matrix=matrix,
            cutoffs=np.asarray([workflow.STRING_REFERENCE_SCORE]),
            score_unit="STRING combined score",
            cutoff_label="reference score",
            note="STRING uses score odds relative to 0.041; the line is a reference, not a Gaussian Tq.",
        )
    )

    # Both HPA profile definitions.
    hpa_profiles = pd.read_csv(project / workflow.HPA_PROFILES_RELATIVE, sep="\t").set_index("symbol").reindex(symbols)
    for stream_id, high_confidence in (("hpa_primary", False), ("hpa_high_confidence", True)):
        definition = definitions[stream_id]
        eligible = workflow.edge_stream_eligibility_matrix(
            project, definition, symbols, graph_metadata=universe
        )
        binary, observed, thresholds = workflow._hpa_binary_profiles(hpa_profiles, high_confidence)
        unit = np.zeros_like(binary, dtype=float)
        norms = np.linalg.norm(binary, axis=1)
        valid = observed & (norms > 0)
        unit[valid] = binary[valid] / norms[valid, np.newaxis]
        distributions.append(
            _matrix_distribution(
                definition=definition,
                symbols=symbols,
                eligible=eligible,
                raw_matrix=unit @ unit.T,
                cutoffs=thresholds[valid],
                score_unit="binary-profile cosine similarity",
                cutoff_label="endpoint Tq",
                note="The red band summarizes node-specific Tq values used before reciprocal symmetrization.",
            )
        )

    # OmniPath curation effort.
    definition = definitions["omnipath_core"]
    eligible = workflow.edge_stream_eligibility_matrix(
        project, definition, symbols, graph_metadata=universe
    )
    matrix = workflow._read_aligned_matrix(project / workflow.OMNIPATH_EFFORT_MATRIX_RELATIVE, symbols)
    distributions.append(
        _matrix_distribution(
            definition=definition,
            symbols=symbols,
            eligible=eligible,
            raw_matrix=matrix,
            cutoffs=np.asarray([6.0]),
            score_unit="curation effort",
            cutoff_label="Tq",
            note="All protein pairs are eligible; an absent OmniPath record is plotted as zero.",
        )
    )

    # STITCH messenger--protein confidence scores.
    definition = definitions["stitch_secondary_messenger"]
    eligible = workflow.edge_stream_eligibility_matrix(
        project, definition, symbols, graph_metadata=universe
    )
    matrix = np.zeros((len(symbols), len(symbols)), dtype=float)
    index = {symbol: position for position, symbol in enumerate(symbols)}
    evidence = pd.read_csv(project / workflow.STITCH_UNIVERSE_EDGES_RELATIVE, sep="\t")
    for row in evidence[["node_a", "node_b", "stitch_score"]].itertuples(index=False):
        i, j = index.get(str(row.node_a)), index.get(str(row.node_b))
        if i is not None and j is not None:
            score = float(row.stitch_score)
            matrix[i, j] = max(matrix[i, j], score)
            matrix[j, i] = matrix[i, j]
    distributions.append(
        _matrix_distribution(
            definition=definition,
            symbols=symbols,
            eligible=eligible,
            raw_matrix=matrix,
            cutoffs=np.asarray([workflow.STITCH_REFERENCE_SCORE]),
            score_unit="STITCH combined score",
            cutoff_label="reference score",
            note="STITCH uses score odds relative to its 0.150 reporting floor, not a Gaussian Tq.",
        )
    )

    # The derived closure source has an anchor cutoff instead of an empirical Tq.
    definition = definitions["scaffold_triadic_closure"]
    preclosure = pd.read_csv(PRECLOSURE_MATRIX, sep="\t", index_col=0).reindex(index=symbols, columns=symbols).fillna(0.5).to_numpy(float)
    protein = universe["node_type"].astype(str).str.casefold().eq("protein").to_numpy(bool)
    scaffold = universe["classes"].astype(str).map(
        lambda value: "adaptor_scaffold" in {item.strip() for item in value.split(";") if item.strip()}
    ).to_numpy(bool)
    eligible = (
        np.logical_and.outer(protein, scaffold)
        | np.logical_and.outer(scaffold, protein)
    )
    np.fill_diagonal(eligible, False)
    distributions.append(
        _matrix_distribution(
            definition=definition,
            symbols=symbols,
            eligible=eligible,
            raw_matrix=preclosure,
            cutoffs=np.asarray([0.9]),
            score_unit="pre-closure edge probability",
            cutoff_label="anchor cutoff",
            note="Derived one-pass rule: two protein--scaffold probabilities must each be strictly above 0.90.",
        )
    )
    return distributions


def _transform(values: np.ndarray, cutoffs: np.ndarray) -> tuple[np.ndarray, np.ndarray, str]:
    values = _finite_nonnegative(values)
    cutoffs = _finite_nonnegative(cutoffs)
    positive = values[values > 0]
    median = float(np.median(positive)) if len(positive) else 0.0
    q99 = float(np.quantile(positive, 0.99)) if len(positive) else 0.0
    use_log = bool(len(positive) and (q99 > 100 * max(median, 1e-12) or q99 > 10_000))
    if use_log:
        return np.log10(1.0 + values), np.log10(1.0 + cutoffs), "log10(1 + raw score)"
    return values, cutoffs, "raw score"


def _histogram(values: np.ndarray, bins: int = 48) -> tuple[np.ndarray, np.ndarray]:
    lower = float(np.min(values)) if len(values) else 0.0
    upper = float(np.max(values)) if len(values) else 1.0
    if math.isclose(lower, upper):
        upper = lower + max(abs(lower) * 0.05, 1.0)
    counts, edges = np.histogram(values, bins=bins, range=(lower, upper))
    fractions = counts.astype(float) / max(len(values), 1)
    return fractions, edges


def _render_panel(
    draw: ImageDraw.ImageDraw,
    stream: StreamDistribution,
    box: tuple[int, int, int, int],
) -> None:
    left, top, right, bottom = box
    draw.text((left + 10, top + 6), stream.label, font=PANEL_FONT, fill=FOREGROUND)
    draw.text((left + 10, top + 34), stream.population_note, font=NOTE_FONT, fill=MUTED)
    transformed, transformed_cutoffs, axis_label = _transform(stream.scores, stream.cutoffs)
    fractions, edges = _histogram(transformed)

    plot_left, plot_right = left + 78, right - 24
    plot_top, plot_bottom = top + 72, bottom - 94
    width, height = plot_right - plot_left, plot_bottom - plot_top
    x_min, x_max = float(edges[0]), float(edges[-1])
    y_max = max(float(fractions.max()) * 1.08, 0.01)
    draw.rectangle((plot_left, plot_top, plot_right, plot_bottom), outline=FRAME, width=1)
    for step in range(5):
        fraction = step / 4
        y = plot_bottom - fraction * height
        draw.line((plot_left, y, plot_right, y), fill=GRID, width=1)
        draw.text((plot_left - 7, y), f"{fraction * y_max:.2f}", font=TICK_FONT, fill=FOREGROUND, anchor="rm")

    def x_pixel(value: float) -> float:
        return plot_left + (value - x_min) / (x_max - x_min) * width

    if len(transformed_cutoffs):
        q10, median, q90 = np.quantile(transformed_cutoffs, [0.1, 0.5, 0.9])
        if q90 > q10 + 1e-12:
            band_left = max(plot_left, min(plot_right, x_pixel(float(q10))))
            band_right = max(plot_left, min(plot_right, x_pixel(float(q90))))
            draw.rectangle((band_left, plot_top, band_right, plot_bottom), fill=CUTOFF_BAND)
        line_x = max(plot_left, min(plot_right, x_pixel(float(median))))
        draw.line((line_x, plot_top, line_x, plot_bottom), fill=CUTOFF, width=3)

    for index, fraction in enumerate(fractions):
        x0 = x_pixel(float(edges[index]))
        x1 = x_pixel(float(edges[index + 1]))
        y = plot_bottom - float(fraction) / y_max * height
        draw.rectangle((x0 + 1, y, max(x0 + 2, x1 - 1), plot_bottom), fill=BAR)

    if x_min <= 0 <= x_max and np.count_nonzero(stream.scores == 0):
        zero_x = x_pixel(0.0)
        draw.line((zero_x, plot_top, zero_x, plot_bottom), fill=ZERO, width=2)

    for step in range(5):
        fraction = step / 4
        x = plot_left + fraction * width
        value = x_min + fraction * (x_max - x_min)
        draw.text((x, plot_bottom + 7), _compact(value), font=TICK_FONT, fill=FOREGROUND, anchor="ma")
    draw.text(((plot_left + plot_right) / 2, plot_bottom + 30), axis_label, font=AXIS_FONT, fill=FOREGROUND, anchor="ma")
    draw.text((left + 10, bottom - 52), stream.score_unit, font=SMALL_BOLD, fill=FOREGROUND)
    if len(stream.cutoffs):
        q10, median, q90 = np.quantile(stream.cutoffs, [0.1, 0.5, 0.9])
        cutoff_text = f"{stream.cutoff_label}: median {_compact(float(median))}"
        if q90 > q10 + 1e-12:
            cutoff_text += f" (10-90%: {_compact(float(q10))}-{_compact(float(q90))})"
    else:
        cutoff_text = f"{stream.cutoff_label}: not applicable"
    draw.text((left + 10, bottom - 31), cutoff_text, font=NOTE_FONT, fill=CUTOFF)


def _render_figure(
    streams: list[StreamDistribution],
    *,
    output: Path,
    title: str,
    columns: int,
) -> None:
    panel_width, panel_height = 650, 390
    rows = math.ceil(len(streams) / columns)
    top_margin = 126
    image = Image.new("RGB", (columns * panel_width, top_margin + rows * panel_height), BACKGROUND)
    draw = ImageDraw.Draw(image)
    draw.text((26, 18), title, font=TITLE_FONT, fill=FOREGROUND)
    draw.text(
        (28, 64),
        "Blue bars: raw score inputs. Red: default cutoff/reference (band = 10th-90th percentile). Purple: zero.",
        font=SUBTITLE_FONT,
        fill=MUTED,
    )
    draw.text(
        (28, 89),
        "Eligible nondetections are included at zero where supported; distributions are descriptive and do not replace the inference rules.",
        font=SUBTITLE_FONT,
        fill=MUTED,
    )
    for index, stream in enumerate(streams):
        row, column = divmod(index, columns)
        box = (
            column * panel_width,
            top_margin + row * panel_height,
            (column + 1) * panel_width,
            top_margin + (row + 1) * panel_height,
        )
        _render_panel(draw, stream, box)
    image.save(output, format="PNG", dpi=(300, 300), optimize=True)


def _write_summary(streams: list[StreamDistribution], output: Path) -> None:
    rows: list[dict] = []
    for stream in streams:
        scores = stream.scores
        cutoffs = stream.cutoffs
        rows.append(
            {
                "stage": stream.stage,
                "stream_id": stream.stream_id,
                "stream_label": stream.label,
                "score_count": len(scores),
                "score_zero_count": int(np.count_nonzero(scores == 0)),
                "score_minimum": float(np.min(scores)) if len(scores) else np.nan,
                "score_median": float(np.median(scores)) if len(scores) else np.nan,
                "score_q75": float(np.quantile(scores, 0.75)) if len(scores) else np.nan,
                "score_maximum": float(np.max(scores)) if len(scores) else np.nan,
                "cutoff_label": stream.cutoff_label,
                "cutoff_count": len(cutoffs),
                "cutoff_q10": float(np.quantile(cutoffs, 0.1)) if len(cutoffs) else np.nan,
                "cutoff_median": float(np.median(cutoffs)) if len(cutoffs) else np.nan,
                "cutoff_q90": float(np.quantile(cutoffs, 0.9)) if len(cutoffs) else np.nan,
                "population_note": stream.population_note,
                "inference_note": stream.inference_note,
            }
        )
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    output = args.output_dir.resolve()
    individual = output / "individual"
    individual.mkdir(parents=True, exist_ok=True)
    registry = workflow.load_registry(project)
    node = _node_distributions(project, registry)
    edge = _edge_distributions(project, registry)

    node_output = output / "node_evidence_score_distributions_with_cutoffs.png"
    edge_output = output / "edge_evidence_score_distributions_with_cutoffs.png"
    _render_figure(node, output=node_output, title="Node evidence scores and default cutoffs", columns=3)
    _render_figure(edge, output=edge_output, title="Edge evidence scores and default cutoffs", columns=2)
    print(node_output)
    print(edge_output)

    for stream in node + edge:
        path = individual / f"{stream.stage}_{_safe_filename(stream.stream_id)}_score_distribution.png"
        _render_figure([stream], output=path, title=f"{stream.label}: scores and cutoff", columns=1)
        print(path)
    _write_summary(node + edge, output / "evidence_score_distribution_summary.tsv")
    (output / "README.md").write_text(
        "# Evidence-score distributions\n\n"
        "These PNGs show the raw quantitative inputs and the default cutoff/reference "
        "for all 12 node streams and all 8 edge streams. Eligible nondetections are "
        "included at zero where continuous negative evidence can be applied. A red "
        "band denotes the 10th-90th percentile when the workflow uses item-specific "
        "Tq values; the red line is their median. STRING and STITCH use odds reference "
        "scores rather than Gaussian-kernel Tq values. Scaffold closure uses an anchor "
        "probability cutoff and has no empirical Tq. For KinasePredictor, the plot uses "
        "the maximum raw motif score per eligible pair only as a readable display; the "
        "actual workflow multiplies all retained site-level prediction factors.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
