#!/usr/bin/env python3
"""Audit evidence-stream information, agreement, and Tq sensitivity.

The analysis deliberately separates two questions:

1. How much evidence does a stream contain by itself?  This is measured in
   absolute log2 Bayes-factor units (bits of evidence).
2. How much does an enabled stream matter to the integrated posterior?  This
   is measured by removing it from the default model and recording posterior
   shifts and 0.5-decision flips.

Node analyses use the complete node-factor catalog.  Edge analyses use the
fixed, audited 891-node seed universe and its 396,495 unique undirected pairs.
Derived scaffold closure is excluded because it has no independent raw score
or Tq; it is calculated from the already-integrated graph.
"""

from __future__ import annotations

import argparse
import copy
import html
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
GUI_DIR = PROJECT_ROOT / "gui"
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

from workflow_engine import (  # noqa: E402
    FACTOR_EPSILON,
    UNIVERSE_RELATIVE,
    default_configuration,
    edge_stream_eligibility_matrix,
    edge_stream_factor_table,
    load_node_factor_catalog,
    load_registry,
    node_stream_observation_masks,
    node_stream_values,
    stable_expit,
)


TQ_MULTIPLIERS = (0.50, 0.75, 1.00, 1.25, 1.50, 2.00)
CONTINUOUS_FLOOR = 1e-3
DECISION_CUTOFF = 0.5
LOG_TOLERANCE = 1e-12


def bernoulli_kl_bits(posterior: np.ndarray, reference: np.ndarray | float) -> np.ndarray:
    """Return D_KL(Bernoulli(posterior) || Bernoulli(reference)) in bits."""
    epsilon = 1e-15
    p = np.clip(np.asarray(posterior, dtype=float), epsilon, 1.0 - epsilon)
    q = np.clip(np.asarray(reference, dtype=float), epsilon, 1.0 - epsilon)
    return p * np.log2(p / q) + (1.0 - p) * np.log2((1.0 - p) / (1.0 - q))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=PROJECT_ROOT,
        help="Project root containing gui/evidence_registry.json.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: results/sensitivity_analysis/evidence_information_2026-08-24).",
    )
    return parser.parse_args()


def safe_spearman(left: np.ndarray, right: np.ndarray) -> float:
    finite = np.isfinite(left) & np.isfinite(right)
    if finite.sum() < 3:
        return math.nan
    left_f = left[finite]
    right_f = right[finite]
    if np.ptp(left_f) <= LOG_TOLERANCE or np.ptp(right_f) <= LOG_TOLERANCE:
        return math.nan
    left_rank = pd.Series(left_f).rank(method="average").to_numpy(float)
    right_rank = pd.Series(right_f).rank(method="average").to_numpy(float)
    return float(np.corrcoef(left_rank, right_rank)[0, 1])


def evidence_metrics(
    stream_id: str,
    label: str,
    stage: str,
    bayes_factors: np.ndarray,
    eligible: np.ndarray,
    *,
    enabled_by_default: bool,
    weight: float,
    mode: str,
) -> dict[str, Any]:
    if (bayes_factors <= 0).any():
        raise ValueError(f"{stream_id} contains non-positive Bayes factors")
    bits = weight * np.log2(bayes_factors)
    standalone_posterior = stable_expit(weight * np.log(bayes_factors))
    posterior_kl = bernoulli_kl_bits(standalone_posterior, 0.5)
    positive = eligible & (bits > LOG_TOLERANCE)
    negative = eligible & (bits < -LOG_TOLERANCE)
    nonneutral = positive | negative
    eligible_count = int(eligible.sum())
    absolute_bits = np.abs(bits)
    return {
        "stage": stage,
        "stream_id": stream_id,
        "stream_label": label,
        "mode": mode,
        "enabled_by_default": enabled_by_default,
        "weight": weight,
        "hypothesis_count": int(len(bayes_factors)),
        "eligible_count": eligible_count,
        "eligible_fraction": eligible_count / len(bayes_factors),
        "positive_count": int(positive.sum()),
        "negative_count": int(negative.sum()),
        "nonneutral_count": int(nonneutral.sum()),
        "nonneutral_fraction_of_eligible": (
            float(nonneutral.sum() / eligible_count) if eligible_count else math.nan
        ),
        "total_absolute_bits": float(absolute_bits[eligible].sum()),
        "mean_absolute_bits_per_eligible": (
            float(absolute_bits[eligible].mean()) if eligible_count else math.nan
        ),
        "median_absolute_bits_per_eligible": (
            float(np.median(absolute_bits[eligible])) if eligible_count else math.nan
        ),
        "mean_absolute_bits_per_nonneutral": (
            float(absolute_bits[nonneutral].mean()) if nonneutral.any() else math.nan
        ),
        "signed_net_bits": float(bits[eligible].sum()),
        "mean_signed_bits_per_eligible": (
            float(bits[eligible].mean()) if eligible_count else math.nan
        ),
        "maximum_absolute_bits": (
            float(absolute_bits[eligible].max()) if eligible_count else math.nan
        ),
        "total_posterior_kl_bits_from_prior_0_5": float(posterior_kl[eligible].sum()),
        "mean_posterior_kl_bits_per_eligible": (
            float(posterior_kl[eligible].mean()) if eligible_count else math.nan
        ),
        "median_posterior_kl_bits_per_eligible": (
            float(np.median(posterior_kl[eligible])) if eligible_count else math.nan
        ),
        "maximum_posterior_kl_bits": (
            float(posterior_kl[eligible].max()) if eligible_count else math.nan
        ),
    }


def posterior_change_metrics(
    stage: str,
    stream_id: str,
    label: str,
    full_posterior: np.ndarray,
    reduced_posterior: np.ndarray,
) -> dict[str, Any]:
    delta = full_posterior - reduced_posterior
    conditional_kl = bernoulli_kl_bits(full_posterior, reduced_posterior)
    full_positive = full_posterior > DECISION_CUTOFF
    reduced_positive = reduced_posterior > DECISION_CUTOFF
    return {
        "stage": stage,
        "stream_id": stream_id,
        "stream_label": label,
        "hypothesis_count": int(len(full_posterior)),
        "sum_absolute_posterior_shift": float(np.abs(delta).sum()),
        "mean_absolute_posterior_shift": float(np.abs(delta).mean()),
        "median_absolute_posterior_shift": float(np.median(np.abs(delta))),
        "maximum_absolute_posterior_shift": float(np.abs(delta).max()),
        "mean_signed_posterior_shift": float(delta.mean()),
        "decision_flip_count": int((full_positive != reduced_positive).sum()),
        "positive_to_nonpositive_count": int((full_positive & ~reduced_positive).sum()),
        "nonpositive_to_positive_count": int((~full_positive & reduced_positive).sum()),
        "spearman_posterior": safe_spearman(full_posterior, reduced_posterior),
        "sum_conditional_kl_bits": float(conditional_kl.sum()),
        "mean_conditional_kl_bits": float(conditional_kl.mean()),
        "median_conditional_kl_bits": float(np.median(conditional_kl)),
        "maximum_conditional_kl_bits": float(conditional_kl.max()),
    }


def pairwise_agreement(
    stage: str,
    streams: dict[str, dict[str, Any]],
    *,
    mode: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    ids = list(streams)
    for left_index, left_id in enumerate(ids):
        left = streams[left_id]
        left_bits = np.log2(left["bayes_factors"])
        left_eligible = left["eligible"]
        for right_id in ids[left_index + 1 :]:
            right = streams[right_id]
            right_bits = np.log2(right["bayes_factors"])
            right_eligible = right["eligible"]
            common = left_eligible & right_eligible
            left_sign = np.sign(np.where(np.abs(left_bits) > LOG_TOLERANCE, left_bits, 0.0))
            right_sign = np.sign(np.where(np.abs(right_bits) > LOG_TOLERANCE, right_bits, 0.0))
            both_nonneutral = common & (left_sign != 0) & (right_sign != 0)
            concordant = both_nonneutral & (left_sign == right_sign)
            discordant = both_nonneutral & (left_sign == -right_sign)
            union_positive = common & ((left_sign > 0) | (right_sign > 0))
            intersect_positive = common & (left_sign > 0) & (right_sign > 0)
            rows.append(
                {
                    "stage": stage,
                    "mode": mode,
                    "stream_a_id": left_id,
                    "stream_a_label": left["label"],
                    "stream_b_id": right_id,
                    "stream_b_label": right["label"],
                    "common_eligible_count": int(common.sum()),
                    "both_nonneutral_count": int(both_nonneutral.sum()),
                    "concordant_sign_count": int(concordant.sum()),
                    "discordant_sign_count": int(discordant.sum()),
                    "sign_concordance_fraction": (
                        float(concordant.sum() / both_nonneutral.sum())
                        if both_nonneutral.any()
                        else math.nan
                    ),
                    "positive_hit_jaccard": (
                        float(intersect_positive.sum() / union_positive.sum())
                        if union_positive.any()
                        else math.nan
                    ),
                    "spearman_log2_bf_common_eligible": safe_spearman(
                        left_bits[common], right_bits[common]
                    ),
                }
            )
    return pd.DataFrame(rows)


def node_stream_bundle(
    project: Path,
    factors: pd.DataFrame,
    definition: dict[str, Any],
    state: dict[str, Any],
    *,
    continuous: bool,
) -> dict[str, Any]:
    values = node_stream_values(
        project,
        factors,
        definition,
        state,
        continuous_negative=continuous,
        minimum_bayes_factor=CONTINUOUS_FLOOR,
    )
    neutral = float(definition["neutral_value"])
    bayes_factors = values / neutral
    eligible, observed = node_stream_observation_masks(factors, definition)
    if not eligible.any():
        eligible = np.abs(bayes_factors - 1.0) > FACTOR_EPSILON
    return {
        "label": definition["label"],
        "bayes_factors": bayes_factors,
        "eligible": eligible,
        "observed": observed,
    }


def edge_upper_bundle(
    project: Path,
    definition: dict[str, Any],
    state: dict[str, Any],
    symbols: list[str],
    upper_i: np.ndarray,
    upper_j: np.ndarray,
    *,
    continuous: bool,
) -> dict[str, Any]:
    n = len(symbols)
    table = edge_stream_factor_table(
        project,
        definition,
        state,
        symbols,
        include_negative=continuous,
        continuous_negative=continuous,
        minimum_bayes_factor=CONTINUOUS_FLOOR,
    )
    matrix = np.ones((n, n), dtype=float)
    symbol_to_index = {symbol: index for index, symbol in enumerate(symbols)}
    if not table.empty:
        left = table["node_a"].astype(str).map(symbol_to_index)
        right = table["node_b"].astype(str).map(symbol_to_index)
        valid = left.notna() & right.notna()
        li = left.loc[valid].astype(int).to_numpy()
        ri = right.loc[valid].astype(int).to_numpy()
        bf = pd.to_numeric(table.loc[valid, "bayes_factor"], errors="raise").to_numpy(float)
        canonical_i = np.minimum(li, ri)
        canonical_j = np.maximum(li, ri)
        # Factor tables should already contain one row per pair.  Multiplication
        # keeps the helper correct if a future source emits independent rows.
        np.multiply.at(matrix, (canonical_i, canonical_j), bf)
        matrix[canonical_j, canonical_i] = matrix[canonical_i, canonical_j]
    eligibility = edge_stream_eligibility_matrix(project, definition, symbols)
    if continuous:
        missing_eligible = eligibility & (np.abs(matrix - 1.0) <= FACTOR_EPSILON)
        matrix[missing_eligible] = CONTINUOUS_FLOOR
    return {
        "label": definition["label"],
        "bayes_factors": matrix[upper_i, upper_j],
        "eligible": eligibility[upper_i, upper_j],
    }


def integrated_posterior(
    bundles: dict[str, dict[str, Any]],
    active_ids: Iterable[str],
    weights: dict[str, float],
) -> np.ndarray:
    first = next(iter(bundles.values()))
    log_odds = np.zeros(len(first["bayes_factors"]), dtype=float)
    for stream_id in active_ids:
        log_odds += weights[stream_id] * np.log(bundles[stream_id]["bayes_factors"])
    return stable_expit(log_odds)


def conflict_table(
    stage: str,
    identifiers: pd.DataFrame,
    streams: dict[str, dict[str, Any]],
    posterior: np.ndarray,
    *,
    maximum_rows: int | None,
) -> tuple[pd.DataFrame, int]:
    ids = list(streams)
    bits = np.vstack([np.log2(streams[stream_id]["bayes_factors"]) for stream_id in ids])
    eligible = np.vstack([streams[stream_id]["eligible"] for stream_id in ids])
    positive = eligible & (bits > LOG_TOLERANCE)
    negative = eligible & (bits < -LOG_TOLERANCE)
    positive_count = positive.sum(axis=0)
    negative_count = negative.sum(axis=0)
    conflict = (positive_count > 0) & (negative_count > 0)
    indices = np.flatnonzero(conflict)
    total = int(len(indices))
    if not total:
        return pd.DataFrame(), 0
    eligible_bits = np.where(eligible, bits, np.nan)
    maximum = np.nanmax(eligible_bits[:, indices], axis=0)
    minimum = np.nanmin(eligible_bits[:, indices], axis=0)
    rows = identifiers.iloc[indices].reset_index(drop=True).copy()
    rows["stage"] = stage
    rows["positive_stream_count"] = positive_count[indices]
    rows["negative_stream_count"] = negative_count[indices]
    rows["posterior_default_active_continuous"] = posterior[indices]
    rows["maximum_positive_bits"] = maximum
    rows["minimum_negative_bits"] = minimum
    rows["evidence_spread_bits"] = maximum - minimum
    rows["summed_log2_bf_all_streams"] = np.where(eligible[:, indices], bits[:, indices], 0.0).sum(axis=0)
    nonneutral_count = positive_count[indices] + negative_count[indices]
    positive_fraction = positive_count[indices] / nonneutral_count
    negative_fraction = negative_count[indices] / nonneutral_count
    rows["sign_consensus_fraction"] = np.abs(
        positive_count[indices] - negative_count[indices]
    ) / nonneutral_count
    rows["sign_entropy_bits"] = -(
        positive_fraction * np.log2(positive_fraction)
        + negative_fraction * np.log2(negative_fraction)
    )
    rows["positive_streams"] = [
        ";".join(ids[index] for index in np.flatnonzero(positive[:, hypothesis]))
        for hypothesis in indices
    ]
    rows["negative_streams"] = [
        ";".join(ids[index] for index in np.flatnonzero(negative[:, hypothesis]))
        for hypothesis in indices
    ]
    rows = rows.sort_values(
        ["evidence_spread_bits", "positive_stream_count", "negative_stream_count"],
        ascending=[False, False, False],
    )
    if maximum_rows is not None:
        rows = rows.head(maximum_rows)
    return rows, total


def multistream_consensus_table(
    identifiers: pd.DataFrame,
    streams: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    """Summarize higher-order sign agreement across all supplied streams."""
    ids = list(streams)
    bits = np.vstack([np.log2(streams[stream_id]["bayes_factors"]) for stream_id in ids])
    eligible = np.vstack([streams[stream_id]["eligible"] for stream_id in ids])
    positive = eligible & (bits > LOG_TOLERANCE)
    negative = eligible & (bits < -LOG_TOLERANCE)
    positive_count = positive.sum(axis=0)
    negative_count = negative.sum(axis=0)
    eligible_count = eligible.sum(axis=0)
    nonneutral_count = positive_count + negative_count
    consensus = np.full(len(identifiers), np.nan, dtype=float)
    has_signal = nonneutral_count > 0
    consensus[has_signal] = (
        np.abs(positive_count[has_signal] - negative_count[has_signal])
        / nonneutral_count[has_signal]
    )
    entropy = np.full(len(identifiers), np.nan, dtype=float)
    entropy[has_signal] = 0.0
    conflict = (positive_count > 0) & (negative_count > 0)
    positive_fraction = positive_count[conflict] / nonneutral_count[conflict]
    negative_fraction = negative_count[conflict] / nonneutral_count[conflict]
    entropy[conflict] = -(
        positive_fraction * np.log2(positive_fraction)
        + negative_fraction * np.log2(negative_fraction)
    )
    output = identifiers.reset_index(drop=True).copy()
    output["eligible_stream_count"] = eligible_count
    output["positive_stream_count"] = positive_count
    output["negative_stream_count"] = negative_count
    output["nonneutral_stream_count"] = nonneutral_count
    output["conflicted"] = conflict
    output["sign_consensus_fraction"] = consensus
    output["sign_entropy_bits"] = entropy
    output["positive_streams"] = [
        ";".join(ids[index] for index in np.flatnonzero(positive[:, hypothesis]))
        for hypothesis in range(len(identifiers))
    ]
    output["negative_streams"] = [
        ";".join(ids[index] for index in np.flatnonzero(negative[:, hypothesis]))
        for hypothesis in range(len(identifiers))
    ]
    return output


def tq_result_row(
    stage: str,
    stream_id: str,
    label: str,
    control_label: str,
    multiplier: float,
    bayes_factors: np.ndarray,
    eligible: np.ndarray,
    *,
    integrated_active: bool,
    posterior: np.ndarray | None,
    baseline_posterior: np.ndarray | None,
) -> dict[str, Any]:
    bits = np.log2(bayes_factors)
    nonneutral = eligible & (np.abs(bits) > LOG_TOLERANCE)
    row: dict[str, Any] = {
        "stage": stage,
        "stream_id": stream_id,
        "stream_label": label,
        "control_label": control_label,
        "multiplier": multiplier,
        "integrated_active": integrated_active,
        "eligible_count": int(eligible.sum()),
        "standalone_positive_count": int((eligible & (bits > LOG_TOLERANCE)).sum()),
        "standalone_negative_count": int((eligible & (bits < -LOG_TOLERANCE)).sum()),
        "standalone_nonneutral_count": int(nonneutral.sum()),
        "standalone_total_absolute_bits": float(np.abs(bits[eligible]).sum()),
        "standalone_mean_absolute_bits_per_eligible": (
            float(np.abs(bits[eligible]).mean()) if eligible.any() else math.nan
        ),
    }
    if posterior is None or baseline_posterior is None:
        row.update(
            {
                "integrated_posterior_mean": math.nan,
                "integrated_selected_count_gt_0_5": math.nan,
                "mean_absolute_posterior_shift_from_multiplier_1": math.nan,
                "maximum_absolute_posterior_shift_from_multiplier_1": math.nan,
                "decision_flip_count_vs_multiplier_1": math.nan,
                "spearman_vs_multiplier_1": math.nan,
            }
        )
    else:
        delta = posterior - baseline_posterior
        row.update(
            {
                "integrated_posterior_mean": float(posterior.mean()),
                "integrated_selected_count_gt_0_5": int((posterior > DECISION_CUTOFF).sum()),
                "mean_absolute_posterior_shift_from_multiplier_1": float(np.abs(delta).mean()),
                "maximum_absolute_posterior_shift_from_multiplier_1": float(np.abs(delta).max()),
                "decision_flip_count_vs_multiplier_1": int(
                    ((posterior > DECISION_CUTOFF) != (baseline_posterior > DECISION_CUTOFF)).sum()
                ),
                "spearman_vs_multiplier_1": safe_spearman(posterior, baseline_posterior),
            }
        )
    return row


def svg_text(value: object) -> str:
    return html.escape(str(value), quote=True)


def write_information_svg(
    metrics: pd.DataFrame,
    loo: pd.DataFrame,
    stage: str,
    output_path: Path,
) -> None:
    current = metrics.loc[metrics["mode"].eq("current_positive_or_neutral")].copy()
    weight_sorted = current.sort_values("total_absolute_bits", ascending=False)
    kl_sorted = current.sort_values(
        "total_posterior_kl_bits_from_prior_0_5", ascending=False
    )
    loo_sorted = loo.sort_values("sum_conditional_kl_bits", ascending=False)
    width = 2400
    left_x, middle_x, right_x = 20, 820, 1620
    panel_width = 760
    label_width = 315
    row_height = 34
    height = max(len(weight_sorted), len(kl_sorted), len(loo_sorted)) * row_height + 150
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        f'<title>{svg_text(stage.title())} evidence information and leave-one-out influence</title>',
        '<style>text{font-family:Arial,sans-serif;fill:#222;font-size:14px}.title{font-size:20px;font-weight:600}.axis{font-size:13px}.bar-a{fill:#3b82f6}.bar-b{fill:#f59e0b}.grid{stroke:#d1d5db;stroke-width:1}</style>',
        f'<text class="title" x="{left_x}" y="28">{svg_text(stage.title())}: total weight of evidence</text>',
        f'<text class="axis" x="{left_x + label_width}" y="52">Σ |weight × log2 BF| (bits)</text>',
        f'<text class="title" x="{middle_x}" y="28">{svg_text(stage.title())}: realized KL information</text>',
        f'<text class="axis" x="{middle_x + label_width}" y="52">Σ D_KL(posterior || 0.5) (bits)</text>',
        f'<text class="title" x="{right_x}" y="28">Default model: conditional information</text>',
        f'<text class="axis" x="{right_x + label_width}" y="52">Σ D_KL(full || without stream) (bits)</text>',
    ]
    for x, frame, value_column, bar_class in (
        (left_x, weight_sorted, "total_absolute_bits", "bar-a"),
        (middle_x, kl_sorted, "total_posterior_kl_bits_from_prior_0_5", "bar-b"),
        (right_x, loo_sorted, "sum_conditional_kl_bits", "bar-a"),
    ):
        maximum = float(frame[value_column].max()) if len(frame) else 1.0
        maximum = maximum if maximum > 0 else 1.0
        bar_width = panel_width - label_width - 70
        for row_index, row in enumerate(frame.itertuples(index=False)):
            y = 78 + row_index * row_height
            label = getattr(row, "stream_label")
            value = float(getattr(row, value_column))
            scaled = bar_width * value / maximum
            parts.extend(
                [
                    f'<text x="{x}" y="{y + 16}">{svg_text(label)}</text>',
                    f'<rect class="{bar_class}" x="{x + label_width}" y="{y}" width="{scaled:.3f}" height="20"/>',
                    f'<text x="{x + label_width + scaled + 6:.3f}" y="{y + 16}">{value:.4g}</text>',
                ]
            )
    parts.append("</svg>")
    output_path.write_text("\n".join(parts), encoding="utf-8")


def agreement_color(value: float) -> str:
    if not math.isfinite(value):
        return "#e5e7eb"
    # Red -> yellow -> green, kept explicit so the SVG is fully standalone.
    if value <= 0.5:
        fraction = value / 0.5
        red, green, blue = 239, int(68 + (203 - 68) * fraction), int(68 + (82 - 68) * fraction)
    else:
        fraction = (value - 0.5) / 0.5
        red, green, blue = int(239 + (34 - 239) * fraction), int(203 + (197 - 203) * fraction), int(82 + (94 - 82) * fraction)
    return f"#{red:02x}{green:02x}{blue:02x}"


def write_agreement_svg(agreement: pd.DataFrame, labels: dict[str, str], output_path: Path) -> None:
    ids = list(labels)
    matrix = np.full((len(ids), len(ids)), np.nan)
    np.fill_diagonal(matrix, 1.0)
    index = {stream_id: position for position, stream_id in enumerate(ids)}
    for row in agreement.itertuples(index=False):
        left = index[row.stream_a_id]
        right = index[row.stream_b_id]
        matrix[left, right] = row.sign_concordance_fraction
        matrix[right, left] = row.sign_concordance_fraction
    cell = 52
    left = 340
    top = 300
    width = left + cell * len(ids) + 180
    height = top + cell * len(ids) + 90
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        '<title>Sign agreement under continuous negative evidence</title>',
        '<style>text{font-family:Arial,sans-serif;fill:#222;font-size:13px}.title{font-size:20px;font-weight:600}.cell{stroke:#fff;stroke-width:1}</style>',
        '<text class="title" x="20" y="28">Sign agreement under continuous negative evidence</text>',
        '<text x="20" y="52">Fraction concordant among hypotheses where both streams are non-neutral</text>',
    ]
    for index, stream_id in enumerate(ids):
        x = left + index * cell + cell / 2
        y = top - 8
        parts.append(
            f'<text x="{x:.1f}" y="{y}" transform="rotate(-55 {x:.1f} {y})" text-anchor="start">{svg_text(labels[stream_id])}</text>'
        )
        parts.append(
            f'<text x="{left - 8}" y="{top + index * cell + cell * 0.68:.1f}" text-anchor="end">{svg_text(labels[stream_id])}</text>'
        )
    for row in range(len(ids)):
        for column in range(len(ids)):
            value = matrix[row, column]
            x = left + column * cell
            y = top + row * cell
            label = "NA" if not math.isfinite(value) else f"{value:.2f}"
            parts.extend(
                [
                    f'<rect class="cell" x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{agreement_color(value)}"/>',
                    f'<text x="{x + cell / 2}" y="{y + cell * 0.62:.1f}" text-anchor="middle">{label}</text>',
                ]
            )
    parts.append("</svg>")
    output_path.write_text("\n".join(parts), encoding="utf-8")


def write_tq_svg(tq: pd.DataFrame, stage: str, output_path: Path) -> None:
    active = tq.loc[tq["integrated_active"]].copy()
    width, height = 1600, 700
    plot_top, plot_height, plot_width = 90, 430, 560
    left_x, right_x = 110, 900
    colors = ("#2563eb", "#ea580c", "#16a34a", "#9333ea", "#0891b2", "#be123c", "#4f46e5")
    x_min, x_max = min(TQ_MULTIPLIERS), max(TQ_MULTIPLIERS)
    metrics = (
        (left_x, "decision_flip_count_vs_multiplier_1", "0.5-decision flips vs multiplier 1"),
        (right_x, "mean_absolute_posterior_shift_from_multiplier_1", "Mean absolute posterior shift"),
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        f'<title>{svg_text(stage.title())} Tq and reference sensitivity</title>',
        '<style>text{font-family:Arial,sans-serif;fill:#222;font-size:14px}.title{font-size:20px;font-weight:600}.axis{stroke:#374151;stroke-width:1}.grid{stroke:#d1d5db;stroke-width:1}.series{fill:none;stroke-width:2}.point{stroke:#fff;stroke-width:1}</style>',
        f'<text class="title" x="20" y="28">{svg_text(stage.title())} Tq/reference sensitivity</text>',
    ]
    for panel_x, metric, y_label in metrics:
        y_max = float(active[metric].max()) if len(active) else 1.0
        y_max = y_max if y_max > 0 else 1.0
        parts.extend(
            [
                f'<line class="axis" x1="{panel_x}" y1="{plot_top}" x2="{panel_x}" y2="{plot_top + plot_height}"/>',
                f'<line class="axis" x1="{panel_x}" y1="{plot_top + plot_height}" x2="{panel_x + plot_width}" y2="{plot_top + plot_height}"/>',
                f'<text x="{panel_x + plot_width / 2}" y="{plot_top + plot_height + 48}" text-anchor="middle">Tq or reference multiplier</text>',
                f'<text x="{panel_x - 72}" y="{plot_top + plot_height / 2}" transform="rotate(-90 {panel_x - 72} {plot_top + plot_height / 2})" text-anchor="middle">{svg_text(y_label)}</text>',
            ]
        )
        for tick in TQ_MULTIPLIERS:
            x = panel_x + (tick - x_min) / (x_max - x_min) * plot_width
            parts.append(f'<text x="{x:.2f}" y="{plot_top + plot_height + 22}" text-anchor="middle">{tick:g}</text>')
        for tick_index in range(6):
            fraction = tick_index / 5
            y = plot_top + plot_height - fraction * plot_height
            value = y_max * fraction
            parts.extend(
                [
                    f'<line class="grid" x1="{panel_x}" y1="{y:.2f}" x2="{panel_x + plot_width}" y2="{y:.2f}"/>',
                    f'<text x="{panel_x - 8}" y="{y + 5:.2f}" text-anchor="end">{value:.3g}</text>',
                ]
            )
        for series_index, (label, group) in enumerate(active.groupby("stream_label", sort=False)):
            group = group.sort_values("multiplier")
            color = colors[series_index % len(colors)]
            coordinates = []
            for row in group.itertuples(index=False):
                x = panel_x + (float(row.multiplier) - x_min) / (x_max - x_min) * plot_width
                value = float(getattr(row, metric))
                y = plot_top + plot_height - value / y_max * plot_height
                coordinates.append((x, y))
            points = " ".join(f"{x:.2f},{y:.2f}" for x, y in coordinates)
            parts.append(f'<polyline class="series" stroke="{color}" points="{points}"/>')
            for x, y in coordinates:
                parts.append(f'<circle class="point" cx="{x:.2f}" cy="{y:.2f}" r="4" fill="{color}"/>')
    legend_y = 615
    for index, label in enumerate(active["stream_label"].drop_duplicates()):
        column = index % 3
        row = index // 3
        x = 160 + column * 490
        y = legend_y + row * 28
        color = colors[index % len(colors)]
        parts.extend(
            [
                f'<line x1="{x}" y1="{y}" x2="{x + 28}" y2="{y}" stroke="{color}" stroke-width="3"/>',
                f'<text x="{x + 38}" y="{y + 5}">{svg_text(label)}</text>',
            ]
        )
    parts.append("</svg>")
    output_path.write_text("\n".join(parts), encoding="utf-8")


def markdown_table(frame: pd.DataFrame, columns: list[str], rows: int = 10) -> str:
    selected = frame.loc[:, columns].head(rows).copy()
    for column in selected.select_dtypes(include=["float"]).columns:
        selected[column] = selected[column].map(lambda value: "NA" if pd.isna(value) else f"{value:.6g}")
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = ["| " + " | ".join(map(str, row)) + " |" for row in selected.itertuples(index=False, name=None)]
    return "\n".join([header, divider, *body])


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    args = parse_args()
    project = args.project_root.resolve()
    output = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else project / "results/sensitivity_analysis/evidence_information_2026-08-24"
    )
    output.mkdir(parents=True, exist_ok=True)

    registry = load_registry(project)
    config = default_configuration(registry)
    node_definitions = {item["id"]: item for item in registry["node_streams"]}
    edge_definitions = {
        item["id"]: item for item in registry["edge_streams"] if not item.get("derived")
    }
    node_factors = load_node_factor_catalog(project, registry)
    node_ids = node_factors["gene_symbol"].astype(str).rename("gene_symbol").to_frame()
    symbols = (
        pd.read_csv(project / UNIVERSE_RELATIVE, sep="\t", dtype=str)["symbol"]
        .astype(str)
        .tolist()
    )
    upper_i, upper_j = np.triu_indices(len(symbols), 1)
    edge_ids = pd.DataFrame(
        {
            "node_a": np.asarray(symbols, dtype=object)[upper_i],
            "node_b": np.asarray(symbols, dtype=object)[upper_j],
        }
    )
    print(f"Node hypotheses: {len(node_factors):,}; edge hypotheses: {len(edge_ids):,}", flush=True)

    node_default: dict[str, dict[str, Any]] = {}
    node_continuous: dict[str, dict[str, Any]] = {}
    node_metric_rows: list[dict[str, Any]] = []
    for stream_id, definition in node_definitions.items():
        state = config["node_streams"][stream_id]
        default_bundle = node_stream_bundle(project, node_factors, definition, state, continuous=False)
        continuous_bundle = node_stream_bundle(project, node_factors, definition, state, continuous=True)
        node_default[stream_id] = default_bundle
        node_continuous[stream_id] = continuous_bundle
        for mode, bundle in (
            ("current_positive_or_neutral", default_bundle),
            (f"continuous_negative_floor_{CONTINUOUS_FLOOR:g}", continuous_bundle),
        ):
            node_metric_rows.append(
                evidence_metrics(
                    stream_id,
                    definition["label"],
                    "node",
                    bundle["bayes_factors"],
                    bundle["eligible"],
                    enabled_by_default=bool(state["enabled"]),
                    weight=float(state["weight"]),
                    mode=mode,
                )
            )
        print(f"Loaded node stream: {definition['label']}", flush=True)

    edge_default: dict[str, dict[str, Any]] = {}
    edge_continuous: dict[str, dict[str, Any]] = {}
    edge_metric_rows: list[dict[str, Any]] = []
    for stream_id, definition in edge_definitions.items():
        state = config["edge_streams"][stream_id]
        default_bundle = edge_upper_bundle(
            project, definition, state, symbols, upper_i, upper_j, continuous=False
        )
        continuous_bundle = edge_upper_bundle(
            project, definition, state, symbols, upper_i, upper_j, continuous=True
        )
        edge_default[stream_id] = default_bundle
        edge_continuous[stream_id] = continuous_bundle
        for mode, bundle in (
            ("current_positive_or_neutral", default_bundle),
            (f"continuous_negative_floor_{CONTINUOUS_FLOOR:g}", continuous_bundle),
        ):
            edge_metric_rows.append(
                evidence_metrics(
                    stream_id,
                    definition["label"],
                    "edge",
                    bundle["bayes_factors"],
                    bundle["eligible"],
                    enabled_by_default=bool(state["enabled"]),
                    weight=float(state["weight"]),
                    mode=mode,
                )
            )
        print(f"Loaded edge stream: {definition['label']}", flush=True)

    node_metrics = pd.DataFrame(node_metric_rows)
    edge_metrics = pd.DataFrame(edge_metric_rows)

    node_active = [stream_id for stream_id, state in config["node_streams"].items() if state["enabled"]]
    edge_active = [
        stream_id
        for stream_id, state in config["edge_streams"].items()
        if state["enabled"] and stream_id in edge_definitions
    ]
    node_weights = {stream_id: float(config["node_streams"][stream_id]["weight"]) for stream_id in node_definitions}
    edge_weights = {stream_id: float(config["edge_streams"][stream_id]["weight"]) for stream_id in edge_definitions}
    node_full = integrated_posterior(node_default, node_active, node_weights)
    edge_full = integrated_posterior(edge_default, edge_active, edge_weights)
    node_loo = pd.DataFrame(
        [
            posterior_change_metrics(
                "node",
                stream_id,
                node_definitions[stream_id]["label"],
                node_full,
                integrated_posterior(node_default, [item for item in node_active if item != stream_id], node_weights),
            )
            for stream_id in node_active
        ]
    )
    edge_loo = pd.DataFrame(
        [
            posterior_change_metrics(
                "edge",
                stream_id,
                edge_definitions[stream_id]["label"],
                edge_full,
                integrated_posterior(edge_default, [item for item in edge_active if item != stream_id], edge_weights),
            )
            for stream_id in edge_active
        ]
    )

    node_agreement_default = pairwise_agreement("node", node_default, mode="current_positive_or_neutral")
    edge_agreement_default = pairwise_agreement("edge", edge_default, mode="current_positive_or_neutral")
    node_agreement_continuous = pairwise_agreement(
        "node", node_continuous, mode=f"continuous_negative_floor_{CONTINUOUS_FLOOR:g}"
    )
    edge_agreement_continuous = pairwise_agreement(
        "edge", edge_continuous, mode=f"continuous_negative_floor_{CONTINUOUS_FLOOR:g}"
    )

    node_continuous_full = integrated_posterior(node_continuous, node_active, node_weights)
    edge_continuous_full = integrated_posterior(edge_continuous, edge_active, edge_weights)
    node_multistream_consensus = multistream_consensus_table(node_ids, node_continuous)
    edge_multistream_consensus = multistream_consensus_table(edge_ids, edge_continuous)
    node_conflicts, node_conflict_total = conflict_table(
        "node", node_ids, node_continuous, node_continuous_full, maximum_rows=None
    )
    edge_conflicts, edge_conflict_total = conflict_table(
        "edge", edge_ids, edge_continuous, edge_continuous_full, maximum_rows=5000
    )
    node_active_continuous = {
        stream_id: node_continuous[stream_id] for stream_id in node_active
    }
    edge_active_continuous = {
        stream_id: edge_continuous[stream_id] for stream_id in edge_active
    }
    node_active_conflicts, node_active_conflict_total = conflict_table(
        "node",
        node_ids,
        node_active_continuous,
        node_continuous_full,
        maximum_rows=None,
    )
    edge_active_conflicts, edge_active_conflict_total = conflict_table(
        "edge",
        edge_ids,
        edge_active_continuous,
        edge_continuous_full,
        maximum_rows=5000,
    )

    node_tq_rows: list[dict[str, Any]] = []
    for stream_id, definition in node_definitions.items():
        state = config["node_streams"][stream_id]
        for multiplier in TQ_MULTIPLIERS:
            varied_state = copy.deepcopy(state)
            varied_state["tq_multiplier"] = multiplier
            varied = node_stream_bundle(project, node_factors, definition, varied_state, continuous=False)
            integrated = None
            if stream_id in node_active:
                varied_bundles = dict(node_default)
                varied_bundles[stream_id] = varied
                integrated = integrated_posterior(varied_bundles, node_active, node_weights)
            node_tq_rows.append(
                tq_result_row(
                    "node", stream_id, definition["label"], definition["normalization"]["control_label"],
                    multiplier, varied["bayes_factors"], varied["eligible"],
                    integrated_active=stream_id in node_active,
                    posterior=integrated,
                    baseline_posterior=node_full if integrated is not None else None,
                )
            )
        print(f"Swept node Tq: {definition['label']}", flush=True)

    edge_tq_rows: list[dict[str, Any]] = []
    for stream_id, definition in edge_definitions.items():
        state = config["edge_streams"][stream_id]
        for multiplier in TQ_MULTIPLIERS:
            varied_state = copy.deepcopy(state)
            varied_state["tq_multiplier"] = multiplier
            varied = edge_upper_bundle(
                project, definition, varied_state, symbols, upper_i, upper_j, continuous=False
            )
            integrated = None
            if stream_id in edge_active:
                varied_bundles = dict(edge_default)
                varied_bundles[stream_id] = varied
                integrated = integrated_posterior(varied_bundles, edge_active, edge_weights)
            edge_tq_rows.append(
                tq_result_row(
                    "edge", stream_id, definition["label"], definition["normalization"]["control_label"],
                    multiplier, varied["bayes_factors"], varied["eligible"],
                    integrated_active=stream_id in edge_active,
                    posterior=integrated,
                    baseline_posterior=edge_full if integrated is not None else None,
                )
            )
        print(f"Swept edge Tq/reference: {definition['label']}", flush=True)

    node_tq = pd.DataFrame(node_tq_rows)
    edge_tq = pd.DataFrame(edge_tq_rows)

    tables = {
        "node_stream_information.tsv": node_metrics,
        "edge_stream_information_seed891.tsv": edge_metrics,
        "node_leave_one_out.tsv": node_loo,
        "edge_leave_one_out_seed891.tsv": edge_loo,
        "node_pairwise_agreement_current.tsv": node_agreement_default,
        "edge_pairwise_agreement_current_seed891.tsv": edge_agreement_default,
        "node_pairwise_agreement_continuous.tsv": node_agreement_continuous,
        "edge_pairwise_agreement_continuous_seed891.tsv": edge_agreement_continuous,
        "node_conflicts_continuous.tsv": node_conflicts,
        "edge_conflicts_continuous_top5000_seed891.tsv": edge_conflicts,
        "node_conflicts_continuous_default_active.tsv": node_active_conflicts,
        "edge_conflicts_continuous_default_active_top5000_seed891.tsv": edge_active_conflicts,
        "node_tq_sensitivity.tsv": node_tq,
        "edge_tq_sensitivity_seed891.tsv": edge_tq,
    }
    for filename, table in tables.items():
        table.to_csv(output / filename, sep="\t", index=False)
    node_multistream_consensus.to_csv(
        output / "node_multistream_consensus_continuous.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    edge_multistream_consensus.to_csv(
        output / "edge_multistream_consensus_continuous_seed891.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    write_information_svg(node_metrics, node_loo, "node", output / "node_information_contribution.svg")
    write_information_svg(edge_metrics, edge_loo, "edge", output / "edge_information_contribution_seed891.svg")
    write_agreement_svg(
        node_agreement_continuous,
        {stream_id: node_definitions[stream_id]["label"] for stream_id in node_definitions},
        output / "node_agreement_continuous.svg",
    )
    write_agreement_svg(
        edge_agreement_continuous,
        {stream_id: edge_definitions[stream_id]["label"] for stream_id in edge_definitions},
        output / "edge_agreement_continuous_seed891.svg",
    )
    write_tq_svg(node_tq, "node", output / "node_tq_sensitivity.svg")
    write_tq_svg(edge_tq, "edge", output / "edge_tq_sensitivity_seed891.svg")

    node_intrinsic_rank = (
        node_metrics.loc[node_metrics["mode"].eq("current_positive_or_neutral")]
        .sort_values("total_absolute_bits", ascending=False)
        .reset_index(drop=True)
    )
    edge_intrinsic_rank = (
        edge_metrics.loc[edge_metrics["mode"].eq("current_positive_or_neutral")]
        .sort_values("total_absolute_bits", ascending=False)
        .reset_index(drop=True)
    )
    node_loo_rank = node_loo.sort_values("sum_absolute_posterior_shift", ascending=False).reset_index(drop=True)
    edge_loo_rank = edge_loo.sort_values("sum_absolute_posterior_shift", ascending=False).reset_index(drop=True)
    node_tq_max = (
        node_tq.loc[node_tq["integrated_active"]]
        .groupby(["stream_id", "stream_label"], as_index=False)
        .agg(
            maximum_decision_flips=("decision_flip_count_vs_multiplier_1", "max"),
            maximum_mean_absolute_shift=("mean_absolute_posterior_shift_from_multiplier_1", "max"),
        )
        .sort_values("maximum_decision_flips", ascending=False)
    )
    edge_tq_max = (
        edge_tq.loc[edge_tq["integrated_active"]]
        .groupby(["stream_id", "stream_label"], as_index=False)
        .agg(
            maximum_decision_flips=("decision_flip_count_vs_multiplier_1", "max"),
            maximum_mean_absolute_shift=("mean_absolute_posterior_shift_from_multiplier_1", "max"),
        )
        .sort_values("maximum_decision_flips", ascending=False)
    )
    node_agreement_ranked = (
        node_agreement_continuous.loc[
            node_agreement_continuous["both_nonneutral_count"].ge(100)
        ]
        .sort_values("sign_concordance_fraction")
        .reset_index(drop=True)
    )
    edge_agreement_ranked = (
        edge_agreement_continuous.loc[
            edge_agreement_continuous["both_nonneutral_count"].ge(100)
        ]
        .sort_values("sign_concordance_fraction")
        .reset_index(drop=True)
    )
    summary = {
        "analysis_date": "2026-08-24",
        "node_hypotheses": int(len(node_factors)),
        "edge_nodes": int(len(symbols)),
        "unique_undirected_edge_hypotheses": int(len(edge_ids)),
        "continuous_negative_floor_for_conflict_analysis": CONTINUOUS_FLOOR,
        "tq_multipliers": list(TQ_MULTIPLIERS),
        "node_default_active_streams": node_active,
        "edge_default_active_streams": edge_active,
        "node_conflict_count_continuous": node_conflict_total,
        "edge_conflict_count_continuous": edge_conflict_total,
        "node_conflict_count_continuous_default_active": node_active_conflict_total,
        "edge_conflict_count_continuous_default_active": edge_active_conflict_total,
        "node_hypotheses_with_two_or_more_nonneutral_streams": int(
            node_multistream_consensus["nonneutral_stream_count"].ge(2).sum()
        ),
        "edge_hypotheses_with_two_or_more_nonneutral_streams": int(
            edge_multistream_consensus["nonneutral_stream_count"].ge(2).sum()
        ),
        "node_mean_sign_entropy_bits_where_two_or_more_nonneutral": float(
            node_multistream_consensus.loc[
                node_multistream_consensus["nonneutral_stream_count"].ge(2),
                "sign_entropy_bits",
            ].mean()
        ),
        "edge_mean_sign_entropy_bits_where_two_or_more_nonneutral": float(
            edge_multistream_consensus.loc[
                edge_multistream_consensus["nonneutral_stream_count"].ge(2),
                "sign_entropy_bits",
            ].mean()
        ),
        "top_node_intrinsic_total_bits": node_intrinsic_rank.iloc[0]["stream_id"],
        "top_edge_intrinsic_total_bits": edge_intrinsic_rank.iloc[0]["stream_id"],
        "top_node_intrinsic_total_posterior_kl_bits": (
            node_intrinsic_rank.sort_values(
                "total_posterior_kl_bits_from_prior_0_5", ascending=False
            ).iloc[0]["stream_id"]
        ),
        "top_edge_intrinsic_total_posterior_kl_bits": (
            edge_intrinsic_rank.sort_values(
                "total_posterior_kl_bits_from_prior_0_5", ascending=False
            ).iloc[0]["stream_id"]
        ),
        "top_node_leave_one_out": node_loo_rank.iloc[0]["stream_id"],
        "top_edge_leave_one_out": edge_loo_rank.iloc[0]["stream_id"],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = f"""# Evidence-stream information and sensitivity analysis

Date: 2026-08-24

## Scope

- Node stage: {len(node_factors):,} candidate hypotheses and all {len(node_definitions)} quantitative node streams.
- Edge stage: the fixed {len(symbols):,}-node seed universe, giving {len(edge_ids):,} unique undirected pair hypotheses, and all {len(edge_definitions)} primary quantitative edge streams.
- The default integrated model uses {len(node_active)} node streams and {len(edge_active)} edge streams. Optional streams are scored standalone but are not silently added to the default model.
- Scaffold-mediated closure is excluded from information and Tq rankings because it is a derived binary rule with no independent score distribution or Tq.

## Definitions

For hypothesis *i* and stream *s*, the signed weight of evidence is `weight_s * log2(BF_is)` and its direction-free magnitude is `abs(weight_s * log2(BF_is))` bits. A Bayes factor of 2 supplies +1 bit of log-odds evidence; 0.5 supplies -1 bit; 1 is neutral. This quantity is exactly the change in base-2 log posterior odds caused by that stream. It is unbounded and is best understood as evidence strength, not as bounded uncertainty reduction.

The analysis now also reports a Shannon/Kullback-Leibler quantity. Starting from prior 0.5, one stream alone gives `p_s = BF_s^weight / (1 + BF_s^weight)`. Its realized information gain is `D_KL(Bernoulli(p_s) || Bernoulli(0.5)) = 1 - H_binary(p_s)`. This is bounded from 0 bits at p=0.5 to 1 bit as p approaches 0 or 1. It measures how far the posterior moved from the equal-credence prior. It is not expected mutual information because the project does not specify a probability distribution over all possible future datasets.

For leave-one-stream-out analysis, the conditional information contribution is `D_KL(Bernoulli(p_full) || Bernoulli(p_without_s))`. This measures how much the full posterior distribution differs from the posterior obtained without that stream, conditional on all other default streams.

Two rankings answer different questions:

1. **Standalone information** sums absolute evidence bits across eligible hypotheses. This rewards both coverage and strength. Mean bits per eligible hypothesis is also reported so a broad stream is not automatically interpreted as stronger per observation.
2. **Leave-one-stream-out influence** removes one default-enabled stream while holding the prior, other streams, weights, and Tq/reference values fixed. The primary statistic is the summed absolute posterior shift; decision flips across probability 0.5 are also counted.

Agreement is evaluated by the sign of log2(BF). The current positive-or-neutral implementation cannot generate true sign conflicts. Therefore a second diagnostic uses continuous negative evidence with BF floor {CONTINUOUS_FLOOR:g}: eligible low/nondetected measurements can produce BF below 1. A conflict requires at least one positive and one negative stream for the same hypothesis. This scenario is a sensitivity analysis, not a claim that database absence proves biological absence.

Tq sensitivity changes one stream at a time over multipliers {', '.join(map(str, TQ_MULTIPLIERS))}, holding everything else fixed. STRING and STITCH use a reference-odds multiplier rather than a Gaussian Tq, and their rows are labeled `Ref ×`.

## Standalone node information ranking

{markdown_table(node_intrinsic_rank, ['stream_label', 'eligible_count', 'nonneutral_count', 'total_absolute_bits', 'mean_absolute_bits_per_eligible', 'total_posterior_kl_bits_from_prior_0_5', 'mean_posterior_kl_bits_per_eligible'], rows=len(node_intrinsic_rank))}

## Marginal node influence in the default integrated model

{markdown_table(node_loo_rank, ['stream_label', 'sum_absolute_posterior_shift', 'mean_absolute_posterior_shift', 'sum_conditional_kl_bits', 'mean_conditional_kl_bits', 'decision_flip_count', 'spearman_posterior'], rows=len(node_loo_rank))}

## Standalone edge information ranking

{markdown_table(edge_intrinsic_rank, ['stream_label', 'eligible_count', 'nonneutral_count', 'total_absolute_bits', 'mean_absolute_bits_per_eligible', 'total_posterior_kl_bits_from_prior_0_5', 'mean_posterior_kl_bits_per_eligible'], rows=len(edge_intrinsic_rank))}

## Marginal edge influence in the default integrated model

{markdown_table(edge_loo_rank, ['stream_label', 'sum_absolute_posterior_shift', 'mean_absolute_posterior_shift', 'sum_conditional_kl_bits', 'mean_conditional_kl_bits', 'decision_flip_count', 'spearman_posterior'], rows=len(edge_loo_rank))}

## Continuous-negative agreement and conflict

- Nodes with at least one positive and one negative stream: {node_conflict_total:,}.
- Unique undirected edges with at least one positive and one negative stream: {edge_conflict_total:,}.
- Restricting the comparison to the default-enabled model gives {node_active_conflict_total:,} conflicted nodes and {edge_active_conflict_total:,} conflicted unique edges.
- Every conflict is listed for nodes. Edge files are limited to the 5,000 largest evidence spreads; exact totals remain reported above and in `summary.json`.
- Pairwise tables report sign concordance only among hypotheses where both streams are non-neutral, plus positive-hit Jaccard overlap and Spearman correlation of log2 BFs over jointly eligible hypotheses.
- Higher-order agreement is also recorded without enumerating and repeatedly testing every subset. For every node and unique edge, the multistream tables report positive, negative, and non-neutral stream counts; sign consensus `abs(n_positive - n_negative) / n_nonneutral`; and binary sign entropy. Consensus is 1 for unanimity and 0 for an exactly balanced split. Entropy is 0 bits for unanimity and 1 bit for an evenly divided sign vote.

Lowest node-stream sign concordance (restricted to pairs with at least 100 jointly non-neutral nodes):

{markdown_table(node_agreement_ranked, ['stream_a_label', 'stream_b_label', 'both_nonneutral_count', 'sign_concordance_fraction', 'positive_hit_jaccard', 'spearman_log2_bf_common_eligible'], rows=min(5, len(node_agreement_ranked)))}

Highest node-stream sign concordance:

{markdown_table(node_agreement_ranked.sort_values('sign_concordance_fraction', ascending=False), ['stream_a_label', 'stream_b_label', 'both_nonneutral_count', 'sign_concordance_fraction', 'positive_hit_jaccard', 'spearman_log2_bf_common_eligible'], rows=min(5, len(node_agreement_ranked)))}

Lowest edge-stream sign concordance (same minimum-count rule):

{markdown_table(edge_agreement_ranked, ['stream_a_label', 'stream_b_label', 'both_nonneutral_count', 'sign_concordance_fraction', 'positive_hit_jaccard', 'spearman_log2_bf_common_eligible'], rows=min(5, len(edge_agreement_ranked)))}

Highest edge-stream sign concordance:

{markdown_table(edge_agreement_ranked.sort_values('sign_concordance_fraction', ascending=False), ['stream_a_label', 'stream_b_label', 'both_nonneutral_count', 'sign_concordance_fraction', 'positive_hit_jaccard', 'spearman_log2_bf_common_eligible'], rows=min(5, len(edge_agreement_ranked)))}

## Tq/reference sensitivity ranking

Node streams, ranked by the maximum number of probability-0.5 decision flips anywhere in the one-at-a-time sweep:

{markdown_table(node_tq_max, ['stream_label', 'maximum_decision_flips', 'maximum_mean_absolute_shift'], rows=len(node_tq_max))}

Edge streams:

{markdown_table(edge_tq_max, ['stream_label', 'maximum_decision_flips', 'maximum_mean_absolute_shift'], rows=len(edge_tq_max))}

## Interpretation cautions

- Total information is not the same as reliability. It quantifies how much the supplied Bayes factors move odds, not whether the evidence is biologically correct.
- Evidence streams within the same dependence group are correlated. Their standalone values are comparable, but enabling correlated alternatives together would overstate independent evidence.
- Leave-one-out influence is conditional on the current default streams, prior 0.5, weight 1, and multiplier 1. A stream can contain substantial standalone information yet add little marginal posterior movement when another stream already supports the same hypotheses.
- A thresholded decision-flip count depends on the 0.5 cutoff. The continuous posterior-shift statistics are less brittle and should be considered alongside flips.
- Continuous negative evidence treats low or absent observations as evidence against only within each source's explicit eligibility map. Results are sensitive to that modeling assumption and the BF floor.
"""
    (output / "METHODS_AND_RESULTS.md").write_text(report, encoding="utf-8")
    print(f"Wrote analysis to {output}", flush=True)


if __name__ == "__main__":
    main()
