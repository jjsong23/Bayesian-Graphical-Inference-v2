#!/usr/bin/env python3
"""Validate ranked signaling paths against a dDAVP phosphoproteomic time course.

The Bayesian path probability remains the primary ranking.  This module appends
temporal-order diagnostics without changing that rank.  Replicate mode uses
per-replicate log2(V/C), intensity-moderated variances, and Monte-Carlo response
times.  Mean mode is retained for comparison with the original contributed
implementation.

The implementation intentionally depends only on NumPy and pandas.  It includes
the small amount of Student-t and rank-correlation math needed by the analysis,
so enabling temporal validation does not introduce an undeclared SciPy runtime
dependency.
"""

from __future__ import annotations

import argparse
import math
import re
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import pandas as pd


DEFAULT_PATH_COL = "path_symbols"
DEFAULT_PATH_DELIM = "->"
DEFAULT_TIMEPOINTS: dict[str, float] = {
    "1 minute": 1.0,
    "2 minute": 2.0,
    "5 minute": 5.0,
    "15 minute": 15.0,
}
DEFAULT_OBSERVED: dict[str, str] = {
    "1 minute": "Observed 1 minute",
    "2 minute": "Observed 2 minute",
    "5 minute": "Observed 5 minute",
    "15 minute": "Observed 15 minute",
}
DEFAULT_GENE_COL = "Gene Symbol"
DEFAULT_SHEET = "All Phosphosites"

RAW_SHEET_MIN: dict[str, float] = {
    "1 min": 1.0,
    "2 min": 2.0,
    "5 min": 5.0,
    "15 min": 15.0,
}
RAW_SHEET_IDX: dict[str, int] = {
    "1 min": 1,
    "2 min": 2,
    "5 min": 3,
    "15 min": 4,
}
RAW_UNIPROT_COL = "Protein UniProt ID"
RAW_POS_COL = "Position in Protein (A)"
RAW_GENE_COL = "Gene Symbol"
P_ADJUST_METHODS = ("within_gene_bonferroni", "none")

ProgressCallback = Callable[[str, float], None]
CancellationCallback = Callable[[], bool]


def _norm(symbol: object) -> str:
    """Normalize a gene symbol for case-insensitive matching."""

    text = str(symbol).strip()
    if not text or text.casefold() in {"nan", "none", "null", "na"}:
        return ""
    return text.upper()


def parse_path(path: object, delimiter: str = DEFAULT_PATH_DELIM) -> list[str]:
    """Split an arrow-delimited path while tolerating surrounding whitespace."""

    pattern = r"\s*" + re.escape(delimiter.strip()) + r"\s*"
    return [item.strip() for item in re.split(pattern, str(path)) if item.strip()]


def _check_cancel(cancel_check: CancellationCallback | None) -> None:
    if cancel_check is not None and cancel_check():
        raise RuntimeError("temporal validation cancelled")


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta using a stable continued fraction."""

    if not a > 0 or not b > 0:
        raise ValueError("beta parameters must be positive")
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0

    def continued_fraction(aa: float, bb: float, xx: float) -> float:
        maximum_iterations = 300
        epsilon = 3e-14
        tiny = np.finfo(float).tiny / epsilon
        qab = aa + bb
        qap = aa + 1.0
        qam = aa - 1.0
        c = 1.0
        d = 1.0 - qab * xx / qap
        if abs(d) < tiny:
            d = tiny
        d = 1.0 / d
        result = d
        for iteration in range(1, maximum_iterations + 1):
            twice = 2 * iteration
            coefficient = (
                iteration * (bb - iteration) * xx
                / ((qam + twice) * (aa + twice))
            )
            d = 1.0 + coefficient * d
            if abs(d) < tiny:
                d = tiny
            c = 1.0 + coefficient / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            result *= d * c

            coefficient = -(
                (aa + iteration)
                * (qab + iteration)
                * xx
                / ((aa + twice) * (qap + twice))
            )
            d = 1.0 + coefficient * d
            if abs(d) < tiny:
                d = tiny
            c = 1.0 + coefficient / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            delta = d * c
            result *= delta
            if abs(delta - 1.0) <= epsilon:
                return result
        raise ArithmeticError("incomplete-beta continued fraction did not converge")

    log_front = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        value = front * continued_fraction(a, b, x) / a
    else:
        value = 1.0 - front * continued_fraction(b, a, 1.0 - x) / b
    return float(min(1.0, max(0.0, value)))


def student_t_two_sided_p(t_statistic: float, degrees_of_freedom: float) -> float:
    """Two-sided Student-t p-value for a nonnegative absolute t statistic."""

    t_value = abs(float(t_statistic))
    degrees = float(degrees_of_freedom)
    if not degrees > 0:
        return float("nan")
    if math.isnan(t_value):
        return float("nan")
    if math.isinf(t_value):
        return 0.0
    x = degrees / (degrees + t_value * t_value)
    return _regularized_incomplete_beta(degrees / 2.0, 0.5, x)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = (start + end - 1) / 2.0
        start = end
    return ranks


def _kendall_tau_b(order_times: Iterable[float]) -> float:
    values = np.asarray(list(order_times), dtype=float)
    if len(values) < 2:
        return float("nan")
    concordant = discordant = tied_time = 0
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            difference = values[right] - values[left]
            if difference > 0:
                concordant += 1
            elif difference < 0:
                discordant += 1
            else:
                tied_time += 1
    untied = concordant + discordant
    denominator = math.sqrt(untied * (untied + tied_time))
    if denominator == 0:
        return float("nan")
    return float((concordant - discordant) / denominator)


def _spearman_rho(order_times: Iterable[float]) -> float:
    values = np.asarray(list(order_times), dtype=float)
    if len(values) < 2 or np.all(values == values[0]):
        return float("nan")
    positions = np.arange(len(values), dtype=float)
    ranks = _average_ranks(values)
    centered_positions = positions - positions.mean()
    centered_ranks = ranks - ranks.mean()
    denominator = math.sqrt(
        float(np.square(centered_positions).sum() * np.square(centered_ranks).sum())
    )
    return (
        float(np.dot(centered_positions, centered_ranks) / denominator)
        if denominator > 0
        else float("nan")
    )


@dataclass
class GeneResponse:
    times: np.ndarray
    values: np.ndarray
    peak: float
    argmax_time: float
    center_of_mass_time: float


def load_gene_responses(
    phospho_path: str | Path,
    *,
    sheet: str = DEFAULT_SHEET,
    gene_col: str = DEFAULT_GENE_COL,
    timepoints: dict[str, float] | None = None,
    observed: dict[str, str] | None = None,
    collapse: str = "best_site",
) -> dict[str, GeneResponse]:
    """Load mean LFCs and collapse each gene to one absolute-LFC trajectory."""

    timepoints = timepoints or DEFAULT_TIMEPOINTS
    observed = observed or DEFAULT_OBSERVED
    frame = pd.read_excel(phospho_path, sheet_name=sheet)
    required = {gene_col, *timepoints, *observed.values()}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(
            "mean temporal workbook is missing columns: " + ", ".join(sorted(missing))
        )
    keys = frame[gene_col].map(_norm)
    minute_columns = list(timepoints)
    minutes = np.asarray([timepoints[column] for column in minute_columns], dtype=float)
    output: dict[str, GeneResponse] = {}
    for gene, subset in frame.loc[keys.ne("")].groupby(keys[keys.ne("")]):
        matrix = np.full((len(subset), len(minute_columns)), np.nan)
        for column_index, column in enumerate(minute_columns):
            values = np.abs(pd.to_numeric(subset[column], errors="coerce").to_numpy())
            is_observed = (
                subset[observed[column]]
                .astype(str)
                .str.strip()
                .str.casefold()
                .eq("yes")
                .to_numpy()
            )
            values[~is_observed] = np.nan
            matrix[:, column_index] = values
        if collapse == "best_site":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                peaks = np.nanmax(matrix, axis=1)
            if np.all(np.isnan(peaks)):
                continue
            envelope = matrix[int(np.nanargmax(peaks)), :]
        elif collapse == "per_timepoint_max":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                envelope = np.nanmax(matrix, axis=0)
        else:
            raise ValueError(f"unknown site-collapse method: {collapse!r}")
        usable = np.isfinite(envelope)
        if not usable.any():
            continue
        times = minutes[usable]
        values = envelope[usable]
        total = float(values.sum())
        output[gene] = GeneResponse(
            times=times,
            values=values,
            peak=float(values.max()),
            argmax_time=float(times[int(np.argmax(values))]),
            center_of_mass_time=(
                float(np.dot(times, values) / total) if total > 0 else float("nan")
            ),
        )
    return output


def _concordance(order_times: list[float]) -> tuple[float, int]:
    if len(order_times) < 2:
        return float("nan"), 0
    concordant = inversions = 0
    for left in range(len(order_times)):
        for right in range(left + 1, len(order_times)):
            if order_times[left] <= order_times[right]:
                concordant += 1
            else:
                inversions += 1
    return concordant / (concordant + inversions), inversions


def score_path_mean(
    symbols: list[str],
    responses: dict[str, GeneResponse],
    minimum_peak: float = 0.0,
) -> dict[str, object]:
    peak_times: list[float] = []
    center_times: list[float] = []
    absent = subthreshold = 0
    trace: list[str] = []
    for symbol in symbols:
        response = responses.get(_norm(symbol))
        if response is None:
            absent += 1
            trace.append(f"{symbol}:absent")
        elif response.peak < minimum_peak:
            subthreshold += 1
            trace.append(
                f"{symbol}:subthreshold({response.center_of_mass_time:.1f},"
                f"{response.peak:.3g})"
            )
        else:
            peak_times.append(response.argmax_time)
            center_times.append(response.center_of_mass_time)
            trace.append(
                f"{symbol}:{response.center_of_mass_time:.1f}({response.peak:.3g})"
            )
    concordance, inversions = _concordance(peak_times)
    return {
        "temporal_n_scored": len(center_times),
        "temporal_n_absent": absent,
        "temporal_n_subthreshold": subthreshold,
        "temporal_peak_time_concordance": concordance,
        "temporal_peak_time_inversions": inversions,
        "temporal_kendall_tau": _kendall_tau_b(center_times),
        "temporal_spearman_rho": _spearman_rho(center_times),
        "temporal_response_time_trace": ";".join(trace),
    }


def rank_paths_mean(
    paths: pd.DataFrame,
    responses: dict[str, GeneResponse],
    *,
    path_col: str = DEFAULT_PATH_COL,
    delimiter: str = DEFAULT_PATH_DELIM,
    minimum_peak: float = 0.0,
) -> pd.DataFrame:
    scores = [
        score_path_mean(parse_path(path, delimiter), responses, minimum_peak)
        for path in paths[path_col]
    ]
    return pd.concat([paths.reset_index(drop=True), pd.DataFrame(scores)], axis=1)


@dataclass
class SiteTrajectory:
    gene: str
    uniprot: str
    site: str
    times: np.ndarray
    log2_ratios: np.ndarray
    intensity: np.ndarray


@dataclass
class ReplicateGeneResponse:
    gene: str
    site: str
    uniprot: str
    passed: bool
    peak_t: float
    peak_p_raw: float
    peak_p: float
    p_adjust_method: str
    gene_test_count: int
    center_of_mass_point: float
    response_time_draws: np.ndarray
    observed_timepoint_count: int


def load_site_trajectories(
    raw_xlsx: str | Path,
    *,
    sheets_min: dict[str, float] | None = None,
    sheets_idx: dict[str, int] | None = None,
    progress: ProgressCallback | None = None,
    cancel_check: CancellationCallback | None = None,
) -> list[SiteTrajectory]:
    """Merge raw per-timepoint sheets into one trajectory per phosphosite."""

    sheets_min = sheets_min or RAW_SHEET_MIN
    sheets_idx = sheets_idx or RAW_SHEET_IDX
    raw_xlsx = Path(raw_xlsx)
    if not raw_xlsx.exists():
        raise FileNotFoundError(f"temporal phosphoproteomic workbook not found: {raw_xlsx}")
    data: dict[tuple[str, str], dict[str, object]] = {}
    for sheet_number, (sheet, minutes) in enumerate(sheets_min.items(), start=1):
        _check_cancel(cancel_check)
        if progress:
            progress(f"Reading temporal phosphoproteomics: {sheet}", sheet_number / len(sheets_min))
        index = sheets_idx[sheet]
        frame = pd.read_excel(raw_xlsx, sheet_name=sheet)
        control_columns = [f"C{index}_Replicate{replicate}_1" for replicate in (1, 2, 3)]
        treated_columns = [f"V{index}_Replicate{replicate}_1" for replicate in (1, 2, 3)]
        required = {
            RAW_UNIPROT_COL,
            RAW_POS_COL,
            RAW_GENE_COL,
            *control_columns,
            *treated_columns,
        }
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(
                f"sheet {sheet!r} is missing columns: " + ", ".join(sorted(missing))
            )
        controls = frame[control_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        treated = frame[treated_columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        observed = (controls > 0) & (treated > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            log2_ratios = np.where(observed, np.log2(treated / controls), np.nan)
        reporter_values = np.concatenate([controls, treated], axis=1)
        reporter_values = np.where(reporter_values > 0, reporter_values, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            intensities = np.nanmedian(reporter_values, axis=1)
        for row_index, row in frame.iterrows():
            gene = _norm(row[RAW_GENE_COL])
            uniprot = str(row[RAW_UNIPROT_COL]).strip()
            site = str(row[RAW_POS_COL]).strip()
            if not gene or not uniprot or uniprot.casefold() == "nan" or not site:
                continue
            key = (uniprot, site)
            record = data.setdefault(
                key,
                {"gene": gene, "uniprot": uniprot, "site": site, "timepoints": {}},
            )
            timepoint_data = record["timepoints"]
            assert isinstance(timepoint_data, dict)
            if minutes in timepoint_data:
                raise ValueError(f"duplicate temporal record for {uniprot} site {site} at {minutes} min")
            timepoint_data[minutes] = (
                log2_ratios[row_index, :].copy(),
                float(intensities[row_index]),
            )

    trajectories: list[SiteTrajectory] = []
    for record in data.values():
        timepoint_data = record["timepoints"]
        assert isinstance(timepoint_data, dict)
        ordered_times = sorted(timepoint_data)
        trajectories.append(
            SiteTrajectory(
                gene=str(record["gene"]),
                uniprot=str(record["uniprot"]),
                site=str(record["site"]),
                times=np.asarray(ordered_times, dtype=float),
                log2_ratios=np.vstack([timepoint_data[time][0] for time in ordered_times]),
                intensity=np.asarray([timepoint_data[time][1] for time in ordered_times], dtype=float),
            )
        )
    if not trajectories:
        raise ValueError("temporal workbook yielded no valid phosphosite trajectories")
    return trajectories


def _timepoint_statistics(values: np.ndarray) -> tuple[float, float, int]:
    observed = values[np.isfinite(values)]
    count = len(observed)
    if count == 0:
        return float("nan"), float("nan"), 0
    if count == 1:
        return float(observed[0]), float("nan"), 1
    return float(observed.mean()), float(observed.var(ddof=1)), count


def fit_intensity_variance_trend(
    log_intensity: np.ndarray,
    raw_variance: np.ndarray,
    *,
    bin_count: int = 25,
    minimum_bin_observations: int = 20,
) -> tuple[Callable[[float], float], tuple[np.ndarray, np.ndarray, float]]:
    """Fit an empirical intensity-dependent prior variance."""

    usable = (
        np.isfinite(log_intensity)
        & np.isfinite(raw_variance)
        & (raw_variance >= 0)
    )
    x_values = log_intensity[usable]
    variances = raw_variance[usable]
    if not len(variances):
        raise ValueError("cannot fit temporal variance trend without replicate variances")
    global_median = max(float(np.median(variances)), np.finfo(float).eps)
    edges = np.unique(np.quantile(x_values, np.linspace(0, 1, bin_count + 1)))
    centers: list[float] = []
    trend_values: list[float] = []
    for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
        if index == len(edges) - 2:
            selected = (x_values >= lower) & (x_values <= upper)
        else:
            selected = (x_values >= lower) & (x_values < upper)
        if int(selected.sum()) >= minimum_bin_observations:
            centers.append(float(np.median(x_values[selected])))
            trend_values.append(
                max(float(np.median(variances[selected])), np.finfo(float).eps)
            )
    center_array = np.asarray(centers, dtype=float)
    trend_array = np.asarray(trend_values, dtype=float)

    def prior_variance(log_value: float) -> float:
        if not np.isfinite(log_value) or not len(center_array):
            return global_median
        return max(
            float(np.interp(log_value, center_array, trend_array)),
            np.finfo(float).eps,
        )

    return prior_variance, (center_array, trend_array, global_median)


def _moderated_variance(
    raw_variance: float,
    observed_count: int,
    prior_variance: float,
    prior_df: float,
) -> tuple[float, float]:
    residual_df = observed_count - 1
    if residual_df <= 0 or not np.isfinite(raw_variance):
        return prior_variance, prior_df
    moderated = (
        prior_df * prior_variance + residual_df * raw_variance
    ) / (prior_df + residual_df)
    return max(float(moderated), np.finfo(float).eps), prior_df + residual_df


def _site_timepoint_frame(
    trajectory: SiteTrajectory,
    prior_variance: Callable[[float], float],
    prior_df: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float] | None:
    times: list[float] = []
    means: list[float] = []
    standard_errors: list[float] = []
    t_statistics: list[float] = []
    p_values: list[float] = []
    for timepoint_index, time in enumerate(trajectory.times):
        mean, raw_variance, count = _timepoint_statistics(
            trajectory.log2_ratios[timepoint_index]
        )
        if count == 0:
            continue
        intensity = trajectory.intensity[timepoint_index]
        log_intensity = (
            math.log10(intensity)
            if np.isfinite(intensity) and intensity > 0
            else float("nan")
        )
        moderated_variance, degrees = _moderated_variance(
            raw_variance,
            count,
            prior_variance(log_intensity),
            prior_df,
        )
        standard_error = math.sqrt(moderated_variance / count)
        t_statistic = abs(mean) / standard_error if standard_error > 0 else float("inf")
        times.append(float(time))
        means.append(mean)
        standard_errors.append(standard_error)
        t_statistics.append(t_statistic)
        p_values.append(student_t_two_sided_p(t_statistic, degrees))
    if not times:
        return None
    return (
        np.asarray(times, dtype=float),
        np.asarray(means, dtype=float),
        np.asarray(standard_errors, dtype=float),
        float(np.nanmax(t_statistics)),
        float(np.nanmin(p_values)),
    )


def build_replicate_responses(
    sites: list[SiteTrajectory],
    *,
    prior_df: float = 4.0,
    alpha: float = 0.05,
    n_boot: int = 2000,
    seed: int = 0,
    p_adjust_method: str = "within_gene_bonferroni",
    progress: ProgressCallback | None = None,
    cancel_check: CancellationCallback | None = None,
) -> tuple[
    dict[str, ReplicateGeneResponse],
    tuple[np.ndarray, np.ndarray, float],
]:
    """Build uncertainty-aware, one-site-per-gene temporal responses.

    ``within_gene_bonferroni`` corrects the representative site's raw peak p
    value for every tested site-timepoint belonging to that gene.  ``none`` is
    retained solely to reproduce the contributed implementation.
    """

    if not prior_df > 0:
        raise ValueError("temporal prior_df must be positive")
    if not 0 < alpha <= 1:
        raise ValueError("temporal alpha must be in (0, 1]")
    if n_boot < 100:
        raise ValueError("temporal Monte-Carlo draws must be at least 100")
    if p_adjust_method not in P_ADJUST_METHODS:
        raise ValueError(
            "unknown temporal p-adjustment method: " + str(p_adjust_method)
        )

    log_intensities: list[float] = []
    raw_variances: list[float] = []
    for trajectory in sites:
        for timepoint_index in range(len(trajectory.times)):
            _, variance, count = _timepoint_statistics(
                trajectory.log2_ratios[timepoint_index]
            )
            intensity = trajectory.intensity[timepoint_index]
            if count >= 2 and np.isfinite(intensity) and intensity > 0:
                log_intensities.append(math.log10(intensity))
                raw_variances.append(variance)
    prior_variance, trend = fit_intensity_variance_trend(
        np.asarray(log_intensities), np.asarray(raw_variances)
    )

    best: dict[str, tuple[float, SiteTrajectory, tuple]] = {}
    gene_test_counts: dict[str, int] = {}
    for site_index, trajectory in enumerate(sites):
        if site_index % 250 == 0:
            _check_cancel(cancel_check)
            if progress:
                progress(
                    "Selecting temporal phosphosites by moderated signal",
                    0.6 * site_index / len(sites),
                )
        frame = _site_timepoint_frame(trajectory, prior_variance, prior_df)
        if frame is None:
            continue
        gene = _norm(trajectory.gene)
        if not gene:
            continue
        gene_test_counts[gene] = gene_test_counts.get(gene, 0) + len(frame[0])
        peak_t = float(frame[3])
        if gene not in best or peak_t > best[gene][0]:
            best[gene] = (peak_t, trajectory, frame)

    generator = np.random.default_rng(seed)
    output: dict[str, ReplicateGeneResponse] = {}
    for gene_index, (gene, (peak_t, trajectory, frame)) in enumerate(best.items()):
        if gene_index % 100 == 0:
            _check_cancel(cancel_check)
            if progress:
                progress(
                    "Propagating temporal response-time uncertainty",
                    0.6 + 0.4 * gene_index / len(best),
                )
        times, means, standard_errors, _, peak_p_raw = frame
        test_count = gene_test_counts[gene]
        peak_p = (
            min(1.0, peak_p_raw * test_count)
            if p_adjust_method == "within_gene_bonferroni"
            else peak_p_raw
        )
        draws = generator.normal(
            loc=means,
            scale=standard_errors,
            size=(n_boot, len(times)),
        )
        absolute_draws = np.abs(draws)
        weight_sums = absolute_draws.sum(axis=1)
        weight_sums = np.where(weight_sums > 0, weight_sums, np.finfo(float).eps)
        response_time_draws = np.dot(absolute_draws, times) / weight_sums
        absolute_means = np.abs(means)
        mean_sum = float(absolute_means.sum())
        center_of_mass = (
            float(np.dot(absolute_means, times) / mean_sum)
            if mean_sum > 0
            else float("nan")
        )
        output[gene] = ReplicateGeneResponse(
            gene=gene,
            site=trajectory.site,
            uniprot=trajectory.uniprot,
            passed=bool(peak_p < alpha),
            peak_t=float(peak_t),
            peak_p_raw=float(peak_p_raw),
            peak_p=float(peak_p),
            p_adjust_method=p_adjust_method,
            gene_test_count=int(test_count),
            center_of_mass_point=center_of_mass,
            response_time_draws=response_time_draws,
            observed_timepoint_count=len(times),
        )
    return output, trend


def _kendall_tau_a(order_times: Iterable[float]) -> float:
    times = np.asarray(list(order_times), dtype=float)
    if len(times) < 2:
        return float("nan")
    pair_count = len(times) * (len(times) - 1) / 2
    sign_sum = 0.0
    for left in range(len(times)):
        for right in range(left + 1, len(times)):
            sign_sum += np.sign(times[right] - times[left])
    return float(sign_sum / pair_count)


def score_path_replicate(
    symbols: list[str],
    responses: dict[str, ReplicateGeneResponse],
) -> dict[str, object]:
    """Score temporal ordering over measured genes that pass the chosen gate."""

    draws: list[np.ndarray] = []
    center_times: list[float] = []
    absent = not_significant = 0
    trace: list[str] = []
    for symbol in symbols:
        response = responses.get(_norm(symbol))
        if response is None:
            absent += 1
            trace.append(f"{symbol}:absent")
        elif not response.passed:
            not_significant += 1
            trace.append(f"{symbol}:not-significant(p={response.peak_p:.2g})")
        else:
            draws.append(response.response_time_draws)
            center_times.append(response.center_of_mass_point)
            trace.append(
                f"{symbol}:{response.center_of_mass_point:.1f}(p={response.peak_p:.2g})"
            )

    scored_count = len(draws)
    result: dict[str, object] = {
        "temporal_n_scored": scored_count,
        "temporal_n_absent": absent,
        "temporal_n_not_significant": not_significant,
        "temporal_soft_precedence": float("nan"),
        "temporal_kendall_tau_point": float("nan"),
        "temporal_kendall_tau_mean": float("nan"),
        "temporal_kendall_tau_low": float("nan"),
        "temporal_kendall_tau_high": float("nan"),
        "temporal_response_time_trace": ";".join(trace),
    }
    if scored_count < 2:
        return result

    response_draws = np.column_stack(draws)
    pair_count = scored_count * (scored_count - 1) / 2
    pair_precedence: list[float] = []
    tau_sign_sum = np.zeros(response_draws.shape[0])
    for left in range(scored_count):
        for right in range(left + 1, scored_count):
            difference = response_draws[:, right] - response_draws[:, left]
            pair_precedence.append(
                float(np.mean(difference > 0) + 0.5 * np.mean(difference == 0))
            )
            tau_sign_sum += np.sign(difference)
    tau_draws = tau_sign_sum / pair_count
    result.update(
        {
            "temporal_soft_precedence": float(np.mean(pair_precedence)),
            "temporal_kendall_tau_point": _kendall_tau_a(center_times),
            "temporal_kendall_tau_mean": float(np.mean(tau_draws)),
            "temporal_kendall_tau_low": float(np.percentile(tau_draws, 2.5)),
            "temporal_kendall_tau_high": float(np.percentile(tau_draws, 97.5)),
        }
    )
    return result


def rank_paths_replicate(
    paths: pd.DataFrame,
    responses: dict[str, ReplicateGeneResponse],
    *,
    path_col: str = DEFAULT_PATH_COL,
    delimiter: str = DEFAULT_PATH_DELIM,
    minimum_scored_nodes: int = 3,
) -> pd.DataFrame:
    """Append replicate-aware scores and a separate temporal evidence rank."""

    scores = [
        score_path_replicate(parse_path(path, delimiter), responses)
        for path in paths[path_col]
    ]
    ranked = pd.concat([paths.reset_index(drop=True), pd.DataFrame(scores)], axis=1)
    informative = (
        ranked["temporal_n_scored"].ge(minimum_scored_nodes)
        & ranked["temporal_kendall_tau_mean"].notna()
    )
    ranked["temporal_order_informative"] = informative
    ranked["temporal_evidence_rank"] = pd.Series(pd.NA, index=ranked.index, dtype="Int64")
    if informative.any():
        order_columns = [
            "temporal_kendall_tau_low",
            "temporal_kendall_tau_mean",
            "temporal_soft_precedence",
        ]
        ascending = [False, False, False]
        if "path_probability_product" in ranked.columns:
            order_columns.append("path_probability_product")
            ascending.append(False)
        ordered_index = ranked.loc[informative].sort_values(
            order_columns,
            ascending=ascending,
            kind="mergesort",
        ).index
        ranked.loc[ordered_index, "temporal_evidence_rank"] = np.arange(
            1, len(ordered_index) + 1
        )
    return ranked


def replicate_response_table(
    responses: dict[str, ReplicateGeneResponse],
) -> pd.DataFrame:
    """Serialize response metadata without the large Monte-Carlo arrays."""

    rows: list[dict[str, object]] = []
    for symbol, response in responses.items():
        record = asdict(response)
        draws = np.asarray(record.pop("response_time_draws"), dtype=float)
        record.update(
            {
                "symbol": symbol,
                "response_time_draw_mean": float(np.mean(draws)),
                "response_time_draw_low": float(np.percentile(draws, 2.5)),
                "response_time_draw_high": float(np.percentile(draws, 97.5)),
            }
        )
        rows.append(record)
    return pd.DataFrame(rows).sort_values("symbol").reset_index(drop=True)


def variance_trend_table(
    trend: tuple[np.ndarray, np.ndarray, float],
) -> pd.DataFrame:
    centers, variances, global_median = trend
    if len(centers):
        return pd.DataFrame(
            {
                "log10_intensity_bin_center": centers,
                "prior_variance": variances,
                "global_median_variance": global_median,
            }
        )
    return pd.DataFrame(
        {
            "log10_intensity_bin_center": [np.nan],
            "prior_variance": [global_median],
            "global_median_variance": [global_median],
        }
    )


def temporal_validation_summary(
    paths: pd.DataFrame,
    responses: dict[str, ReplicateGeneResponse],
    *,
    site_count: int,
    prior_df: float,
    alpha: float,
    n_boot: int,
    seed: int,
    p_adjust_method: str,
    minimum_scored_nodes: int,
) -> dict[str, object]:
    informative = paths["temporal_order_informative"].astype(bool)
    rankable = paths["temporal_kendall_tau_mean"].notna()
    return {
        "enabled": True,
        "mode": "replicate",
        "site_trajectory_count": int(site_count),
        "measured_gene_count": int(len(responses)),
        "genes_passing_gate": int(sum(response.passed for response in responses.values())),
        "prior_df": float(prior_df),
        "alpha": float(alpha),
        "p_adjust_method": p_adjust_method,
        "monte_carlo_draws": int(n_boot),
        "random_seed": int(seed),
        "minimum_scored_nodes_for_temporal_rank": int(minimum_scored_nodes),
        "path_count": int(len(paths)),
        "paths_with_at_least_two_scored_nodes": int(rankable.sum()),
        "temporally_informative_path_count": int(informative.sum()),
        "primary_bayesian_rank_changed": False,
        "interpretation": (
            "Temporal evidence is an annotation and secondary evidence rank. "
            "It does not update edge probabilities or reorder the primary Bayesian rank."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("mean", "replicate"), default="replicate")
    parser.add_argument("--paths", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--path-col", default=DEFAULT_PATH_COL)
    parser.add_argument("--delim", default=DEFAULT_PATH_DELIM)
    parser.add_argument("--phospho", help="mean-LFC workbook for mean mode")
    parser.add_argument("--sheet", default=DEFAULT_SHEET)
    parser.add_argument(
        "--collapse",
        choices=("best_site", "per_timepoint_max"),
        default="best_site",
    )
    parser.add_argument("--min-peak", type=float, default=0.0)
    parser.add_argument("--raw-phospho", help="raw replicate workbook")
    parser.add_argument("--prior-df", type=float, default=4.0)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--p-adjust-method",
        choices=P_ADJUST_METHODS,
        default="within_gene_bonferroni",
    )
    parser.add_argument("--min-scored", type=int, default=3)
    parser.add_argument("--responses-out")
    parser.add_argument("--trend-out")
    arguments = parser.parse_args()

    paths = pd.read_csv(arguments.paths, sep="\t")
    if arguments.mode == "mean":
        if not arguments.phospho:
            parser.error("--phospho is required for mean mode")
        responses = load_gene_responses(
            arguments.phospho,
            sheet=arguments.sheet,
            collapse=arguments.collapse,
        )
        result = rank_paths_mean(
            paths,
            responses,
            path_col=arguments.path_col,
            delimiter=arguments.delim,
            minimum_peak=arguments.min_peak,
        )
    else:
        if not arguments.raw_phospho:
            parser.error("--raw-phospho is required for replicate mode")
        sites = load_site_trajectories(arguments.raw_phospho)
        responses, trend = build_replicate_responses(
            sites,
            prior_df=arguments.prior_df,
            alpha=arguments.alpha,
            n_boot=arguments.n_boot,
            seed=arguments.seed,
            p_adjust_method=arguments.p_adjust_method,
        )
        result = rank_paths_replicate(
            paths,
            responses,
            path_col=arguments.path_col,
            delimiter=arguments.delim,
            minimum_scored_nodes=arguments.min_scored,
        )
        if arguments.responses_out:
            replicate_response_table(responses).to_csv(
                arguments.responses_out, sep="\t", index=False
            )
        if arguments.trend_out:
            variance_trend_table(trend).to_csv(
                arguments.trend_out, sep="\t", index=False
            )
    result.to_csv(arguments.out, sep="\t", index=False)
    print(f"wrote {arguments.out} ({len(result)} paths; mode={arguments.mode})")


if __name__ == "__main__":
    main()
