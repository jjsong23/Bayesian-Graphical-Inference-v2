"""Bayes-factor-like likelihood scores for the signaling-node universe.

The background threshold is calculated from every finite numeric value in the
input Series. Only after that threshold has been calculated are results limited
to gene symbols in the signaling universe.
"""

from __future__ import annotations

from collections.abc import Collection

import numpy as np
import pandas as pd


def complement_minimum_factors(
    values: pd.Series,
    thresholds: float | pd.Series,
    *,
    minimum_factor: float = 0.5,
) -> pd.Series:
    """Score values against scalar or row-specific positive thresholds.

    This is the project's reusable complement-of-minimum scoring kernel::

        max(minimum_factor, 1 - exp(-0.5 * (value / threshold) ** 2))

    A threshold Series is aligned to ``values`` by index, enabling matched
    backgrounds such as site-count-specific phosphoprotein thresholds.
    Missing/non-finite values receive ``minimum_factor``. Thresholds for finite
    values must be positive and finite.
    """
    if not isinstance(values, pd.Series):
        raise TypeError("values must be a pandas Series")
    if values.index.has_duplicates:
        raise ValueError("values index must be unique")
    if not 0.0 <= minimum_factor <= 1.0:
        raise ValueError("minimum_factor must be between 0 and 1")

    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    if isinstance(thresholds, pd.Series):
        if thresholds.index.has_duplicates:
            raise ValueError("thresholds index must be unique")
        threshold_values = pd.to_numeric(thresholds, errors="coerce").reindex(numeric.index).astype(float)
    else:
        threshold_values = pd.Series(float(thresholds), index=numeric.index, dtype=float)

    finite_values = np.isfinite(numeric.to_numpy())
    threshold_array = threshold_values.to_numpy(dtype=float)
    invalid_thresholds = finite_values & (~np.isfinite(threshold_array) | (threshold_array <= 0))
    if invalid_thresholds.any():
        invalid_index = numeric.index[np.flatnonzero(invalid_thresholds)][:5].tolist()
        raise ValueError(f"finite values require positive finite thresholds; invalid at {invalid_index}")

    factors = np.full(len(numeric), minimum_factor, dtype=float)
    if finite_values.any():
        z = numeric.to_numpy(dtype=float)[finite_values] / threshold_array[finite_values]
        raw = 1.0 - np.exp(-0.5 * np.square(z))
        factors[finite_values] = np.maximum(minimum_factor, raw)
    return pd.Series(factors, index=numeric.index, name="bayes_factor")


def continuous_complement_bayes_factors(
    values: pd.Series,
    thresholds: float | pd.Series,
    *,
    neutral_likelihood: float = 0.5,
    minimum_bayes_factor: float = 1e-6,
) -> pd.Series:
    """Return two-sided Bayes evidence from the unfloored complement kernel.

    Unlike :func:`complement_minimum_factors`, this function does not force the
    likelihood to be at least 0.5. It evaluates::

        BF = max(minimum_bayes_factor,
                 (1 - exp(-0.5 * (value / threshold) ** 2))
                 / neutral_likelihood)

    Consequently, weak values have BF below 1 and strong values have BF above
    1. Callers may represent an eligible nondetection as ``value=0``. The small
    positive BF floor replaces the mathematical zero so downstream log-odds
    updates remain finite. Out-of-scope rows should not be passed as zeros;
    they should remain neutral at BF=1 in the calling integration layer.
    """
    if not isinstance(values, pd.Series):
        raise TypeError("values must be a pandas Series")
    if values.index.has_duplicates:
        raise ValueError("values index must be unique")
    if not 0.0 < neutral_likelihood < 1.0:
        raise ValueError("neutral_likelihood must be strictly between 0 and 1")
    if not 0.0 < minimum_bayes_factor <= 1.0:
        raise ValueError("minimum_bayes_factor must be in (0, 1]")

    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    if isinstance(thresholds, pd.Series):
        if thresholds.index.has_duplicates:
            raise ValueError("thresholds index must be unique")
        threshold_values = (
            pd.to_numeric(thresholds, errors="coerce")
            .reindex(numeric.index)
            .astype(float)
        )
    else:
        threshold_values = pd.Series(float(thresholds), index=numeric.index)

    value_array = numeric.to_numpy(float)
    threshold_array = threshold_values.to_numpy(float)
    valid = (
        np.isfinite(value_array)
        & np.isfinite(threshold_array)
        & (threshold_array > 0)
    )
    result = np.full(len(numeric), minimum_bayes_factor, dtype=float)
    if np.any(valid):
        z = np.maximum(value_array[valid], 0.0) / threshold_array[valid]
        likelihood = 1.0 - np.exp(-0.5 * np.square(z))
        result[valid] = np.maximum(
            minimum_bayes_factor,
            likelihood / neutral_likelihood,
        )
    return pd.Series(result, index=numeric.index, name="bayes_factor")


def signaling_bayes_factors(
    values: pd.Series,
    signaling_universe: Collection[str],
    *,
    q: float = 0.75,
    minimum_factor: float = 0.5,
) -> pd.Series:
    """Calculate abundance-based factors for signaling-universe members.

    Parameters
    ----------
    values
        A pandas Series indexed by gene symbol. All finite numeric values are
        used as the empirical background distribution for ``T_q``.
    signaling_universe
        Gene symbols eligible for the returned result.
    q
        Quantile used for ``T_q``, expressed from 0 to 1. The default is 0.75.
    minimum_factor
        Floor applied to the complement-of-minimum-Bayes-factor score.

    Returns
    -------
    pandas.Series
        Factors for input rows whose gene symbols occur in the signaling
        universe. Input order is preserved. Metadata including ``T_q`` and the
        background size is stored in ``result.attrs``.

    Notes
    -----
    For numeric value ``x`` and background threshold ``T_q``, the score is::

        max(minimum_factor, 1 - exp(-0.5 * (x / T_q) ** 2))

    Missing or non-numeric candidate values receive ``minimum_factor`` but do
    not contribute to the background distribution.
    """
    if not isinstance(values, pd.Series):
        raise TypeError("values must be a pandas Series")
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must be between 0 and 1")
    if not 0.0 <= minimum_factor <= 1.0:
        raise ValueError("minimum_factor must be between 0 and 1")
    if values.index.has_duplicates:
        duplicates = values.index[values.index.duplicated()].unique().tolist()
        preview = ", ".join(map(str, duplicates[:5]))
        raise ValueError(f"values index contains duplicate gene symbols: {preview}")

    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    finite_mask = np.isfinite(numeric.to_numpy())
    background = numeric.iloc[np.flatnonzero(finite_mask)]
    if background.empty:
        raise ValueError("values contains no finite numeric background values")

    t_q = float(background.quantile(q))
    if not np.isfinite(t_q) or t_q <= 0:
        raise ValueError(f"T_q must be positive and finite; calculated {t_q!r}")

    universe = {str(symbol).strip() for symbol in signaling_universe if str(symbol).strip()}
    normalized_index = pd.Index(values.index.map(lambda symbol: str(symbol).strip()))
    candidate_mask = normalized_index.isin(universe)
    candidates = numeric.iloc[np.flatnonzero(candidate_mask)].copy()
    candidates.index = normalized_index[candidate_mask]

    result = complement_minimum_factors(
        candidates,
        t_q,
        minimum_factor=minimum_factor,
    )
    result.attrs.update(
        {
            "q": q,
            "T_q": t_q,
            "minimum_factor": minimum_factor,
            "background_size": int(background.size),
            "candidate_size": int(result.size),
            "background_source": "all finite numeric values in the input Series",
        }
    )
    return result


def bayes_update(
    prior: pd.Series,
    factors: pd.Series,
    *,
    missing_factor: float = 0.5,
) -> pd.Series:
    """Update and normalize a probability vector with aligned evidence factors.

    Factor values missing for a node are assigned ``missing_factor``. This
    mirrors the project's convention that an unobserved gene is not evidence
    of absence.
    """
    if not isinstance(prior, pd.Series) or not isinstance(factors, pd.Series):
        raise TypeError("prior and factors must be pandas Series")
    if prior.index.has_duplicates or factors.index.has_duplicates:
        raise ValueError("prior and factors must have unique indices")
    if not 0.0 < missing_factor <= 1.0:
        raise ValueError("missing_factor must be greater than 0 and at most 1")

    numeric_prior = pd.to_numeric(prior, errors="coerce").astype(float)
    if numeric_prior.isna().any() or not np.isfinite(numeric_prior.to_numpy()).all():
        raise ValueError("prior must contain only finite numeric values")
    if (numeric_prior < 0).any() or numeric_prior.sum() <= 0:
        raise ValueError("prior must be nonnegative with a positive sum")
    numeric_prior = numeric_prior / numeric_prior.sum()

    aligned = pd.to_numeric(factors, errors="coerce").reindex(numeric_prior.index)
    aligned = aligned.fillna(missing_factor).astype(float)
    if not np.isfinite(aligned.to_numpy()).all() or (aligned <= 0).any():
        raise ValueError("factors must be positive finite values")

    joint = numeric_prior * aligned
    posterior = joint / joint.sum()
    posterior.name = "posterior_probability"
    posterior.attrs.update(
        {
            "missing_factor": missing_factor,
            "observed_factor_count": int(factors.index.intersection(prior.index).size),
        }
    )
    return posterior


def binary_bayes_factor_update(
    prior: pd.Series,
    factors: pd.Series,
    *,
    missing_factor: float = 1.0,
) -> pd.Series:
    """Update independent binary hypotheses with aligned Bayes factors.

    Each row represents its own ``present`` versus ``not present`` hypothesis.
    Probabilities are therefore *not* normalized across rows::

        posterior_odds = prior_odds * Bayes_factor

    A missing or inapplicable observation is neutral by default (BF = 1).
    This is the appropriate update for the current node-presence model and for
    any other collection of independent Bernoulli hypotheses.
    """
    if not isinstance(prior, pd.Series) or not isinstance(factors, pd.Series):
        raise TypeError("prior and factors must be pandas Series")
    if prior.index.has_duplicates or factors.index.has_duplicates:
        raise ValueError("prior and factors must have unique indices")
    if not np.isfinite(float(missing_factor)) or missing_factor <= 0.0:
        raise ValueError("missing_factor must be positive and finite")

    numeric_prior = pd.to_numeric(prior, errors="coerce").astype(float)
    if numeric_prior.isna().any() or not np.isfinite(numeric_prior.to_numpy()).all():
        raise ValueError("prior must contain only finite numeric values")
    if ((numeric_prior <= 0.0) | (numeric_prior >= 1.0)).any():
        raise ValueError("binary prior probabilities must be strictly between 0 and 1")

    aligned = pd.to_numeric(factors, errors="coerce").reindex(numeric_prior.index)
    aligned = aligned.fillna(float(missing_factor)).astype(float)
    if not np.isfinite(aligned.to_numpy()).all() or (aligned <= 0.0).any():
        raise ValueError("Bayes factors must be positive finite values")

    numerator = numeric_prior * aligned
    denominator = numerator + (1.0 - numeric_prior)
    posterior = numerator / denominator
    posterior.name = "posterior_probability"
    posterior.attrs.update(
        {
            "missing_factor": float(missing_factor),
            "observed_factor_count": int(factors.index.intersection(prior.index).size),
            "hypothesis_model": "independent binary present versus not present",
            "normalized_across_rows": False,
        }
    )
    return posterior


def binary_edge_bayes_update(
    prior: pd.Series,
    link_likelihood: pd.Series,
    *,
    no_link_likelihood: float = 0.5,
) -> pd.Series:
    """Update independent binary edge probabilities with one evidence stream.

    Unlike :func:`bayes_update`, this function does not normalize probabilities
    across candidates.  Each row is a separate Bernoulli edge hypothesis::

        P(edge | evidence) =
            P(edge) P(evidence | edge)
            -----------------------------------------------
            P(edge) P(evidence | edge)
            + (1 - P(edge)) P(evidence | no edge)

    The localization workflow uses the complement-of-minimum factor as
    ``P(evidence | edge)`` and 0.5 as the neutral reference likelihood for
    ``P(evidence | no edge)``.  Consequently, a factor of 0.5 leaves a 0.5
    prior unchanged.
    """
    if not isinstance(prior, pd.Series) or not isinstance(link_likelihood, pd.Series):
        raise TypeError("prior and link_likelihood must be pandas Series")
    if prior.index.has_duplicates or link_likelihood.index.has_duplicates:
        raise ValueError("prior and link_likelihood must have unique indices")
    if not 0.0 < no_link_likelihood <= 1.0:
        raise ValueError("no_link_likelihood must be greater than 0 and at most 1")

    numeric_prior = pd.to_numeric(prior, errors="coerce").astype(float)
    likelihood = pd.to_numeric(link_likelihood, errors="coerce").reindex(numeric_prior.index).astype(float)
    if (
        numeric_prior.isna().any()
        or likelihood.isna().any()
        or not np.isfinite(numeric_prior.to_numpy()).all()
        or not np.isfinite(likelihood.to_numpy()).all()
    ):
        raise ValueError("prior and link_likelihood must contain finite numeric values")
    if ((numeric_prior < 0) | (numeric_prior > 1)).any():
        raise ValueError("prior probabilities must be between 0 and 1")
    if ((likelihood <= 0) | (likelihood > 1)).any():
        raise ValueError("link_likelihood values must be greater than 0 and at most 1")

    numerator = numeric_prior * likelihood
    denominator = numerator + (1.0 - numeric_prior) * no_link_likelihood
    posterior = numerator / denominator
    posterior.name = "posterior_edge_probability"
    posterior.attrs["no_link_likelihood"] = no_link_likelihood
    return posterior
