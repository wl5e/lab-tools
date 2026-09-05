"""Core growth curve fitting and parameter extraction."""
import warnings
import numpy as np
from scipy.optimize import curve_fit
from .models import logistic, gompertz

MODEL_MAP = {
    "logistic": logistic,
    "gompertz": gompertz,
}


def _initial_guess(t, y):
    """Generate initial parameter guesses from data."""
    A0 = np.max(y) * 1.05
    # crude derivative using finite differences
    dydt = np.gradient(y, t)
    # find the index of maximum slope
    idx_max = np.argmax(dydt)
    mu0 = max(dydt[idx_max], 1e-6)
    # lag: time where y reaches ~10% of max
    y10 = 0.1 * A0
    try:
        # first point above 10%
        idx10 = np.where(y >= y10)[0][0]
        lag0 = t[idx10] - t[0]
    except IndexError:
        lag0 = t[0]
    lag0 = max(lag0, 0.0)
    return A0, mu0, lag0


def fit_growth_curve(t, y, model="logistic"):
    """
    Fit a growth model to time‑series data.

    Parameters
    ----------
    t : array_like
        Time values.
    y : array_like
        Observed OD values.
    model : str, optional
        Either 'logistic' or 'gompertz' (default 'logistic').

    Returns
    -------
    dict
        Keys: 'A', 'mu', 'lag', 'r_squared', 'model', 'doubling_time'.
        Values are NaN if fitting fails.
    """
    if model not in MODEL_MAP:
        raise ValueError(f"Unknown model '{model}'. Choose from {list(MODEL_MAP.keys())}")

    func = MODEL_MAP[model]
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    # Drop NaN/Inf pairs
    mask = np.isfinite(t) & np.isfinite(y)
    t = t[mask]
    y = y[mask]

    if len(t) < 4:
        raise ValueError(f"Need at least 4 data points, got {len(t)}.")

    if np.max(y) <= 0:
        raise ValueError("All OD values are zero or negative; cannot fit.")

    tmin, tmax = np.min(t), np.max(t)
    p0 = np.array(_initial_guess(t, y))
    bounds = (
        [0.0, 0.0, 0.0],
        [np.inf, np.inf, tmax]
    )

    try:
        popt, pcov = curve_fit(func, t, y, p0=p0, bounds=bounds, maxfev=10000)
    except Exception as exc:
        warnings.warn(f"Fit failed: {exc}")
        return {
            "A": np.nan,
            "mu": np.nan,
            "lag": np.nan,
            "r_squared": np.nan,
            "model": model,
            "doubling_time": np.nan,
        }

    A, mu, lag = popt
    # Compute R²
    y_pred = func(t, *popt)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    if ss_tot == 0:
        r2 = 1.0 if ss_res == 0 else 0.0
    else:
        r2 = 1.0 - ss_res / ss_tot

    doubling_time = np.log(2) / mu if mu > 0 else np.inf

    return {
        "A": A,
        "mu": mu,
        "lag": lag,
        "r_squared": r2,
        "model": model,
        "doubling_time": doubling_time,
    }
