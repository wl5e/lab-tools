"""Tests for the growth curve fitting module."""
import pytest
import numpy as np
from lab_tools.microbiology.growth_curve.analyzer import fit_growth_curve
from lab_tools.microbiology.growth_curve.models import logistic, gompertz


def create_synthetic_data(model_func, params, noise=0.01, num_points=20):
    t = np.linspace(0, 10, num_points)
    true_y = model_func(t, **params)
    rng = np.random.default_rng(42)
    y = true_y + rng.normal(0, noise, size=t.shape)
    y = np.maximum(y, 0.0)  # prevent negative values
    return t, y


@pytest.mark.parametrize("model", ["logistic", "gompertz"])
def test_fit_recovers_parameters(model):
    model_func = {"logistic": logistic, "gompertz": gompertz}[model]
    true_params = {"A": 1.2, "mu": 0.7, "lag": 1.5}
    t, y = create_synthetic_data(model_func, true_params, noise=0.005, num_points=40)
    res = fit_growth_curve(t, y, model=model)
    assert res["model"] == model
    assert not np.isnan(res["A"])
    assert not np.isnan(res["mu"])
    assert not np.isnan(res["lag"])
    assert not np.isnan(res["r_squared"])
    assert not np.isnan(res["doubling_time"])
    # Check that parameters are roughly correct (tolerance liberal due to noise)
    assert np.isclose(res["A"], true_params["A"], rtol=0.2)
    assert np.isclose(res["mu"], true_params["mu"], rtol=0.3)
    assert np.isclose(res["lag"], true_params["lag"], rtol=0.3)
    # R² should be very high
    assert res["r_squared"] > 0.95


def test_fit_fails_on_constant_zero():
    t = np.arange(5, dtype=float)
    y = np.zeros_like(t)
    with pytest.raises(ValueError, match="zero or negative"):
        fit_growth_curve(t, y)


def test_fit_fails_on_few_points():
    t = np.array([0, 1, 2])
    y = np.array([0.1, 0.2, 0.3])
    with pytest.raises(ValueError, match="at least 4"):
        fit_growth_curve(t, y)


def test_fit_handles_nan_values():
    t = np.array([0, 1, 2, 3, 4, 5], dtype=float)
    y = np.array([0.1, np.nan, 0.3, 0.6, 0.8, 1.0])
    res = fit_growth_curve(t, y, model="logistic")
    assert not np.isnan(res["A"])  # should succeed after dropping nan
