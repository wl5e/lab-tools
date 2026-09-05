"""Tests for growth model functions."""
import numpy as np
from lab_tools.microbiology.growth_curve.models import logistic, gompertz


def test_logistic_shape():
    t = np.linspace(0, 10, 50)
    y = logistic(t, A=1.0, mu=0.8, lag=2.0)
    assert len(y) == len(t)
    assert np.all(y >= 0)
    assert np.all(y <= 1.0 + 1e-9)
    # logistic should be increasing after lag
    after_lag = t > 2.5
    assert np.all(np.diff(y[after_lag]) >= -1e-9)


def test_gompertz_shape():
    t = np.linspace(0, 10, 50)
    y = gompertz(t, A=1.0, mu=0.8, lag=2.0)
    assert len(y) == len(t)
    assert np.all(y >= 0)
    assert np.all(y <= 1.0 + 1e-9)
    after_lag = t > 2.5
    assert np.all(np.diff(y[after_lag]) >= -1e-9)
