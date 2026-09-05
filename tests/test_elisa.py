import pytest
import numpy as np

from lab_tools.analytical.elisa import StandardCurve, predict_unknowns


@pytest.fixture
def standard_data():
    """Create a clean 5-point standard dataset following a typical 4PL shape."""
    conc = np.array([1000.0, 500.0, 100.0, 10.0, 1.0])
    # Generated from 4PL with A=0.1, D=3.0, C=50, B=1.0
    od = 3.0 + (0.1 - 3.0) / (1.0 + (conc / 50.0) ** 1.0)
    return conc, od


@pytest.fixture
def linear_data():
    """Log‑log linear dataset."""
    conc = np.array([100, 10, 1, 0.1])
    m, b = 0.8, 0.5
    od = np.exp(b) * (conc ** m)
    return conc, od


def test_4pl_fit(standard_data):
    conc, od = standard_data
    curve = StandardCurve(conc, od, model="4pl")
    curve.fit()
    assert curve.r_squared > 0.999
    assert "A" in curve.params
    # Check parameters close to true values
    assert abs(curve.params["D"] - 3.0) < 0.1
    assert abs(curve.params["A"] - 0.1) < 0.1


def test_linear_fit(linear_data):
    conc, od = linear_data
    curve = StandardCurve(conc, od, model="linear")
    curve.fit()
    assert curve.r_squared > 0.999
    assert "slope" in curve.params
    assert abs(curve.params["slope"] - 0.8) < 0.01


def test_predict_concentration_4pl(standard_data):
    conc, od = standard_data
    curve = StandardCurve(conc, od, model="4pl")
    curve.fit()
    # Test known concentration
    test_od = 1.5  # roughly corresponds to conc ~ ?
    pred_conc = curve.predict_conc(np.array([test_od]))[0]
    assert np.isfinite(pred_conc)
    assert pred_conc > 0


def test_predict_concentration_linear(linear_data):
    conc, od = linear_data
    curve = StandardCurve(conc, od, model="linear")
    curve.fit()
    test_od = 5.0
    pred = curve.predict_conc(np.array([test_od]))[0]
    expected = (test_od / np.exp(0.5)) ** (1 / 0.8)
    assert np.isclose(pred, expected, rtol=0.01)


def test_outlier_detection(standard_data):
    conc, od = standard_data
    # Introduce an outlier
    od[2] += 2.0
    curve = StandardCurve(conc, od, model="4pl", outlier_threshold=1.2)
    curve.fit()
    assert len(curve.outliers) >= 1


def test_error_handling():
    with pytest.raises(ValueError):
        StandardCurve([1, 2], [1, 2], model="4pl")  # too few points
    with pytest.raises(ValueError):
        StandardCurve([1, 2, 3], [1, 2, 3], model="invalid")
    with pytest.raises(ValueError):
        StandardCurve([1, 2], [1, 2, 3], model="linear")


def test_predict_unknowns(standard_data):
    conc, od = standard_data
    curve = StandardCurve(conc, od, model="4pl")
    curve.fit()
    unknown_ods = np.array([2.0, 0.5])
    sample_ids = np.array(["S1", "S2"])
    result = predict_unknowns(curve, unknown_ods, sample_ids)
    assert len(result) == 2
    assert "Concentration" in result[0]
    assert float(result[0]["Concentration"]) > 0
