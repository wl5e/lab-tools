"""Tests for the D-value and Z-value calculator."""

import math

import pytest

from lab_tools.sterility.d_z_f0 import (
    compute_f0,
    compute_d_value,
    compute_z_value,
    linear_regression,
    parse_profile_rows,
    read_csv,
)


def test_linear_regression_perfect():
    x = [0, 5, 10]
    y = [6, 4, 2]  # slope = -0.4
    slope, intercept, r2, se_slope, se_intercept, sigma = linear_regression(x, y)
    assert slope == pytest.approx(-0.4)
    assert intercept == pytest.approx(6.0)
    assert r2 == pytest.approx(1.0)
    assert se_slope == pytest.approx(0.0, abs=1e-10)  # perfect fit
    assert se_intercept == pytest.approx(0.0, abs=1e-10)


def test_linear_regression_variance():
    x = [1, 2, 3, 4]
    y = [2, 4, 5, 7]  # noisy
    slope, intercept, r2, se_slope, se_intercept, sigma = linear_regression(x, y)
    # just check no exception and reasonable values
    assert 1.0 < slope < 2.0
    assert r2 > 0.9


def test_linear_regression_insufficient_points():
    with pytest.raises(ValueError, match="least two"):
        linear_regression([1], [2])
    with pytest.raises(ValueError, match="least two"):
        linear_regression([], [])


def test_linear_regression_no_variance():
    with pytest.raises(ValueError, match="variance"):
        linear_regression([5, 5, 5], [1, 2, 3])


def test_compute_d_value():
    # D=2 min => logN = logN0 - t/D, slope = -1/2
    times = [0, 2, 4, 6]
    logN = [6.0, 5.0, 4.0, 3.0]
    result = compute_d_value(times, logN)
    assert result["D_value_min"] == pytest.approx(2.0, rel=1e-3)
    assert result["r_squared"] == pytest.approx(1.0)
    assert result["se_D"] is not None
    assert result["confidence_interval_95"][0] <= 2.0 <= result["confidence_interval_95"][1]
    assert result["unit"] == "minutes"


def test_compute_d_value_positive_slope():
    times = [0, 1, 2]
    logN = [1, 2, 3]  # growth, not death
    with pytest.raises(ValueError, match="non-negative"):
        compute_d_value(times, logN)


def test_compute_d_value_insufficient():
    with pytest.raises(ValueError):
        compute_d_value([1], [2])
    with pytest.raises(ValueError, match="same length"):
        compute_d_value([1, 2], [3])


def test_compute_d_value_slope_near_zero():
    # almost no change
    times = [0, 10]
    logN = [6, 5.999]
    result = compute_d_value(times, logN)  # slope slightly negative
    # D huge
    assert result["D_value_min"] > 1000


def test_compute_z_value():
    # Known reference: D(121°C)=2 min, Z=10°C => D(131°C)=0.2 min
    temps = [121, 131]
    d_vals = [2.0, 0.2]
    result = compute_z_value(temps, d_vals)
    assert result["Z_value_C"] == pytest.approx(10.0, rel=1e-3)
    assert result["r_squared"] == pytest.approx(1.0)
    assert result["slope"] < 0


def test_compute_z_value_non_negative_slope():
    temps = [100, 110]
    d_vals = [1, 2]  # D increases with temperature, wrong
    with pytest.raises(ValueError, match="non-negative"):
        compute_z_value(temps, d_vals)


def test_compute_z_value_insufficient():
    with pytest.raises(ValueError):
        compute_z_value([121], [2.0])


def test_compute_z_value_length_mismatch():
    with pytest.raises(ValueError, match="same length"):
        compute_z_value([121, 131], [2.0])


def test_compute_f0_reference_hold():
    assert compute_f0([(121.1, 15.0)], z=10.0) == pytest.approx(15.0)


def test_compute_f0_stepped_profile():
    profile = [(121.1, 10.0), (131.1, 10.0)]
    assert compute_f0(profile, z=10.0) == pytest.approx(110.0)


def test_parse_profile_rows_accepts_time_temperature_csv(tmp_path):
    csv_path = tmp_path / "profile.csv"
    csv_path.write_text("time,temperature" + chr(10) + "0,121.1" + chr(10) + "15,121.1" + chr(10))
    rows = read_csv(str(csv_path))
    profile = parse_profile_rows(rows)
    assert profile == [(121.1, 0.0), (121.1, 15.0)]
    assert compute_f0(profile, z=10.0) == pytest.approx(15.0)


def test_parse_profile_rows_rejects_malformed_numeric():
    rows = [{"time": "zero", "temperature": "121.1"}]
    with pytest.raises(ValueError, match="Row 1: time and temperature must be numeric"):
        parse_profile_rows(rows)


def test_parse_profile_rows_rejects_missing_column():
    rows = [{"time": "0"}]
    with pytest.raises(ValueError, match="missing column 'temperature'"):
        parse_profile_rows(rows)


def test_parse_profile_rows_rejects_negative_time():
    rows = [{"time": "-1", "temperature": "121.1"}]
    with pytest.raises(ValueError, match="time must be non-negative"):
        parse_profile_rows(rows)


def test_f0_cli_prints_f0(tmp_path, capsys):
    csv_path = tmp_path / "profile.csv"
    csv_path.write_text("time,temperature\n0,121.1\n15,121.1\n", encoding="utf-8")
    from lab_tools.cli import main as cli_main

    rc = cli_main(["d-z-f0", "f0", "--input", str(csv_path), "--z", "10"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "F0: 15.000 minutes" in out
