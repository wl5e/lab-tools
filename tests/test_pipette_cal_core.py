"""Tests for core calculation functions."""

import pytest
from lab_tools.lab_ops.pipette_cal.core import (
    _water_density,
    get_z_factor,
    mass_to_volume,
    get_tolerances,
    analyze_calibration,
)


def test_water_density_valid():
    # At 20°C, literature value ~0.9982 g/mL
    rho = _water_density(20.0)
    assert abs(rho - 0.99825) < 0.001


def test_water_density_extremes():
    # Within 15-30°C should work
    assert _water_density(15.0) > 0.999
    assert _water_density(30.0) < 1.000


def test_water_density_out_of_range():
    with pytest.raises(ValueError):
        _water_density(14.9)
    with pytest.raises(ValueError):
        _water_density(30.1)


def test_z_factor_typical():
    z = get_z_factor(20.0, 1013.25)
    # Should be ~1.0028
    assert 1.0025 < z < 1.0035


def test_z_factor_pressure_effect():
    z_high = get_z_factor(20.0, 1050.0)
    z_low = get_z_factor(20.0, 980.0)
    # Higher pressure increases air density, which reduces the denominator (rho_water - rho_air)
    # more than it reduces the numerator (1 - rho_air/rho_brass), so Z increases slightly.
    assert z_high > z_low


def test_mass_to_volume():
    v = mass_to_volume(100.0, 20.0)
    assert v == pytest.approx(100.28, abs=0.5)


def test_mass_to_volume_negative_mass():
    with pytest.raises(ValueError):
        mass_to_volume(0, 20.0)
    with pytest.raises(ValueError):
        mass_to_volume(-5, 20.0)


def test_get_tolerances():
    assert get_tolerances(0.5) == (8.0, 8.0)
    assert get_tolerances(1.5) == (3.0, 3.0)
    assert get_tolerances(10.0) == (1.0, 0.8)
    assert get_tolerances(100.0) == (0.6, 0.3)
    assert get_tolerances(500.0) == (0.6, 0.2)
    assert get_tolerances(1000.0) == (0.5, 0.13)


def test_get_tolerances_invalid():
    with pytest.raises(ValueError):
        get_tolerances(0.09)
    with pytest.raises(ValueError):
        get_tolerances(99999)


def test_analyze_calibration_pass():
    # Perfect data should pass tolerances
    raw = [
        {
            "PipetteID": "P20",
            "NominalVolume": 20.0,
            "TargetVolume": 20.0,
            "Weight_mg": 20.0,
            "Temperature_C": 20.0,
            "AirPressure_hPa": 1013.25,
        },
        {
            "PipetteID": "P20",
            "NominalVolume": 20.0,
            "TargetVolume": 20.0,
            "Weight_mg": 20.0,
            "Temperature_C": 20.0,
            "AirPressure_hPa": 1013.25,
        },
        {
            "PipetteID": "P20",
            "NominalVolume": 20.0,
            "TargetVolume": 20.0,
            "Weight_mg": 20.0,
            "Temperature_C": 20.0,
            "AirPressure_hPa": 1013.25,
        },
    ]
    results = analyze_calibration(raw)
    assert len(results) == 1
    assert results[0]["PassFail"] == "PASS"


def test_analyze_calibration_fail_accuracy():
    # 10% error on 100 µL nominal -> should fail
    raw = []
    for _ in range(5):
        raw.append(
            {
                "PipetteID": "P100",
                "NominalVolume": 100.0,
                "TargetVolume": 100.0,
                "Weight_mg": 90.0,  # large error
                "Temperature_C": 20.0,
            }
        )
    results = analyze_calibration(raw)
    assert results[0]["PassFail"] == "FAIL"


def test_analyze_calibration_too_few_measurements():
    raw = [
        {
            "PipetteID": "P100",
            "NominalVolume": 100.0,
            "TargetVolume": 100.0,
            "Weight_mg": 100.0,
            "Temperature_C": 20.0,
        },
        {
            "PipetteID": "P100",
            "NominalVolume": 100.0,
            "TargetVolume": 100.0,
            "Weight_mg": 100.0,
            "Temperature_C": 20.0,
        },
    ]
    results = analyze_calibration(raw)
    assert "FAIL" in results[0]["PassFail"]  # Error message embedded
