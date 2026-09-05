"""Tests for CSV input module."""

import pytest
import tempfile
from pathlib import Path
from lab_tools.lab_ops.pipette_cal.io import read_calibration_data


def _write_csv(content: str) -> str:
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="", encoding="utf-8-sig")
    tmp.write(content)
    tmp.close()
    return tmp.name


def test_read_valid_file():
    csv = _write_csv(
        "PipetteID,NominalVolume,TargetVolume,Weight_mg,Temperature_C,AirPressure_hPa\n"
        "P10,10,5,5.0,21.0,1015\n"
    )
    data = read_calibration_data(csv)
    Path(csv).unlink()
    assert len(data) == 1
    row = data[0]
    assert row["PipetteID"] == "P10"
    assert row["NominalVolume"] == 10.0
    assert row["TargetVolume"] == 5.0
    assert row["Weight_mg"] == 5.0
    assert row["Temperature_C"] == 21.0
    assert row["AirPressure_hPa"] == 1015.0


def test_missing_optional_pressure():
    csv = _write_csv(
        "PipetteID,NominalVolume,TargetVolume,Weight_mg,Temperature_C\n"
        "P10,10,5,5.0,21.0\n"
    )
    data = read_calibration_data(csv)
    Path(csv).unlink()
    assert data[0]["AirPressure_hPa"] == 1013.25


def test_missing_required_column():
    csv = _write_csv(
        "PipetteID,NominalVolume,Weight_mg,Temperature_C\n"
        "P10,10,5,21.0\n"
    )
    with pytest.raises(KeyError):
        read_calibration_data(csv)
    Path(csv).unlink()


def test_invalid_volume():
    csv = _write_csv(
        "PipetteID,NominalVolume,TargetVolume,Weight_mg,Temperature_C\n"
        "P10,10,15,15,21.0\n"
    )
    with pytest.raises(ValueError, match="Invalid nominal"):
        read_calibration_data(csv)
    Path(csv).unlink()


def test_non_numeric_weight():
    csv = _write_csv(
        "PipetteID,NominalVolume,TargetVolume,Weight_mg,Temperature_C\n"
        "P10,10,5,abc,21.0\n"
    )
    with pytest.raises(ValueError):
        read_calibration_data(csv)
    Path(csv).unlink()


def test_file_not_found():
    with pytest.raises(FileNotFoundError):
        read_calibration_data("nonexistent_xyz.csv")
