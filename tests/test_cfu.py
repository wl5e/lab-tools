"""Tests for the Dilution Plate CFU Calculator."""

import json
import math
import os
import tempfile

import pytest

from lab_tools.microbiology.cfu import (
    CFUResult,
    DEFAULT_COUNTABLE_RANGE,
    DilutionCFUCalculator,
    DilutionPlate,
    YEAST_MOLD_RANGE,
    _parse_count_value,
    generate_example_csv,
    main,
    parse_csv_input,
    parse_json_input,
)


class TestParseCountValue:
    """Tests for the _parse_count_value helper."""

    def test_numeric(self):
        count, tntc, tftc = _parse_count_value("42")
        assert count == 42
        assert not tntc
        assert not tftc

    def test_tntc_upper(self):
        count, tntc, tftc = _parse_count_value("TNTC")
        assert tntc
        assert not tftc

    def test_tntc_lower(self):
        count, tntc, tftc = _parse_count_value("tntc")
        assert tntc

    def test_tftc(self):
        count, tntc, tftc = _parse_count_value("TFTC")
        assert tftc
        assert not tntc

    def test_zero_is_tftc(self):
        count, tntc, tftc = _parse_count_value("0")
        assert count == 0
        assert tftc

    def test_greater_than_max(self):
        count, tntc, tftc = _parse_count_value(">300")
        assert tntc

    def test_less_than_min(self):
        count, tntc, tftc = _parse_count_value("<30")
        assert tftc

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            _parse_count_value("abc")

    def test_negative_raises(self):
        with pytest.raises(ValueError):
            _parse_count_value("-5")


class TestDilutionCFUCalculator:
    """Tests for the core calculator."""

    def setup_method(self):
        self.calc = DilutionCFUCalculator()

    def test_is_countable(self):
        assert self.calc.is_countable(150)
        assert not self.calc.is_countable(10)
        assert not self.calc.is_countable(350)
        assert self.calc.is_countable(30)
        assert self.calc.is_countable(300)

    def test_basic_calculation_two_dilutions(self):
        plates = [
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=245, replicate_id="A"),
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=231, replicate_id="B"),
            DilutionPlate(dilution_factor=0.001, volume_plated_ml=0.1, colony_count=42, replicate_id="A"),
            DilutionPlate(dilution_factor=0.001, volume_plated_ml=0.1, colony_count=38, replicate_id="B"),
        ]
        result = self.calc.calculate_from_plates(plates)
        assert result.plates_used == 4
        assert result.cfu_per_ml > 0
        assert result.uncertainty_percent > 0
        assert len(result.dilutions_used) == 2

    def test_with_sample_mass(self):
        plates = [
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=150, replicate_id="A"),
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=160, replicate_id="B"),
        ]
        result = self.calc.calculate_from_plates(plates, sample_mass_g=1.0, dilution_blank_ml=9.0)
        assert result.cfu_per_g is not None
        assert result.cfu_per_g > result.cfu_per_ml

    def test_tntc_excluded(self):
        plates = [
            DilutionPlate(dilution_factor=0.1, volume_plated_ml=0.1, colony_count=0, is_tntc=True, replicate_id="A"),
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=100, replicate_id="A"),
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=110, replicate_id="B"),
        ]
        result = self.calc.calculate_from_plates(plates)
        assert result.plates_used == 2
        assert any("TNTC" in n for n in result.notes)

    def test_all_tntc_raises(self):
        plates = [
            DilutionPlate(dilution_factor=0.1, volume_plated_ml=0.1, colony_count=0, is_tntc=True),
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=0, is_tntc=True),
        ]
        with pytest.raises(ValueError, match="No usable plates"):
            self.calc.calculate_from_plates(plates)

    def test_no_plates_raises(self):
        with pytest.raises(ValueError, match="No plates provided"):
            self.calc.calculate_from_plates([])

    def test_fallback_outside_countable_range(self):
        plates = [
            DilutionPlate(dilution_factor=0.0001, volume_plated_ml=0.1, colony_count=8, replicate_id="A"),
            DilutionPlate(dilution_factor=0.0001, volume_plated_ml=0.1, colony_count=12, replicate_id="B"),
        ]
        result = self.calc.calculate_from_plates(plates)
        assert result.plates_used == 2
        assert any("WARNING" in n for n in result.notes)

    def test_single_plate_uncertainty(self):
        plates = [
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=100, replicate_id="A"),
        ]
        result = self.calc.calculate_from_plates(plates)
        assert result.plates_used == 1
        # Poisson-based CV ≈ 1/sqrt(100) = 10%
        assert 9.0 < result.uncertainty_percent < 11.0

    def test_yeast_mold_range(self):
        calc = DilutionCFUCalculator(countable_range=YEAST_MOLD_RANGE)
        assert calc.is_countable(50)
        assert not calc.is_countable(10)

    def test_format_result(self):
        plates = [
            DilutionPlate(dilution_factor=0.01, volume_plated_ml=0.1, colony_count=100, replicate_id="A"),
        ]
        result = self.calc.calculate_from_plates(plates)
        output = self.calc.format_result(result)
        assert "CFU/mL" in output
        assert "Uncertainty" in output


class TestPlatesFromSerialDilution:
    """Tests for the convenience constructor."""

    def test_standard_series(self):
        counts = [-1, -1, 245, 231, 42, 38, 3, 5]
        plates = DilutionCFUCalculator.plates_from_serial_dilution(
            counts, initial_dilution=0.1, num_dilutions=4, replicates=2
        )
        assert len(plates) == 8
        assert plates[0].is_tntc
        assert plates[2].colony_count == 245
        assert plates[6].is_tftc


class TestCSVInput:
    """Tests for CSV parsing."""

    def test_valid_csv(self):
        content = "dilution_factor,volume_plated_ml,colony_count,replicate_id\n" \
                  "0.01,0.1,245,A\n" \
                  "0.01,0.1,TNTC,B\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(content)
            f.flush()
            plates = parse_csv_input(f.name)
        os.unlink(f.name)
        assert len(plates) == 2
        assert plates[0].colony_count == 245
        assert plates[1].is_tntc

    def test_missing_file(self):
        with pytest.raises(FileNotFoundError):
            parse_csv_input("/nonexistent/path.csv")


class TestJSONInput:
    """Tests for JSON parsing."""

    def test_valid_json(self):
        data = [
            {"dilution_factor": 0.01, "volume_plated_ml": 0.1, "colony_count": "245", "replicate_id": "A"},
            {"dilution_factor": 0.001, "volume_plated_ml": 0.1, "colony_count": "TNTC", "replicate_id": "B"},
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            f.flush()
            plates = parse_json_input(f.name)
        os.unlink(f.name)
        assert len(plates) == 2
        assert plates[0].colony_count == 245
        assert plates[1].is_tntc


class TestExampleCSV:
    """Test example CSV generation."""

    def test_generate(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            pass
        try:
            generate_example_csv(f.name)
            with open(f.name, "r") as fh:
                content = fh.read()
            assert "dilution_factor" in content
            assert "TNTC" in content
        finally:
            os.unlink(f.name)


class TestCLI:
    """Integration tests for the command-line interface."""

    def test_example_command(self, capsys):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            tmpname = f.name
        try:
            ret = main(["example", tmpname])
            assert ret == 0
            captured = capsys.readouterr()
            assert "Example CSV" in captured.out
            assert os.path.exists(tmpname)
        finally:
            os.unlink(tmpname)

    def test_quick_command(self, capsys):
        ret = main([
            "quick",
            "--counts", "245,231,42,38",
            "--dilutions", "0.01,0.001",
            "--replicates", "2",
        ])
        assert ret == 0
        captured = capsys.readouterr()
        assert "CFU/mL" in captured.out

    def test_quick_mismatched_counts(self, capsys):
        ret = main([
            "quick",
            "--counts", "245,231",
            "--dilutions", "0.01,0.001",
            "--replicates", "2",
        ])
        assert ret == 1
        captured = capsys.readouterr()
        assert "Expected" in captured.err

    def test_calculate_json_output(self, capsys):
        data = [
            {"dilution_factor": 0.01, "volume_plated_ml": 0.1, "colony_count": "100", "replicate_id": "A"},
            {"dilution_factor": 0.01, "volume_plated_ml": 0.1, "colony_count": "110", "replicate_id": "B"},
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            f.flush()
            tmpname = f.name
        try:
            ret = main(["calculate", tmpname, "--json-output"])
            assert ret == 0
            captured = capsys.readouterr()
            result = json.loads(captured.out)
            assert "cfu_per_ml" in result
            assert result["plates_used"] == 2
        finally:
            os.unlink(tmpname)

    def test_no_command_shows_help(self, capsys):
        ret = main([])
        assert ret == 0
        captured = capsys.readouterr()
        assert "usage" in captured.out.lower() or "calculate" in captured.out

    def test_invalid_file_extension(self, capsys):
        ret = main(["calculate", "data.xlsx"])
        assert ret == 1
        captured = capsys.readouterr()
        assert "extension" in captured.err
