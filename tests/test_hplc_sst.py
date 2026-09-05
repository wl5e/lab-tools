import io
import math
import statistics
import tempfile
from pathlib import Path

import pytest

from lab_tools.analytical.hplc_sst import (
    calc_plates_half_height,
    calc_plates_base_width,
    calc_tailing_from_asymmetry,
    calc_resolution_base,
    calc_resolution_half,
    process_data,
    read_csv,
    validate_columns,
    to_float,
)


# ---------- Unit tests for calculation functions ----------
def test_plates_half_height():
    assert math.isclose(calc_plates_half_height(2.5, 0.12), 2404.5139, rel_tol=1e-4)
    with pytest.raises(ZeroDivisionError):
        calc_plates_half_height(2.0, 0.0)


def test_plates_base_width():
    assert math.isclose(calc_plates_base_width(2.5, 0.2), 2500.0)
    with pytest.raises(ZeroDivisionError):
        calc_plates_base_width(2.0, 0.0)


def test_tailing_from_asymmetry():
    assert math.isclose(calc_tailing_from_asymmetry(1.1), 1.05)
    assert math.isclose(calc_tailing_from_asymmetry(2.0), 1.5)


def test_resolution_base():
    assert math.isclose(calc_resolution_base(2.5, 4.0, 0.2, 0.25), 6.6667, rel_tol=1e-4)


def test_resolution_half():
    # 1.18 * (4.0-2.5) / (0.12+0.15) = 1.18*1.5/0.27 ≈ 6.5556
    assert math.isclose(calc_resolution_half(2.5, 4.0, 0.12, 0.15), 6.5556, rel_tol=1e-4)


# ---------- Input validation tests ----------
def test_validate_columns_missing():
    rows = [{"a": "1", "b": "2"}]
    with pytest.raises(SystemExit):
        validate_columns(rows, ["c"], "Test")


def test_to_float_invalid():
    with pytest.raises(SystemExit):
        to_float("abc", "col", 0)


def test_to_float_valid():
    assert to_float("3.14", "pi", 5) == 3.14


# ---------- Integration tests using temporary CSV files ----------
SAMPLE_CSV = """peak_name,retention_time,width_base,width_half,asymmetry
Compound A,2.5,0.2,0.12,1.1
Compound B,4.0,0.25,0.15,1.2
Solvent,1.0,0.1,0.06,1.0
"""

REPLICATE_CSV = """peak_name,retention_time,width_base,width_half,asymmetry,injection
Compound A,2.4,0.2,0.12,1.1,1
Compound A,2.6,0.2,0.12,1.1,2
Compound A,2.5,0.2,0.12,1.1,3
"""

INVALID_CSV = """peak_name,retention_time,width_base
Compound A,2.5,0.2
NotANumber,hello,0.2
"""


def _write_temp_csv(content: str) -> Path:
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, encoding='utf-8')
    tmp.write(content)
    tmp.close()
    return Path(tmp.name)


def test_basic_run_method_base():
    csv_file = _write_temp_csv(SAMPLE_CSV)
    try:
        rows = read_csv(str(csv_file))
        # default method base, tailing asymmetry, resolution base
        results = process_data(rows, method='base', tailing_method='asymmetry',
                               resolution_method='base', injection_col=None)
        assert len(results) == 3
        # Solvent first (tR=1.0), Compound A second, Compound B third
        assert results[0]['name'] == 'Solvent'
        assert math.isclose(results[0]['tr_mean'], 1.0)
        assert math.isclose(results[0]['plates'], calc_plates_base_width(1.0, 0.1))  # 1600
        assert math.isclose(results[0]['tailing'], calc_tailing_from_asymmetry(1.0))  # 1.0
        assert results[0]['resolution'] is not None  # resolution between Solvent and A
        assert math.isclose(results[0]['resolution'],
                            calc_resolution_base(1.0, 2.5, 0.1, 0.2), rel_tol=1e-4)
        assert math.isnan(results[0]['tr_rsd'])  # single injection

        # Compound A
        assert results[1]['name'] == 'Compound A'
        # resolution between A and B
        assert math.isclose(results[1]['resolution'],
                            calc_resolution_base(2.5, 4.0, 0.2, 0.25), rel_tol=1e-4)
        # Compound B has no next peak, resolution None
        assert results[2]['resolution'] is None
    finally:
        csv_file.unlink()


def test_method_half_height():
    csv_file = _write_temp_csv(SAMPLE_CSV)
    try:
        rows = read_csv(str(csv_file))
        results = process_data(rows, method='half', tailing_method='asymmetry',
                               resolution_method='half', injection_col=None)
        assert len(results) == 3
        # Check plates using half height
        assert math.isclose(results[0]['plates'],
                            calc_plates_half_height(1.0, 0.06), rel_tol=1e-4)
        # resolution half
        assert math.isclose(results[0]['resolution'],
                            calc_resolution_half(1.0, 2.5, 0.06, 0.12), rel_tol=1e-4)
    finally:
        csv_file.unlink()


def test_replicate_injections_rsd():
    csv_file = _write_temp_csv(REPLICATE_CSV)
    try:
        rows = read_csv(str(csv_file))
        results = process_data(rows, method='base', tailing_method='asymmetry',
                               resolution_method='base', injection_col='injection')
        assert len(results) == 1
        assert results[0]['name'] == 'Compound A'
        assert results[0]['tr_count'] == 3
        tr_vals = [2.4, 2.6, 2.5]
        mean = statistics.mean(tr_vals)
        rsd = (statistics.stdev(tr_vals) / mean) * 100.0
        assert math.isclose(results[0]['tr_mean'], mean, rel_tol=1e-4)
        assert math.isclose(results[0]['tr_rsd'], rsd, rel_tol=1e-4)
    finally:
        csv_file.unlink()


def test_invalid_numeric_data_exits():
    csv_file = _write_temp_csv(INVALID_CSV)
    try:
        rows = read_csv(str(csv_file))
        with pytest.raises(SystemExit):
            process_data(rows, method='base', tailing_method='asymmetry',
                         resolution_method='base', injection_col=None)
    finally:
        csv_file.unlink()


def test_missing_required_column_exits():
    rows = [{"peak_name": "A", "retention_time": "1.0"}]  # missing width_base
    with pytest.raises(SystemExit):
        process_data(rows, method='base', tailing_method='asymmetry',
                     resolution_method='base', injection_col=None)
