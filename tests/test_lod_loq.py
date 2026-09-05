"""Tests for LOD/LOQ calculator."""
import math
import csv
import pytest
import tempfile
from lab_tools.analytical.lod_loq import parse_csv, linear_regression, calculate_lod_loq


# -------------------------------------------------------------------
# calculate_lod_loq
# -------------------------------------------------------------------
def test_calculate_lod_loq_typical():
    """Basic LOD / LOQ from a typical calibration slope and Sy."""
    lod, loq = calculate_lod_loq(slope=0.0941, sy=0.0030)
    assert lod == pytest.approx(0.105, abs=0.001)
    assert loq == pytest.approx(0.319, abs=0.001)


def test_calculate_lod_loq_zero_slope_raises():
    with pytest.raises(ValueError, match="Slope must be positive"):
        calculate_lod_loq(0.0, 0.01)


def test_calculate_lod_loq_negative_slope_raises():
    with pytest.raises(ValueError, match="Slope must be positive"):
        calculate_lod_loq(-0.05, 0.01)


# -------------------------------------------------------------------
# linear_regression
# -------------------------------------------------------------------
def test_regression_perfect_line():
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 6, 8, 10]   # y = 2x
    reg = linear_regression(x, y)
    assert reg['slope'] == pytest.approx(2.0)
    assert reg['intercept'] == pytest.approx(0.0)
    assert reg['r_squared'] == pytest.approx(1.0)
    assert reg['sy'] == pytest.approx(0.0, abs=1e-9)


def test_regression_known_data():
    # from calibration example in README
    conc = [0, 2, 4, 6, 8, 10]
    signal = [0.002, 0.185, 0.377, 0.561, 0.751, 0.943]
    reg = linear_regression(conc, signal)
    assert reg['n'] == 6
    assert reg['slope'] == pytest.approx(0.0941, abs=0.001)
    assert reg['intercept'] == pytest.approx(-0.000667, abs=0.001)
    assert reg['r_squared'] >= 0.999
    assert reg['sy'] == pytest.approx(0.0030, abs=0.001)


def test_regression_insufficient_points():
    with pytest.raises(ValueError, match="Need >=3 points"):
        linear_regression([1, 2], [3, 4])


def test_regression_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        linear_regression([1, 2, 3], [4, 5])


def test_regression_constant_x():
    # all x values equal -> zero variance
    x = [5, 5, 5]
    y = [1, 2, 3]
    with pytest.raises(ValueError, match="near-zero variance"):
        linear_regression(x, y)


def test_regression_constant_y():
    # all y equal -> perfect horizontal fit, slope 0, Sy 0
    x = [1, 2, 3]
    y = [4, 4, 4]
    reg = linear_regression(x, y)
    assert reg['slope'] == pytest.approx(0.0)
    assert reg['intercept'] == pytest.approx(4.0)
    # r_squared: ss_total=0 leads to perfect fit? We set r_squared=1.0 if ss_total==0
    assert reg['r_squared'] == 1.0
    assert reg['sy'] == pytest.approx(0.0)


# -------------------------------------------------------------------
# parse_csv
# -------------------------------------------------------------------

SAMPLE_CSV = "conc,area\n0,0.002\n2,0.185\n4,0.377\n6,0.561\n8,0.751\n10,0.943\n"


def test_parse_valid_csv():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(SAMPLE_CSV)
        f.flush()
        conc, signal = parse_csv(f.name)
    assert conc == [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    assert signal == [0.002, 0.185, 0.377, 0.561, 0.751, 0.943]


def test_parse_no_header():
    # data without header, --no-header
    content = "0,0.002\n2,0.185\n4,0.377\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        f.flush()
        conc, signal = parse_csv(f.name, has_header=False)
    assert conc == [0.0, 2.0, 4.0]
    assert signal == [0.002, 0.185, 0.377]


def test_parse_custom_delimiter():
    content = "conc;area\n0;0.002\n2;0.185\n4;0.377\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        f.flush()
        conc, signal = parse_csv(f.name, delimiter=';')
    assert conc == [0.0, 2.0, 4.0]


def test_parse_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_csv("nonexistent_file.csv")


def test_parse_too_few_columns():
    content = "conc\n0\n2\n4\n"  # single column
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        f.flush()
        with pytest.raises(ValueError, match="expected 2 columns"):
            parse_csv(f.name)


def test_parse_non_numeric():
    content = "conc,area\n0,hello\n2,0.185\n4,0.377\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        f.flush()
        with pytest.raises(ValueError, match="non-numeric"):
            parse_csv(f.name)


def test_parse_insufficient_points():
    content = "conc,area\n0,0.002\n2,0.185\n"  # only 2 points
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(content)
        f.flush()
        with pytest.raises(ValueError, match="Insufficient calibration points"):
            parse_csv(f.name)


def test_parse_empty_file():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("")  # truly empty
        f.flush()
        # parsing an empty file with header expects at least one row
        with pytest.raises(ValueError, match="empty"):
            parse_csv(f.name)
