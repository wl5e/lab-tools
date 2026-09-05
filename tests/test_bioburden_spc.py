import os
import tempfile
import pytest
from lab_tools.microbiology.bioburden_spc.chart import load_data, compute_individuals_limits, detect_violations, format_individual_chart

@pytest.fixture
def sample_csv():
    data = "date,bioburden\n2023-01-01,5\n2023-01-02,8\n2023-01-03,3\n2023-01-04,6\n2023-01-05,9\n2023-01-06,2\n2023-01-07,5\n2023-01-08,8\n2023-01-09,3\n2023-01-10,6\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, encoding='utf-8') as f:
        f.write(data)
        filename = f.name
    yield filename
    os.unlink(filename)

def test_load_data_basic(sample_csv):
    dates, values = load_data(sample_csv)
    assert len(dates) == 10
    assert len(values) == 10
    assert values[0] == 5.0

def test_load_data_missing_column():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("date,bad\n2023-01-01,5")
        fname = f.name
    try:
        with pytest.raises(ValueError, match="must contain.*bioburden"):
            load_data(fname, value_col='bioburden')
    finally:
        os.unlink(fname)

def test_load_data_invalid_number():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("date,bioburden\n2023-01-01,abc")
        fname = f.name
    try:
        with pytest.raises(ValueError, match="Invalid bioburden"):
            load_data(fname)
    finally:
        os.unlink(fname)

def test_compute_limits():
    values = [5, 8, 3]
    limits = compute_individuals_limits(values)
    assert limits['xbar'] == pytest.approx(5.3333333, rel=0.01)
    assert limits['avg_mr'] == pytest.approx(4.0)
    assert limits['sigma'] == pytest.approx(4.0/1.128, rel=0.01)
    assert limits['ucl_x'] == pytest.approx(5.3333 + 2.66*4, rel=0.01)
    assert limits['lcl_x'] == 0.0  # max(0,...)

def test_detect_violations_rule1():
    limits = {
        'xbar': 5.0, 'ucl_x': 6.0, 'lcl_x': 2.0, 'sigma': 1.0,
        'ucl1': 6.0, 'lcl1': 4.0, 'ucl2': 7.0, 'lcl2': 3.0,
        'avg_mr': 1.0, 'ucl_mr': 3.0, 'lcl_mr': 0.0, 'moving_ranges': [1]
    }
    dates = ['d1','d2','d3']
    values = [3, 10, 4]  # 10 > ucl
    v = detect_violations(dates, values, limits)
    assert any("Rule 1" in msg for _, msg in v)

def test_detect_violations_rule2():
    limits = {
        'xbar': 3.0, 'ucl_x': 10.0, 'lcl_x': 0.0, 'sigma': 1.0,
        'ucl1': 4.0, 'lcl1': 2.0, 'ucl2': 5.0, 'lcl2': 1.0,
        'avg_mr': 1.0, 'ucl_mr': 3.0, 'lcl_mr': 0.0, 'moving_ranges': [1,1]
    }
    dates = ['d1','d2','d3','d4','d5','d6','d7']
    values = [4, 5, 4.5, 6, 5, 4.8, 5.2]  # all > xbar (3)
    v = detect_violations(dates, values, limits)
    assert any("Seven points in a row above" in msg for _, msg in v)

def test_detect_violations_rule3():
    limits = {
        'xbar': 5.0, 'ucl_x': 10.0, 'lcl_x': 0.0, 'sigma': 1.0,
        'ucl1': 6.0, 'lcl1': 4.0, 'ucl2': 7.0, 'lcl2': 3.0,
        'avg_mr': 1.0, 'ucl_mr': 3.0, 'lcl_mr': 0.0, 'moving_ranges': [1,1]
    }
    dates = ['d1','d2','d3']
    values = [8, 7.5, 6]  # two >7 (ucl2) above mean
    v = detect_violations(dates, values, limits)
    assert any("Rule 3" in msg for _, msg in v)

def test_detect_violations_rule4():
    limits = {
        'xbar': 5.0, 'ucl_x': 10.0, 'lcl_x': 0.0, 'sigma': 1.0,
        'ucl1': 6.0, 'lcl1': 4.0, 'ucl2': 7.0, 'lcl2': 3.0,
        'avg_mr': 1.0, 'ucl_mr': 3.0, 'lcl_mr': 0.0, 'moving_ranges': [1,1]
    }
    dates = ['d1','d2','d3','d4','d5']
    values = [6.5, 7, 6.2, 6.8, 5.5]  # four of five > ucl1 and above mean
    v = detect_violations(dates, values, limits)
    assert any("Rule 4" in msg for _, msg in v)

def test_format_individual_chart():
    limits = {
        'xbar': 5.0, 'ucl_x': 8.0, 'lcl_x': 2.0, 'sigma': 1.0,
        'avg_mr': 1.0, 'ucl_mr': 3.0, 'lcl_mr': 0.0, 'moving_ranges': [1,1],
        'ucl1': 6.0, 'lcl1': 4.0, 'ucl2': 7.0, 'lcl2': 3.0
    }
    dates = ['d1','d2','d3']
    values = [5, 7, 3]
    violations = [(1, "Test violation")]
    chart = format_individual_chart(dates, values, limits, violations)
    assert '*' in chart
    assert 'X' in chart  # point 2 marked
    assert 'Mean: 5.00' in chart
