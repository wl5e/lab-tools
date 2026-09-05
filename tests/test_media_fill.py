import pytest
import tempfile
import os
import lab_tools.microbiology.media_fill as mfa

def test_wilson_ci_normal():
    rate, lo, hi = mfa.wilson_ci(100, 5, 1.96)
    assert 0.01 <= lo <= rate <= hi <= 0.15
    assert rate == 0.05

def test_wilson_ci_zero_success():
    rate, lo, hi = mfa.wilson_ci(50, 0, 1.96)
    assert rate == 0.0
    assert lo == 0.0
    assert hi > 0.0

def test_wilson_ci_all_success():
    rate, lo, hi = mfa.wilson_ci(50, 50, 1.96)
    assert rate == 1.0
    assert hi == 1.0
    assert lo < 1.0

def test_get_z_valid():
    assert mfa.get_z(0.05) == 1.96

def test_get_z_invalid():
    with pytest.raises(ValueError):
        mfa.get_z(0.03)

def test_parse_csv_valid():
    csv_content = "batch,n,contaminated\nB1,100,2\nB2,150,0\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        fname = f.name
    rows = mfa.parse_csv(fname)
    os.unlink(fname)
    assert len(rows) == 2
    assert rows[0]['batch'] == 'B1'
    assert rows[0]['n'] == 100
    assert rows[0]['x'] == 2

def test_parse_csv_missing_column():
    csv_content = "batch,contaminated\nB1,2\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        fname = f.name
    with pytest.raises(ValueError, match="Missing required"):
        mfa.parse_csv(fname)
    os.unlink(fname)

def test_parse_csv_invalid_number():
    csv_content = "batch,n,contaminated\nB1,100,not_a_number\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        fname = f.name
    with pytest.raises(ValueError):
        mfa.parse_csv(fname)
    os.unlink(fname)

def test_cochran_armitage_no_trend():
    res = [
        mfa.BatchResult('B1', 100, 5, 0.05, 0.02, 0.1),
        mfa.BatchResult('B2', 100, 5, 0.05, 0.02, 0.1),
        mfa.BatchResult('B3', 100, 5, 0.05, 0.02, 0.1)
    ]
    Z, p = mfa.cochran_armitage_trend(res)
    assert abs(Z) < 0.01
    assert p > 0.4

def test_cochran_armitage_increasing():
    res = [
        mfa.BatchResult('B1', 100, 2, 0.02, 0.005, 0.07),
        mfa.BatchResult('B2', 100, 4, 0.04, 0.01, 0.10),
        mfa.BatchResult('B3', 100, 8, 0.08, 0.03, 0.16)
    ]
    Z, p = mfa.cochran_armitage_trend(res)
    assert Z > 1.5
    assert p < 0.1

def test_aggregate():
    res = [
        mfa.BatchResult('B1', 200, 4, 0.02, 0.01, 0.04),
        mfa.BatchResult('B2', 300, 6, 0.02, 0.01, 0.04)
    ]
    agg = mfa.aggregate(res, 0.05)
    assert agg.n == 500
    assert agg.x == 10
    assert agg.rate == 0.02
    assert 0.01 <= agg.ci_lower <= agg.rate <= agg.ci_upper <= 0.04

def test_full_analysis_csv(tmp_path):
    # Create a simple CSV file and run compute_analysis to check output file
    csv_file = tmp_path / "test.csv"
    csv_file.write_text("batch,n,contaminated\nB1,1000,5\nB2,1000,8\n")
    out_file = tmp_path / "out.csv"
    summary = mfa.compute_analysis(str(csv_file), str(out_file), 0.05)
    assert "Aggregate" in summary
    assert out_file.exists()
    with open(out_file) as f:
        lines = f.readlines()
    assert len(lines) >= 4  # header + 2 batch rows + TOTAL + blank + trend
