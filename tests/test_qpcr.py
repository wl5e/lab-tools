"""Tests for the qPCR ΔΔCq analyzer module."""

import csv
import io
import math
import tempfile
import os
import sys
import pytest
from lab_tools.molecular.qpcr import (
    _read_cq_data,
    _compute_means_and_sds,
    _find_control_sample,
    analyze_qpcr,
    export_csv,
    Result
)

SIMPLE_CSV = """Sample,Gene,Cq
WT,GAPDH,18.5
WT,GAPDH,18.7
WT,Target,22.1
WT,Target,21.9
KO,GAPDH,19.0
KO,GAPDH,18.8
KO,Target,20.5
KO,Target,20.3"""

def _write_temp_csv(content: str, **write_kwargs) -> str:
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False, **write_kwargs) as f:
        f.write(content)
        f.flush()
        return f.name


def test_read_cq_data_basic():
    path = _write_temp_csv(SIMPLE_CSV)
    try:
        data = _read_cq_data(path, ',')
        assert set(data.keys()) == {'WT', 'KO'}
        assert data['WT']['GAPDH'] == [18.5, 18.7]
        assert data['KO']['Target'] == [20.5, 20.3]
    finally:
        os.unlink(path)


def test_read_cq_data_missing_column():
    csv = "Sample,Dummy,Cq\nA,X,20"
    path = _write_temp_csv(csv)
    try:
        with pytest.raises(ValueError, match='Missing required'):
            _read_cq_data(path, ',')
    finally:
        os.unlink(path)


def test_read_cq_data_non_numeric():
    csv = "Sample,Gene,Cq\nA,G1,NA\nA,G1,20.0"
    path = _write_temp_csv(csv)
    try:
        data = _read_cq_data(path, ',')
        assert data['A']['G1'] == [20.0]  # NA skipped
    finally:
        os.unlink(path)


def test_compute_means_and_sds():
    raw = {'WT': {'GAPDH': [18.5, 18.7], 'Target': [22.1]}}
    summary = _compute_means_and_sds(raw)
    mean_gap, sd_gap = summary['WT']['GAPDH']
    assert round(mean_gap, 4) == 18.6
    assert round(sd_gap, 4) == 0.1414  # sample std
    mean_target, sd_target = summary['WT']['Target']
    assert mean_target == 22.1
    assert sd_target == 0.0


def test_find_control_sample_unique():
    samples = ['WT_untreated', 'KO', 'WT_treated']
    ctrl = _find_control_sample(samples, 'untreated')
    assert ctrl == 'WT_untreated'


def test_find_control_sample_no_match():
    with pytest.raises(ValueError, match='No sample matches'):
        _find_control_sample(['A', 'B'], 'xyz')


def test_find_control_sample_multiple_matches():
    with pytest.raises(ValueError, match='Multiple samples'):
        _find_control_sample(['WT_1', 'WT_2', 'KO'], 'WT')


def test_analyze_qpcr_full():
    path = _write_temp_csv(SIMPLE_CSV)
    try:
        results = analyze_qpcr(path, 'GAPDH', 'WT', ',')
    finally:
        os.unlink(path)

    # One non-control sample KO, gene Target
    assert len(results) == 1
    r = results[0]
    assert r.sample == 'KO'
    assert r.gene == 'Target'
    # mean Cq KO Target: 20.4, GAPDH: 18.9 -> dCq = 1.5
    assert round(r.mean_cq_gene, 2) == 20.40
    assert round(r.mean_cq_ref, 2) == 18.90
    assert round(r.dCq, 2) == 1.50
    # WT dCq: Target mean 22.0, GAPDH 18.6 -> 3.4
    # ddCq = 1.5 - 3.4 = -1.9
    assert round(r.ddCq, 2) == -1.90
    # fold change = 2^(1.9) ~ 3.7321
    assert round(r.fold_change, 4) == 3.7321
    # error: SD KO Target = 0.1414, KO GAPDH=0.1414, WT Target=0.1414, WT GAPDH=0.1414
    # dCq KO SD = sqrt(0.02+0.02)=0.2, dCq WT SD = sqrt(0.02+0.02)=0.2
    # ddCq SD = sqrt(0.04+0.04)=0.2828
    assert round(r.sd_ddCq, 4) == 0.2828


def test_missing_ref_gene_in_sample():
    csv = "Sample,Gene,Cq\nCtrl,GAPDH,20\nCtrl,GOI,25\nTest,GOI,22"  # Test missing GAPDH
    path = _write_temp_csv(csv)
    try:
        results = analyze_qpcr(path, 'GAPDH', 'Ctrl', ',')
        # Only Ctrl analysed; Test skipped due to missing ref
        assert len(results) == 0
    finally:
        os.unlink(path)


def test_control_not_found():
    csv = "Sample,Gene,Cq\nA,GAPDH,20"
    path = _write_temp_csv(csv)
    try:
        with pytest.raises(ValueError, match='No sample matches'):
            analyze_qpcr(path, 'GAPDH', 'NonExistent', ',')
    finally:
        os.unlink(path)


def test_export_csv_string_output():
    results = [
        Result('KO', 'Target', 20.4, 0.14, 18.9, 0.14, 1.5, 0.2, -1.9, 0.283, 3.7325, 3.0, 4.5)
    ]
    f = io.StringIO()
    export_csv(results, f, delimiter=',')
    output = f.getvalue()
    reader = csv.DictReader(io.StringIO(output))
    rows = list(reader)
    assert len(rows) == 1
    assert rows[0]['Sample'] == 'KO'
    assert rows[0]['Gene'] == 'Target'
    # numerical strings
    assert float(rows[0]['fold_change']) == 3.7325
    assert float(rows[0]['fold_change_low']) == 3.0
    assert float(rows[0]['fold_change_high']) == 4.5
