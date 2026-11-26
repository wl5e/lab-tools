"""Tests for the README per-tool index."""

from pathlib import Path

from lab_tools.cli import TOOLS

README = Path(__file__).resolve().parents[1] / "README.md"


def _index_rows():
    return [
        line
        for line in README.read_text(encoding="utf-8").splitlines()
        if line.startswith("| `")
    ]


def test_index_has_one_row_per_tool():
    rows = _index_rows()
    assert len(rows) == len(TOOLS) == 17


def test_index_lists_every_tool():
    listed = {row.split("`")[1] for row in _index_rows()}
    assert listed == set(TOOLS)


def test_index_rows_have_purpose_and_command():
    for row in _index_rows():
        assert "`lab-tools " in row
