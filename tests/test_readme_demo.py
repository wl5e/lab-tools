"""Ensure every tool has a demo command in the README."""

from pathlib import Path

from lab_tools.cli import TOOLS, main


def test_list_command_prints_tools(capsys):
    assert main(["list"]) == 0
    out = capsys.readouterr().out
    listed = set(out.strip().splitlines())
    assert listed == set(TOOLS)


def test_every_tool_has_demo_in_readme():
    readme = (Path(__file__).parent.parent / "README.md").read_text()
    for tool in TOOLS:
        assert f"lab-tools {tool}" in readme, f"missing demo for {tool}"
