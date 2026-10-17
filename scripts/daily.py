#!/usr/bin/env python3
"""Daily LLM-driven increment for lab-tools.

Reads the first open item in ``BACKLOG.md``, implements it (deterministic
handler or DeepSeek), and commits only when the change is real and verified:

- at least one file changed;
- the local scope gate passes (<= 3 files, a test when production code changes,
  <= 150 changed lines);
- the full test suite is green.

An item that fails three attempts is moved to ``BACKLOG.blocked.md`` and the
run tries the next item (up to ``MAX_ITEMS`` per run). There is deliberately no
"skip" commit: the repository only records work that actually shipped.

Modes:
    --dry-run   print the next items without touching anything
    --apply     apply the change and run tests, but do not commit (local check)
    (default)   apply + test + commit + push (used by GitHub Actions)
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BACKLOG = ROOT / "BACKLOG.md"
BLOCKED = ROOT / "BACKLOG.blocked.md"

MAX_ITEMS = 3
MAX_FILES = 3
MAX_DIFF_LINES = 150

_ITEM_RE = re.compile(r"^- \[ \] `([^`]+)` (.*?)\s*\|\s*handler:\s*(\w*)\s*$")
_IGNORED_PREFIXES = ("BACKLOG", ".pytest_cache", "__pycache__")

# Deterministic handlers, registered by the name used in ``| handler: ...``.
# Handlers take no arguments, edit the working tree and return a commit subject.


def _run(cmd, **kwargs):
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, **kwargs)


def _replace(path: Path, old: str, new: str, count: int = 1) -> None:
    """Replace ``old`` with ``new`` in ``path``; fail loudly on anchor drift."""
    text = path.read_text(encoding="utf-8")
    found = text.count(old)
    if found != count:
        raise RuntimeError(
            f"{path.name}: expected {count} occurrence(s) of anchor, found {found}"
        )
    path.write_text(text.replace(old, new), encoding="utf-8")


def hemocytometer_chamber_depth() -> str:
    """Expose the existing chamber factor through the hemocytometer CLI."""
    module = ROOT / "lab_tools" / "lab_ops" / "hemocytometer.py"

    flags = (
        "    {parser}.add_argument('--chamber-depth', type=float, default=0.1, "
        "help='Chamber depth in mm (default 0.1)')\n"
        "    {parser}.add_argument('--chamber-factor', type=float, default=None, "
        "help='Chamber factor in cells/mL; overrides --chamber-depth')\n"
    )
    _replace(
        module,
        "    count_parser.add_argument('--dilution', type=float, default=1.0, "
        "help='Dilution factor (default 1.0)')\n",
        "    count_parser.add_argument('--dilution', type=float, default=1.0, "
        "help='Dilution factor (default 1.0)')\n"
        + flags.format(parser="count_parser"),
    )
    _replace(
        module,
        "    dil_parser.add_argument('--dilution', type=float, default=1.0, "
        "help='Dilution factor of counted sample')\n",
        "    dil_parser.add_argument('--dilution', type=float, default=1.0, "
        "help='Dilution factor of counted sample')\n"
        + flags.format(parser="dil_parser"),
    )
    _replace(
        module,
        "    vol_parser.add_argument('--dilution', type=float, default=1.0, "
        "help='Dilution factor')\n",
        "    vol_parser.add_argument('--dilution', type=float, default=1.0, "
        "help='Dilution factor')\n"
        + flags.format(parser="vol_parser"),
    )

    _replace(
        module,
        "def main(argv: Optional[List[str]] = None) -> int:",
        "def _chamber_factor(args) -> float:\n"
        "    if args.chamber_factor is not None:\n"
        "        return args.chamber_factor\n"
        "    return args.chamber_depth * 1e5\n"
        "\n"
        "\n"
        "def main(argv: Optional[List[str]] = None) -> int:",
    )
    _replace(
        module,
        "            cc = CellCounter(live, dead, squares, dilution)",
        "            cc = CellCounter(live, dead, squares, dilution, _chamber_factor(args))",
    )
    _replace(
        module,
        "            cc = CellCounter(args.live, args.dead, args.squares, args.dilution)",
        "            cc = CellCounter(args.live, args.dead, args.squares, args.dilution, _chamber_factor(args))",
        count=2,
    )

    tests = ROOT / "tests" / "test_hemocytometer.py"
    tests.write_text(
        tests.read_text(encoding="utf-8")
        + '''

def test_chamber_depth_scales_cells_per_ml():
    shallow = CellCounter(100, 0, 4, 1.0, 1e4).cells_per_ml
    deep = CellCounter(100, 0, 4, 1.0, 2e4).cells_per_ml
    assert deep == 2 * shallow


def test_chamber_depth_cli_flag(capsys):
    from lab_tools.lab_ops.hemocytometer import main

    assert main(["count", "--live", "100", "--squares", "4", "--chamber-depth", "0.2"]) == 0
    out = capsys.readouterr().out
    assert "500,000" in out


def test_chamber_factor_cli_flag(capsys):
    from lab_tools.lab_ops.hemocytometer import main

    assert main(["count", "--live", "100", "--squares", "4", "--chamber-factor", "20000"]) == 0
    out = capsys.readouterr().out
    assert "500,000" in out
''',
        encoding="utf-8",
    )
    return "feat(hemocytometer): expose chamber depth and factor in the CLI"


HANDLERS = {
    "hemocytometer_chamber_depth": hemocytometer_chamber_depth,
}


def open_items():
    items = []
    for line in BACKLOG.read_text(encoding="utf-8").splitlines():
        m = _ITEM_RE.match(line)
        if m:
            items.append((m.group(1), m.group(2), m.group(3)))
    return items


def next_item():
    items = open_items()
    return items[0] if items else None


def mark_done(slug: str) -> None:
    lines = BACKLOG.read_text(encoding="utf-8").splitlines()
    for i, line in enumerate(lines):
        m = _ITEM_RE.match(line)
        if m and m.group(1) == slug:
            lines[i] = line.replace("- [ ]", "- [x]", 1)
            break
    BACKLOG.write_text("\n".join(lines) + "\n", encoding="utf-8")


def block_item(slug: str, title: str, reason: str) -> None:
    """Move an item out of the active queue, with an audit trail."""
    lines = BACKLOG.read_text(encoding="utf-8").splitlines()
    kept = []
    for line in lines:
        m = _ITEM_RE.match(line)
        if m and m.group(1) == slug:
            continue
        kept.append(line)
    BACKLOG.write_text("\n".join(kept).rstrip() + "\n", encoding="utf-8")

    if not BLOCKED.exists():
        BLOCKED.write_text(
            "# Blocked items\n\n"
            "Moved here after three failed attempts. Fix the item (or split it)\n"
            "and move it back to BACKLOG.md to retry.\n\n",
            encoding="utf-8",
        )
    with BLOCKED.open("a", encoding="utf-8") as handle:
        handle.write(
            f"- `{slug}` {title}\n"
            f"  blocked: {reason} ({date.today().isoformat()})\n"
        )


def run_tests():
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    return _run([sys.executable, "-m", "pytest", "-q"], env=env)


def _revert(written) -> None:
    """Undo only the files this attempt wrote (never the whole tree)."""
    for rel, existed in written:
        p = ROOT / rel
        if existed:
            _run(["git", "checkout", "--", rel])
        elif p.exists():
            p.unlink()


def _changed_paths():
    tracked = _run(["git", "diff", "--name-only"]).stdout.splitlines()
    untracked = _run(["git", "ls-files", "--others", "--exclude-standard"]).stdout.splitlines()
    paths = [p for p in dict.fromkeys(tracked + untracked) if p]
    return [p for p in paths if not p.startswith(_IGNORED_PREFIXES)]


def _diff_lines() -> int:
    total = 0
    for line in _run(["git", "diff", "--numstat"]).stdout.splitlines():
        parts = line.split("\t")
        for value in parts[:2]:
            if value.isdigit():
                total += int(value)
    for rel in _run(["git", "ls-files", "--others", "--exclude-standard"]).stdout.splitlines():
        if rel.startswith(_IGNORED_PREFIXES):
            continue
        path = ROOT / rel
        if path.is_file():
            total += len(path.read_text(encoding="utf-8", errors="ignore").splitlines())
    return total


def _local_gate(written) -> tuple[bool, str]:
    """Free, deterministic scope gate."""
    paths = [rel for rel, _ in written] if written else _changed_paths()
    if not paths:
        return False, "no file changed"
    if len(paths) > MAX_FILES:
        return False, f"{len(paths)} files changed (max {MAX_FILES})"
    has_test = any(rel.startswith("tests/") for rel in paths)
    has_production = any(
        rel.startswith("lab_tools/") and rel.endswith(".py") for rel in paths
    )
    if has_production and not has_test:
        return False, "production code changed without a test"
    changed = _diff_lines()
    if changed > MAX_DIFF_LINES:
        return False, f"{changed} changed lines (max {MAX_DIFF_LINES})"
    return True, ""


def _commit_and_push(subject: str, body: str, add_args) -> None:
    env = {**os.environ}
    name = env.get("GIT_AUTHOR_NAME", "Collins Amatu Gorgerat")
    email = env.get("GIT_AUTHOR_EMAIL", "133529715+wl5e@users.noreply.github.com")
    _run(["git", "config", "user.name", name])
    _run(["git", "config", "user.email", email])
    _run(["git", "add", *add_args])
    res = _run(["git", "commit", "-m", subject, "-m", body])
    if res.returncode != 0:
        raise RuntimeError(f"git commit failed: {res.stderr.strip()}")
    push = _run(["git", "push"])
    if push.returncode != 0:
        print(f"git push failed: {push.stderr.strip()}; rebasing and retrying.")
        rebase = _run(["git", "pull", "--rebase"])
        if rebase.returncode != 0:
            raise RuntimeError(f"git pull --rebase failed: {rebase.stderr.strip()}")
        push = _run(["git", "push"])
    if push.returncode != 0:
        raise RuntimeError(f"git push failed after retry: {push.stderr.strip()}")


def commit_and_push(subject: str, slug: str) -> None:
    body = f"Backlog item: {slug}\n\nAutomated daily increment — test suite green."
    _commit_and_push(subject, body, ["-A"])


def commit_blocked(slugs) -> None:
    body = "Blocked items moved to BACKLOG.blocked.md.\n\n" + "\n".join(
        f"- {slug}" for slug in slugs
    )
    _commit_and_push(
        f"chore: block {len(slugs)} backlog item(s) after 3 failed attempts",
        body,
        ["BACKLOG.md", "BACKLOG.blocked.md"],
    )


def _attempt(slug: str, title: str, handler: str, apply_only: bool) -> tuple[bool, str]:
    """Try one item. Returns (shipped, reason_if_not)."""
    if handler:
        fn = HANDLERS.get(handler)
        if fn is None:
            return False, f"unknown handler {handler!r}"
        try:
            subject = fn()
        except Exception as exc:  # noqa: BLE001 - handlers are trusted but may drift
            _revert([])
            return False, f"handler failed: {exc}"
        if not _changed_paths():
            return False, "handler produced no change"
        ok, reason = _local_gate([])
        if not ok:
            _revert([])
            return False, reason
        if run_tests().returncode != 0:
            _revert([])
            return False, "tests failed"
        mark_done(slug)
        if not apply_only:
            commit_and_push(subject, slug)
        return True, ""

    import llm

    written = []
    feedback = None
    for attempt in (1, 2, 3):
        try:
            subject, written = llm.implement(slug, title, feedback)
        except Exception as exc:  # noqa: BLE001 - transient provider/JSON failures
            print(f"LLM attempt {attempt} raised: {exc}")
            _revert(written)
            written = []
            feedback = (
                f"Your previous response could not be used: {exc}. "
                "Return ONLY a valid JSON object with the full content of each "
                "changed file, and keep the change small."
            )
            continue

        if not written:
            print(f"Attempt {attempt}: the model changed nothing; retrying.")
            feedback = (
                "You must change at least one file. There is no 'already satisfied' "
                "answer: implement the item for real, with a test."
            )
            _revert(written)
            written = []
            continue

        workflow_touched = [rel for rel, _ in written if rel.startswith(".github/workflows/")]
        if workflow_touched:
            print(
                f"Attempt {attempt}: change touches {workflow_touched}; "
                "reverting and retrying (GITHUB_TOKEN cannot push workflow files)."
            )
            feedback = (
                "Do NOT modify or create any file under .github/workflows/ — the "
                "automation's push token is rejected by GitHub for workflow files. "
                "Re-implement the item without touching the workflows directory."
            )
            _revert(written)
            written = []
            continue

        if not _changed_paths():
            print("LLM produced no effective change; retrying.")
            feedback = "Your files were identical to the current content. Make a real change."
            _revert(written)
            written = []
            continue

        ok, reason = _local_gate(written)
        if not ok:
            print(f"Attempt {attempt}: local gate rejected the change: {reason}")
            feedback = (
                f"Your change was rejected by the local scope gate: {reason}. "
                "Keep it to at most 3 files, add or update a test when you change "
                "production code, and keep the diff under 150 lines."
            )
            _revert(written)
            written = []
            continue

        res = run_tests()
        if res.returncode != 0:
            print(f"Attempt {attempt}: tests failed; feeding output back and retrying.")
            feedback = res.stdout + res.stderr
            _revert(written)
            written = []
            continue

        mark_done(slug)
        if not apply_only:
            commit_and_push(subject, slug)
            print(f"Committed and pushed: {subject}")
        else:
            print("Tests green. (--apply: not committing.)")
        return True, ""

    return False, "3 attempts failed"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="print the plan only")
    parser.add_argument("--apply", action="store_true", help="apply + test, no commit")
    args = parser.parse_args(argv)

    items = open_items()
    if not items:
        print("No undone backlog item. Nothing to do.")
        return 0

    if args.dry_run:
        for slug, _, handler in items[:MAX_ITEMS]:
            mode = f"handler:{handler}" if handler else "LLM"
            print(f"- {slug} [{mode}]")
        return 0

    if not os.environ.get("DEEPSEEK_API_KEY") and not any(h for _, _, h in items[:MAX_ITEMS]):
        print("DEEPSEEK_API_KEY not set; nothing to do.")
        return 0

    blocked = []
    for slug, title, handler in items[:MAX_ITEMS]:
        print(f"Selected backlog item: {slug}")
        shipped, reason = _attempt(slug, title, handler, args.apply)
        if shipped:
            return 0
        print(f"Blocking {slug}: {reason}")
        block_item(slug, title, reason)
        blocked.append(slug)

    if blocked and not args.apply:
        commit_blocked(blocked)
    print(f"No item shipped this run; blocked: {', '.join(blocked)}.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
