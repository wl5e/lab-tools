#!/usr/bin/env python3
"""Daily LLM-driven increment for lab-tools.

Reads the first not-yet-done item in ``BACKLOG.md``, implements it via DeepSeek
(see ``llm.py``), runs the full test suite, and — only if every test is green —
checks the item off and commits. A commit is therefore *only* made when the
change is real and verified; there is no filler.

Modes:
    --dry-run   print the plan without touching anything
    --apply     apply the change and run tests, but do not commit (local check)
    (default)   apply + test + commit + push (used by GitHub Actions)
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BACKLOG = ROOT / "BACKLOG.md"

_ITEM_RE = re.compile(r"^- \[ \] `([^`]+)` (.*?)\s*\|\s*handler:\s*(\w*)\s*$")


def _run(cmd, **kwargs):
    return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, **kwargs)


def next_item():
    for line in BACKLOG.read_text(encoding="utf-8").splitlines():
        m = _ITEM_RE.match(line)
        if m:
            return m.group(1), m.group(2), m.group(3)
    return None


def mark_done(slug: str) -> None:
    lines = BACKLOG.read_text(encoding="utf-8").splitlines()
    for i, line in enumerate(lines):
        m = _ITEM_RE.match(line)
        if m and m.group(1) == slug:
            lines[i] = line.replace("- [ ]", "- [x]", 1)
            break
    BACKLOG.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_tests():
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    return _run([sys.executable, "-m", "pytest", "-q"], env=env)


def _revert(written) -> None:
    if written:
        for rel, existed in written:
            p = ROOT / rel
            if existed:
                _run(["git", "checkout", "--", rel])
            elif p.exists():
                p.unlink()
    else:
        _run(["git", "checkout", "--", "."])


def commit_and_push(subject: str, slug: str) -> None:
    env = {**os.environ}
    name = env.get("GIT_AUTHOR_NAME", "Collins Amatu Gorgerat")
    email = env.get("GIT_AUTHOR_EMAIL", "133529715+wl5e@users.noreply.github.com")
    _run(["git", "config", "user.name", name])
    _run(["git", "config", "user.email", email])
    _run(["git", "add", "-A"])
    body = f"Backlog item: {slug}\n\nAutomated daily increment — test suite green."
    res = _run(["git", "commit", "-m", subject, "-m", body])
    if res.returncode != 0:
        raise RuntimeError(f"git commit failed: {res.stderr.strip()}")
    _run(["git", "push"])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="print the plan only")
    parser.add_argument("--apply", action="store_true", help="apply + test, no commit")
    args = parser.parse_args(argv)

    item = next_item()
    if item is None:
        print("No undone backlog item. Nothing to do.")
        return 0

    slug, title, _ = item
    print(f"Selected backlog item: {slug}")

    if args.dry_run:
        print("Dry run — no changes made.")
        return 0

    if not os.environ.get("DEEPSEEK_API_KEY"):
        print("DEEPSEEK_API_KEY not set; nothing to do.")
        return 0

    import llm
    written = []
    feedback = None
    for attempt in (1, 2, 3):
        try:
            subject, written = llm.implement(slug, title, feedback)
        except Exception as exc:  # noqa: BLE001 - fail the day cleanly
            print(f"LLM attempt {attempt} raised: {exc}")
            _revert(written)
            return 1
        if not written:
            # The item is already satisfied — check it off, no code change.
            mark_done(slug)
            commit_and_push(f"chore: {slug} already satisfied (no change needed)", slug)
            print(f"Marked {slug} done (already satisfied).")
            return 0
        if not _run(["git", "diff", "--stat"]).stdout.strip():
            print("LLM produced no effective change; skipping (no commit).")
            _revert(written)
            return 1
        res = run_tests()
        if res.returncode == 0:
            break
        print(f"Attempt {attempt}: tests failed; feeding output back and retrying.")
        feedback = res.stdout + res.stderr
        _revert(written)
        written = []
    else:
        print("LLM did not produce a green change after 3 attempts. No commit.")
        return 1

    mark_done(slug)
    if args.apply:
        print("Tests green. (--apply: not committing.)")
        return 0

    commit_and_push(subject, slug)
    print(f"Committed and pushed: {subject}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
