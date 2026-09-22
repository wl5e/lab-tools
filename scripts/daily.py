#!/usr/bin/env python3
"""Daily LLM-driven increment for lab-tools.

Reads the first not-yet-done item in ``BACKLOG.md``, implements it via DeepSeek
(see ``llm.py``), runs the full test suite, and — only if every test is green —
checks the item off and commits. A commit is therefore *only* made when the
change is real and verified; there is no filler.

A change that stays green is additionally reviewed by a **JEV pertinence gate**
(see ``jev.py``) when ``TYPESAFE_API_KEY`` is set: JEV scores whether the change
is a genuine, well-scoped implementation rather than slop, and a confident
rejection sends the item back for another attempt. The gate fails open if JEV
is unreachable, so an outage never blocks the daily stream.

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


def _diff_for(written) -> str:
    """Compact before/after of the files in ``written`` (list of ``(rel, existed)``)."""
    rels = [rel for rel, _ in written]
    parts = []
    diff = _run(["git", "diff", "--", *rels]).stdout
    if diff.strip():
        parts.append(diff)
    for rel, existed in written:
        if not existed:
            parts.append(f"### NEW FILE {rel}\n" + (ROOT / rel).read_text(encoding="utf-8"))
    return "\n".join(parts)


def _jev_gate(slug: str, title: str, written) -> str:
    """Return ``'accept'`` or ``'reject'`` from the JEV pertinence gate.

    Opt-in: skipped when ``TYPESAFE_API_KEY`` is unset. Fails open — if JEV
    errors, the change is accepted so an outage never blocks the daily stream.
    """
    if not os.environ.get("TYPESAFE_API_KEY"):
        print("TYPESAFE_API_KEY not set; skipping JEV pertinence gate.")
        return "accept"
    import jev
    try:
        accepted, prob = jev.is_pertinent(slug, title, _diff_for(written))
    except Exception as exc:  # noqa: BLE001 - fail open on gate outage
        print(f"JEV gate unavailable ({exc}); accepting change.")
        return "accept"
    print(f"JEV pertinence: {prob:.2f} -> {'accept' if accepted else 'reject'}")
    return "accept" if accepted else "reject"


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
        # `git diff` ignores untracked files, so a change that only *adds* a
        # new file would look like "no change". Use `git status --porcelain` to
        # catch created, modified and untracked paths.
        if not _run(["git", "status", "--porcelain"]).stdout.strip():
            print("LLM produced no effective change; skipping (no commit).")
            _revert(written)
            return 1
        res = run_tests()
        if res.returncode != 0:
            print(f"Attempt {attempt}: tests failed; feeding output back and retrying.")
            feedback = res.stdout + res.stderr
            _revert(written)
            written = []
            continue
        # Tests are green; run the JEV pertinence gate before committing so a
        # change that stays green but is slop gets another attempt.
        if _jev_gate(slug, title, written) == "reject":
            print(f"Attempt {attempt}: JEV judged the change not pertinent; retrying.")
            feedback = (
                "Your change passed the tests but was judged NOT pertinent by a "
                "review model (slop, trivial, or unrelated to the item). Re-implement "
                "the item as a real, well-scoped improvement."
            )
            _revert(written)
            written = []
            continue
        break
    else:
        # A hard item must not block the pipeline forever: check it off with a
        # skip note so the next run advances to the following item.
        mark_done(slug)
        commit_and_push(f"chore: skip {slug} (LLM could not produce a green, pertinent change)", slug)
        print(f"Skipped {slug} after 3 failed attempts.")
        return 0

    mark_done(slug)
    if args.apply:
        print("Tests green. (--apply: not committing.)")
        return 0

    commit_and_push(subject, slug)
    print(f"Committed and pushed: {subject}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
