"""DeepSeek-backed implementation of backlog items that lack a handler.

Used by ``daily.py`` when the next undone item has no deterministic handler and
``DEEPSEEK_API_KEY`` is set. The model is given the item description plus a
*scoped* view of the repository — the manifest (package structure) and only the
files of the tool(s) named by the item — and returns the full new content of
only the files it changes. The change is committed only if the test suite
stays green.

Scoping is what keeps cost and latency down: sending one tool's files instead
of the whole 17-tool codebase. The weekly replenish (``propose_items``) keeps
the full view because it needs to see everything to avoid duplicates.

The API is OpenAI-compatible, so swapping to another provider later is a
one-line ``API_URL`` / model change.
"""

from __future__ import annotations

import json
import os
import re
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(__file__).resolve().parent.parent
API_URL = "https://api.deepseek.com/chat/completions"
DEFAULT_MODEL = "deepseek-v4-flash"

# Always sent: a lightweight manifest so the model knows the package layout and
# which tools exist, without the cost of every source file.
MANIFEST_GLOBS = [
    "lab_tools/__init__.py",
    "lab_tools/cli.py",
    "lab_tools/*/__init__.py",
    "pyproject.toml",
]

SYSTEM_PROMPT = (
    "You are an expert Python developer maintaining `lab-tools`, a library of "
    "small, tested CLI tools for pharmaceutical / microbiology / molecular-"
    "biology laboratory data (pure Python or numpy/scipy, no other deps). "
    "Implement exactly ONE well-scoped backlog item.\n"
    "Rules:\n"
    "- Make a MINIMAL, correct change; touch 1-3 files at most.\n"
    "- If you change production code, add or update a test that covers it.\n"
    "- Do NOT break existing behaviour; the existing test suite must keep passing.\n"
    "- Follow the existing style and keep the same public API.\n"
    "- Return ONLY a JSON object, with this exact shape and no surrounding prose:\n"
    '  {"summary": "<conventional commit subject>",\n'
    '   "files": {"<relative/path>": "<full new file content>", ...}}\n'
    '- "files" holds the FULL new content of every file you create or change.\n'
    "- If the item is already satisfied, return an empty files map and say so in summary."
)

PROPOSE_PROMPT = (
    "You are the maintainer of `lab-tools`, a library of small CLI tools for "
    "pharmaceutical / microbiology / molecular-biology laboratory data. "
    "Propose NEW backlog items that would genuinely improve it.\n"
    "Requirements:\n"
    "- Propose 10 items, each concrete and implementable in ONE focused change.\n"
    "- Do NOT duplicate anything already in the backlog or already in the code.\n"
    "- No trivial filler (comment edits, renames, cosmetic tweaks).\n"
    "- Slugs are kebab-case.\n"
    '- Return ONLY JSON: {"items": [{"slug": "...", "title": "...", "why": "..."}]}'
)


def _index_tools() -> Dict[str, List[Path]]:
    """Map each tool name to its source files (module + its tests)."""
    tools: Dict[str, List[Path]] = {}
    for domain in (ROOT / "lab_tools").iterdir():
        if not domain.is_dir() or domain.name.startswith("_"):
            continue
        for entry in sorted(domain.iterdir()):
            if entry.name in ("__init__.py", "__pycache__"):
                continue
            if entry.is_file() and entry.suffix == ".py":
                tools.setdefault(entry.stem, []).append(entry)
            elif entry.is_dir():
                tools.setdefault(entry.name, []).extend(sorted(entry.rglob("*.py")))
    tests = ROOT / "tests"
    if tests.exists():
        for tool, paths in tools.items():
            for tp in sorted(tests.glob(f"test_{tool}*.py")):
                if tp not in paths:
                    paths.append(tp)
    return tools


def _manifest_paths() -> List[Path]:
    paths: List[Path] = []
    for g in MANIFEST_GLOBS:
        for p in sorted(ROOT.glob(g)):
            if p.is_file():
                paths.append(p)
    return paths


def _render(paths: List[Path]) -> str:
    return "\n\n".join(
        f"### {p.relative_to(ROOT)}\n{p.read_text(encoding='utf-8')}" for p in paths
    )


def _gather_context_full() -> str:
    """The whole codebase — used only by the weekly replenish."""
    paths: List[Path] = []
    for pattern in ["lab_tools/**/*.py", "tests/**/*.py"]:
        paths.extend(p for p in sorted(ROOT.glob(pattern)) if "__pycache__" not in str(p))
    for g in ["pyproject.toml", "README.md"]:
        p = ROOT / g
        if p.exists():
            paths.append(p)
    return _render(sorted(paths))


def _gather_context_scoped(slug: str, description: str) -> str:
    """Manifest + only the files of the tool(s) named by the item."""
    paths = _manifest_paths()
    tools = _index_tools()
    slug_compact = re.sub(r"[^a-z0-9]", "", slug.lower())
    desc = description.lower()
    matched = [n for n in tools if n.replace("_", "") in slug_compact or n in desc]
    for n in matched:
        for p in tools[n]:
            if p not in paths:
                paths.append(p)
    return _render(paths)


def _call(messages) -> str:
    key = os.environ["DEEPSEEK_API_KEY"]
    model = os.environ.get("DEEPSEEK_MODEL") or DEFAULT_MODEL
    payload = json.dumps(
        {
            "model": model,
            "messages": messages,
            "temperature": 0.2,
            "response_format": {"type": "json_object"},
        }
    ).encode("utf-8")
    req = urllib.request.Request(
        API_URL,
        data=payload,
        headers={
            "Authorization": f"Bearer {key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"DeepSeek HTTP {exc.code}: {exc.read().decode()[:300]}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"DeepSeek connection error: {exc.reason}") from exc
    return data["choices"][0]["message"]["content"]


def implement(slug: str, description: str, feedback: str | None = None) -> Tuple[str, List[Tuple[str, bool]]]:
    """Implement one item via DeepSeek. Returns (commit_subject, written)."""
    context = _gather_context_scoped(slug, description)
    user = (
        f"Backlog item `{slug}`: {description}\n\n"
        f"Repository files (scoped to this item):\n\n{context}"
    )
    if feedback:
        user += (
            "\n\nYour previous attempt did NOT pass the test suite. "
            "Here is the failure output — fix it and return the corrected JSON:\n\n"
            + feedback
        )

    content = _call(
        [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user},
        ]
    )
    try:
        result = json.loads(content)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"DeepSeek returned invalid JSON: {exc}") from exc

    files = result.get("files") or {}
    if not files:
        # Item already satisfied — no change needed; the caller checks it off.
        return result.get("summary") or f"chore: {slug}", []

    written: List[Tuple[str, bool]] = []
    for rel, text in files.items():
        target = ROOT / rel
        try:
            target.resolve().relative_to(ROOT.resolve())
        except ValueError as exc:
            raise RuntimeError(f"refusing to write outside the repo: {rel}") from exc
        existed = target.exists()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")
        written.append((rel, existed))

    return result.get("summary") or f"feat: {slug}", written


def propose_items(count: int = 10):
    """Ask DeepSeek for ``count`` new backlog items, as (slug, title, why)."""
    context = _gather_context_full()
    backlog = (ROOT / "BACKLOG.md").read_text(encoding="utf-8")
    user = (
        "Current repository files:\n\n" + context
        + "\n\nCurrent backlog (do not duplicate):\n\n" + backlog
    )
    content = _call(
        [
            {"role": "system", "content": PROPOSE_PROMPT},
            {"role": "user", "content": user},
        ]
    )
    try:
        data = json.loads(content)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"DeepSeek returned invalid JSON: {exc}") from exc
    items = data.get("items") or []
    return [(it["slug"], it["title"], it.get("why", "")) for it in items]
