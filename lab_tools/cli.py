"""Unified command-line entry point.

Dispatch ``lab-tools <tool> [args...]`` to the matching tool's ``main``.
Each tool owns its argument parsing; this file just routes to it.
"""

from __future__ import annotations

import sys
from typing import Dict, List, Optional

from lab_tools.microbiology import cfu

# name -> callable(argv) returning an int exit code
TOOLS: Dict[str, object] = {
    "cfu": cfu.main,
}


def _usage() -> str:
    return (
        "usage: lab-tools <tool> [args...]\n\n"
        "tools:\n"
        + "".join(f"  {name}\n" for name in sorted(TOOLS))
    )


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if not argv or argv[0] in ("-h", "--help", "help"):
        print(_usage())
        return 0
    name = argv[0]
    if name not in TOOLS:
        print(f"lab-tools: unknown tool '{name}'\n", file=sys.stderr)
        print(_usage(), file=sys.stderr)
        return 2
    return TOOLS[name](argv[1:])  # type: ignore[operator]


if __name__ == "__main__":
    sys.exit(main())
