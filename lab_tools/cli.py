"""Unified command-line entry point.

Dispatch ``lab-tools <tool> [args...]`` to the matching tool's ``main``.
Each tool owns its argument parsing; this file just routes to it.
"""

from __future__ import annotations

import sys
from typing import Dict, List, Optional

from lab_tools.analytical.elisa import main as elisa_main
from lab_tools.analytical.hplc_sst import main as hplc_sst_main
from lab_tools.analytical.lod_loq import main as lod_loq_main
from lab_tools.lab_ops.hemocytometer import main as hemocytometer_main
from lab_tools.lab_ops.pipette_cal.cli import main as pipette_cal_main
from lab_tools.microbiology.bioburden_spc.chart import main as bioburden_spc_main
from lab_tools.microbiology.cfu import main as cfu_main
from lab_tools.microbiology.growth_curve.main import main as growth_curve_main
from lab_tools.microbiology.media_fill import main as media_fill_main
from lab_tools.microbiology.mpn import main as mpn_main
from lab_tools.molecular.phylogeny import main as phylogeny_main
from lab_tools.molecular.primers import main as primers_main
from lab_tools.molecular.qpcr import main as qpcr_main
from lab_tools.molecular.restriction.cli import main as restriction_main
from lab_tools.sterility.d_z_f0 import main as d_z_f0_main
from lab_tools.sterility.endotoxin.cli import main as endotoxin_main
from lab_tools.sterility.sterility import main as sterility_main

# tool name -> callable(argv) returning an int exit code
TOOLS: Dict[str, object] = {
    "cfu": cfu_main,
    "mpn": mpn_main,
    "media-fill": media_fill_main,
    "growth-curve": growth_curve_main,
    "bioburden-spc": bioburden_spc_main,
    "d-z-f0": d_z_f0_main,
    "sterility": sterility_main,
    "endotoxin": endotoxin_main,
    "qpcr": qpcr_main,
    "primers": primers_main,
    "restriction": restriction_main,
    "phylogeny": phylogeny_main,
    "elisa": elisa_main,
    "lod-loq": lod_loq_main,
    "hplc-sst": hplc_sst_main,
    "hemocytometer": hemocytometer_main,
    "pipette-cal": pipette_cal_main,
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
