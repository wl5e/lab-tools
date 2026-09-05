#!/usr/bin/env python3
"""
Hemocytometer Cell Counter

Command-line tool for GMP laboratories to calculate cell concentration,
viability, total viable cells, and dilution volumes from haemocytometer counts.
"""

import argparse
import sys
from typing import List, Optional, Tuple


class CellCounter:
    """Handles all hemocytometer calculations."""

    # Standard Neubauer chamber: 0.1 µL volume per square -> 10^4 factor for cells/mL
    DEFAULT_CHAMBER_FACTOR = 10**4  # per mL

    def __init__(
        self,
        live_cells: int,
        dead_cells: int,
        squares_counted: int,
        dilution_factor: float = 1.0,
        chamber_factor: float = DEFAULT_CHAMBER_FACTOR,
    ):
        if live_cells < 0 or dead_cells < 0:
            raise ValueError("Cell counts must be non-negative.")
        if squares_counted <= 0:
            raise ValueError("Number of squares counted must be positive.")
        if dilution_factor <= 0:
            raise ValueError("Dilution factor must be positive.")

        self.live = live_cells
        self.dead = dead_cells
        self.squares = squares_counted
        self.dilution = dilution_factor
        self.chamber_factor = chamber_factor

    @property
    def total_cells_counted(self) -> int:
        return self.live + self.dead

    @property
    def viability(self) -> float:
        total = self.total_cells_counted
        if total == 0:
            return 0.0
        return (self.live / total) * 100.0

    @property
    def average_per_square(self) -> float:
        return self.total_cells_counted / self.squares

    @property
    def cells_per_ml(self) -> float:
        """Cells per mL in the original (undiluted) sample."""
        return self.average_per_square * self.dilution * self.chamber_factor

    @property
    def viable_cells_per_ml(self) -> float:
        """Viable cells per mL in the original sample."""
        if self.total_cells_counted == 0:
            return 0.0
        return self.cells_per_ml * (self.live / self.total_cells_counted)

    def total_viable_cells(self, volume_ml: float) -> float:
        """Total viable cells in given volume (mL) of original sample."""
        if volume_ml < 0:
            raise ValueError("Volume must be non-negative.")
        return self.viable_cells_per_ml * volume_ml

    def volume_for_cells(self, desired_cells: float) -> float:
        """
        Volume (mL) of original sample required to obtain desired number
        of viable cells. Returns infinity if no viable cells exist.
        """
        if desired_cells <= 0:
            raise ValueError("Desired cells must be positive.")
        vpm = self.viable_cells_per_ml
        if vpm == 0:
            return float('inf')
        return desired_cells / vpm

    def dilution_volumes(
        self, target_cells_per_ml: float, final_volume_ml: float
    ) -> Tuple[float, float]:
        """
        Calculate volumes to achieve a target concentration using C1V1 = C2V2.
        Returns (volume_of_sample, volume_of_diluent) in mL.
        """
        if target_cells_per_ml <= 0 or final_volume_ml <= 0:
            raise ValueError("Target concentration and final volume must be positive.")
        c1 = self.viable_cells_per_ml
        if c1 == 0:
            raise ValueError("Original sample has zero viable cells; cannot dilute.")
        v1 = (target_cells_per_ml * final_volume_ml) / c1
        if v1 > final_volume_ml:
            raise ValueError(
                "Required sample volume exceeds final volume. "
                "Increase final volume or reduce target concentration."
            )
        v2 = final_volume_ml - v1
        return v1, v2


def interactive_count() -> Tuple[int, int, int, float]:
    """Prompt user for count parameters."""
    print("Enter hemocytometer count data:")
    while True:
        try:
            live = int(input("Live cells: ").strip())
            break
        except ValueError:
            print("Invalid integer. Try again.")
    while True:
        try:
            dead = int(input("Dead cells: ").strip())
            break
        except ValueError:
            print("Invalid integer. Try again.")
    while True:
        try:
            squares = int(input("Squares counted (default 4): ").strip() or "4")
            break
        except ValueError:
            print("Invalid integer. Try again.")
    while True:
        try:
            dil = input("Dilution factor (default 1.0): ").strip()
            dilution = float(dil) if dil else 1.0
            break
        except ValueError:
            print("Invalid float. Try again.")
    return live, dead, squares, dilution


def print_table(rows: list, headers: list):
    """Format a simple table with borders."""
    col_widths = [max(len(str(h)), max(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    border = "┌" + "─" * (sum(col_widths) + 3 * (len(headers) - 1) + 2) + "┐"
    header_line = "│ " + " │ ".join(h.ljust(w) for h, w in zip(headers, col_widths)) + " │"
    separator = "├" + "─" * (sum(col_widths) + 3 * (len(headers) - 1) + 2) + "┤"
    bottom = "└" + "─" * (sum(col_widths) + 3 * (len(headers) - 1) + 2) + "┘"
    print(border)
    print(header_line)
    print(separator)
    for row in rows:
        print("│ " + " │ ".join(str(v).ljust(w) for v, w in zip(row, col_widths)) + " │")
    print(bottom)


def print_count_results(cc: CellCounter, volume: Optional[float] = None):
    """Display counting results in a formatted table."""
    rows = [
        ("Total cells counted", str(cc.total_cells_counted)),
        ("Average per square", f"{cc.average_per_square:.2f}"),
        ("Viability", f"{cc.viability:.2f}%"),
        ("Cells/mL (original)", f"{cc.cells_per_ml:,.2f}"),
        ("Viable cells/mL (orig)", f"{cc.viable_cells_per_ml:,.2f}"),
    ]
    headers = ["Metric", "Value"]
    print_table(rows, headers)

    if volume is not None:
        try:
            total = cc.total_viable_cells(volume)
            rows = [(f"Total viable cells in {volume:.3f} mL", f"{total:,.2f}")]
            print_table(rows, headers)
        except ValueError as e:
            print(f"Note: {e}", file=sys.stderr)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Hemocytometer Cell Counter – GMP-compliant cell counting and dilution calculator."
    )
    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # count command
    count_parser = subparsers.add_parser('count', help='Calculate from hemocytometer counts')
    count_group = count_parser.add_mutually_exclusive_group(required=True)
    count_group.add_argument('--live', type=int, help='Number of live cells counted')
    count_group.add_argument('-i', '--interactive', action='store_true', help='Interactive input mode')
    count_parser.add_argument('--dead', type=int, default=0, help='Number of dead cells counted')
    count_parser.add_argument('--squares', type=int, default=4, help='Number of squares counted (default 4)')
    count_parser.add_argument('--dilution', type=float, default=1.0, help='Dilution factor (default 1.0)')
    count_parser.add_argument('--volume', type=float, help='Original sample volume (mL) for total viable cells')

    # dilute command
    dil_parser = subparsers.add_parser('dilute', help='Calculate volumes to achieve target concentration')
    dil_parser.add_argument('--live', type=int, required=True, help='Live cells counted')
    dil_parser.add_argument('--dead', type=int, default=0, help='Dead cells counted')
    dil_parser.add_argument('--squares', type=int, default=4, help='Number of squares counted')
    dil_parser.add_argument('--dilution', type=float, default=1.0, help='Dilution factor of counted sample')
    dil_parser.add_argument('--target-conc', type=float, required=True, help='Target viable cells/mL')
    dil_parser.add_argument('--final-volume', type=float, required=True, help='Final desired volume (mL)')

    # volume command
    vol_parser = subparsers.add_parser('volume', help='Volume needed for desired number of viable cells')
    vol_parser.add_argument('--live', type=int, required=True, help='Live cells counted')
    vol_parser.add_argument('--dead', type=int, default=0, help='Dead cells counted')
    vol_parser.add_argument('--squares', type=int, default=4, help='Number of squares counted')
    vol_parser.add_argument('--dilution', type=float, default=1.0, help='Dilution factor')
    vol_parser.add_argument('--desired-cells', type=float, required=True, help='Desired number of viable cells')

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    try:
        if args.command == 'count':
            if args.interactive:
                live, dead, squares, dilution = interactive_count()
            else:
                if args.live is None:
                    print("Error: --live is required (or use --interactive).", file=sys.stderr)
                    return 1
                live = args.live
                dead = args.dead
                squares = args.squares
                dilution = args.dilution

            cc = CellCounter(live, dead, squares, dilution)
            print_count_results(cc, args.volume)

        elif args.command == 'dilute':
            cc = CellCounter(args.live, args.dead, args.squares, args.dilution)
            v_sample, v_diluent = cc.dilution_volumes(
                args.target_conc, args.final_volume
            )
            rows = [
                ("Target concentration", f"{args.target_conc:,.0f} cells/mL"),
                ("Final volume", f"{args.final_volume:.3f} mL"),
                ("Sample to add", f"{v_sample:.3f} mL"),
                ("Diluent to add", f"{v_diluent:.3f} mL"),
            ]
            print_table(rows, ["Parameter", "Value"])

        elif args.command == 'volume':
            cc = CellCounter(args.live, args.dead, args.squares, args.dilution)
            vol = cc.volume_for_cells(args.desired_cells)
            if vol == float('inf'):
                print("Cannot obtain desired cells: no viable cells in sample.")
            else:
                rows = [
                    ("Desired viable cells", f"{args.desired_cells:,.0f}"),
                    ("Sample volume required", f"{vol:.3f} mL"),
                ]
                print_table(rows, ["Metric", "Value"])
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
