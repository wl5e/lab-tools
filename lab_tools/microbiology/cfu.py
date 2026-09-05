#!/usr/bin/env python3
"""
Dilution Plate CFU Calculator
==============================
A GMP-compliant tool for calculating colony-forming units (CFU)
from serial dilution plate count data with uncertainty propagation.

Supports spread plate and pour plate methods, handles TNTC/TFTC flags,
and generates formatted reports suitable for pharmaceutical microbiology
documentation.

Author: Collins Amatu Gorgerat
License: MIT
"""

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Countable range per USP <61>/<62> and EP 2.6.12/2.6.13 guidelines
DEFAULT_COUNTABLE_RANGE: Tuple[int, int] = (30, 300)
YEAST_MOLD_RANGE: Tuple[int, int] = (15, 150)


@dataclass
class DilutionPlate:
    """Represents a single plate count at a given dilution level."""
    dilution_factor: float
    volume_plated_ml: float
    colony_count: int
    is_tntc: bool = False
    is_tftc: bool = False
    replicate_id: Optional[str] = None


@dataclass
class CFUResult:
    """Complete result of a CFU calculation with metadata."""
    cfu_per_ml: float
    cfu_per_g: Optional[float]
    uncertainty_percent: float
    plates_used: int
    dilutions_used: List[float] = field(default_factory=list)
    method: str = "spread_plate"
    countable_range: Tuple[int, int] = DEFAULT_COUNTABLE_RANGE
    notes: List[str] = field(default_factory=list)


class DilutionCFUCalculator:
    """
    Calculate CFU concentrations from serial dilution plate count data.

    Implements USP/EP harmonized methodology for microbial enumeration:
    - Countable range filtering (30-300 CFU for bacteria, 15-150 for yeast/mold)
    - Weighted averaging across dilution levels
    - Uncertainty propagation using coefficient of variation
    - TNTC/TFTC handling with fallback logic
    """

    def __init__(
        self,
        countable_range: Tuple[int, int] = DEFAULT_COUNTABLE_RANGE
    ) -> None:
        self.countable_min, self.countable_max = countable_range

    def is_countable(self, count: int) -> bool:
        """Check if a colony count falls within the valid countable range."""
        return self.countable_min <= count <= self.countable_max

    def calculate_from_plates(
        self,
        plates: List[DilutionPlate],
        sample_mass_g: Optional[float] = None,
        dilution_blank_ml: float = 9.0,
        method: str = "spread_plate",
    ) -> CFUResult:
        """
        Calculate CFU/mL and optionally CFU/g from plate count data.

        Args:
            plates: List of DilutionPlate objects from one or more dilutions.
            sample_mass_g: Sample mass in grams. If provided, CFU/g is calculated.
            dilution_blank_ml: Volume (mL) of diluent in the primary dilution blank.
            method: Plating method used ("spread_plate" or "pour_plate").

        Returns:
            CFUResult with calculated concentrations and uncertainty.

        Raises:
            ValueError: If no usable plates are available.
        """
        if not plates:
            raise ValueError("No plates provided for calculation.")

        notes: List[str] = []
        usable_plates: List[DilutionPlate] = []

        for plate in plates:
            if plate.is_tntc:
                notes.append(
                    f"Dilution {plate.dilution_factor:.1e} "
                    f"(rep {plate.replicate_id or '?'}): TNTC — excluded"
                )
                continue
            if plate.is_tftc:
                notes.append(
                    f"Dilution {plate.dilution_factor:.1e} "
                    f"(rep {plate.replicate_id or '?'}): TFTC — excluded"
                )
                continue
            if self.is_countable(plate.colony_count):
                usable_plates.append(plate)
            elif plate.colony_count > self.countable_max:
                notes.append(
                    f"Dilution {plate.dilution_factor:.1e} "
                    f"({plate.colony_count} CFU): exceeds countable max "
                    f"({self.countable_max}) — excluded"
                )
            else:
                notes.append(
                    f"Dilution {plate.dilution_factor:.1e} "
                    f"({plate.colony_count} CFU): below countable min "
                    f"({self.countable_min}) — excluded"
                )

        # Fallback: if no plates in countable range, use all non-TNTC/TFTC plates
        if not usable_plates:
            fallback = [
                p for p in plates
                if not p.is_tntc and not p.is_tftc and p.colony_count > 0
            ]
            if fallback:
                usable_plates = fallback
                notes.append(
                    "WARNING: No plates within countable range "
                    f"({self.countable_min}–{self.countable_max}). "
                    "Using all available plate data as fallback."
                )
            else:
                raise ValueError(
                    "No usable plates found. All plates are TNTC, TFTC, or zero. "
                    "Consider adjusting dilution series and re-plating."
                )

        # Calculate CFU/mL for each usable plate
        cfu_values: List[float] = []
        dilution_groups: Dict[float, List[float]] = {}

        for plate in usable_plates:
            cfu = plate.colony_count / (plate.volume_plated_ml * plate.dilution_factor)
            cfu_values.append(cfu)
            dilution_groups.setdefault(plate.dilution_factor, []).append(cfu)

        # Weighted average across all usable plates
        mean_cfu = sum(cfu_values) / len(cfu_values)

        # Uncertainty: coefficient of variation (CV%)
        if len(cfu_values) > 1:
            variance = sum((x - mean_cfu) ** 2 for x in cfu_values) / (len(cfu_values) - 1)
            std_dev = math.sqrt(variance)
            cv_percent = (std_dev / mean_cfu) * 100 if mean_cfu > 0 else float("inf")
        else:
            # Single plate: Poisson-based estimate (CV ≈ 1/√n)
            avg_count = sum(p.colony_count for p in usable_plates) / len(usable_plates)
            cv_percent = (1.0 / math.sqrt(avg_count)) * 100 if avg_count > 0 else float("inf")

        # CFU/g calculation
        cfu_per_g: Optional[float] = None
        if sample_mass_g is not None and sample_mass_g > 0:
            total_primary_volume = dilution_blank_ml + sample_mass_g
            cfu_per_g = mean_cfu * total_primary_volume / sample_mass_g

        return CFUResult(
            cfu_per_ml=round(mean_cfu, 2),
            cfu_per_g=round(cfu_per_g, 2) if cfu_per_g is not None else None,
            uncertainty_percent=round(cv_percent, 2),
            plates_used=len(usable_plates),
            dilutions_used=sorted(dilution_groups.keys()),
            method=method,
            countable_range=(self.countable_min, self.countable_max),
            notes=notes,
        )

    @staticmethod
    def plates_from_serial_dilution(
        counts: List[int],
        initial_dilution: float = 0.1,
        dilution_step: float = 10.0,
        num_dilutions: int = 6,
        countable_min: int = DEFAULT_COUNTABLE_RANGE[0],
        volume_plated_ml: float = 0.1,
        replicates: int = 2,
    ) -> List[DilutionPlate]:
        """
        Convenience method: generate plate objects from a standard serial dilution.

        Args:
            counts: Colony counts per plate, ordered from least to most dilute.
            initial_dilution: First dilution factor (e.g., 0.1 for 10⁻¹).
            dilution_step: Step factor between dilutions (default 10 for 1:10).
            num_dilutions: Total number of dilution levels.
            volume_plated_ml: Volume plated per plate in mL.
            replicates: Number of replicate plates per dilution level.

        Returns:
            List of DilutionPlate objects.
        """
        plates: List[DilutionPlate] = []
        for i in range(num_dilutions):
            dilution = initial_dilution * (dilution_step ** i)
            for r in range(replicates):
                idx = i * replicates + r
                if idx >= len(counts):
                    break
                raw = counts[idx]
                is_tntc = raw < 0
                is_tftc = raw >= 0 and raw < countable_min
                plates.append(DilutionPlate(
                    dilution_factor=dilution,
                    volume_plated_ml=volume_plated_ml,
                    colony_count=abs(raw),
                    is_tntc=is_tntc,
                    is_tftc=is_tftc,
                    replicate_id=chr(65 + r),
                ))
        return plates

    def format_result(self, result: CFUResult) -> str:
        """Format a CFUResult as a human-readable GMP-style report."""
        lines: List[str] = []
        sep = "=" * 62
        lines.append(sep)
        lines.append("        DILUTION PLATE CFU ANALYSIS REPORT")
        lines.append(sep)
        lines.append(f"  Method:                  {result.method.replace('_', ' ').title()}")
        lines.append(f"  Countable Range:         {result.countable_range[0]}–{result.countable_range[1]} CFU/plate")
        lines.append(f"  Plates in Analysis:      {result.plates_used}")
        lines.append(f"  Dilutions Used:          {', '.join(f'{d:.1e}' for d in result.dilutions_used)}")
        lines.append("-" * 62)
        lines.append(f"  CFU/mL:                  {result.cfu_per_ml:.2e}")
        if result.cfu_per_g is not None:
            lines.append(f"  CFU/g:                   {result.cfu_per_g:.2e}")
        lines.append(f"  Uncertainty (CV%):       {result.uncertainty_percent:.2f}%")
        lines.append("-" * 62)
        if result.notes:
            lines.append("  Notes:")
            for note in result.notes:
                lines.append(f"    • {note}")
        lines.append(sep)
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Input Parsing
# ---------------------------------------------------------------------------

def _parse_count_value(raw: str) -> Tuple[int, bool, bool]:
    """
    Parse a colony count string into (abs_count, is_tntc, is_tftc).

    Accepts numeric values and the special flags: TNTC, >MAX, TFTC, <MIN.
    """
    val = raw.strip().upper()
    if val in ("TNTC", ">300", ">MAX"):
        return (-1, True, False)
    if val in ("TFTC", "<30", "<MIN"):
        return (0, False, True)
    if val == "0":
        return (0, False, True)
    try:
        count = int(val)
        if count < 0:
            raise ValueError(f"Negative colony count not allowed: {count}")
        return (count, False, False)
    except ValueError:
        raise ValueError(
            f"Invalid colony count value: '{raw}'. "
            "Expected integer, 'TNTC', or 'TFTC'."
        ) from None


def parse_csv_input(filepath: str) -> List[DilutionPlate]:
    """Parse plate count data from a CSV file."""
    plates: List[DilutionPlate] = []
    with open(filepath, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV file appears to be empty or has no header row.")
        for row_num, row in enumerate(reader, start=2):
            try:
                dil_key = next(
                    (k for k in row if k and "dilution" in k.lower() and "factor" in k.lower()),
                    None,
                )
                vol_key = next(
                    (k for k in row if k and "volume" in k.lower() and "plate" in k.lower()),
                    None,
                )
                count_key = next(
                    (k for k in row if k and ("count" in k.lower() or "colony" in k.lower())),
                    None,
                )
                rep_key = next(
                    (k for k in row if k and "replicate" in k.lower()),
                    None,
                )

                if not dil_key or not count_key:
                    raise KeyError("Missing required columns (dilution_factor, colony_count).")

                dilution = float(row[dil_key])
                volume = float(row[vol_key]) if vol_key else 0.1
                count, is_tntc, is_tftc = _parse_count_value(row[count_key])
                replicate = row[rep_key].strip() if rep_key and row.get(rep_key) else None

                plates.append(DilutionPlate(
                    dilution_factor=dilution,
                    volume_plated_ml=volume,
                    colony_count=count,
                    is_tntc=is_tntc,
                    is_tftc=is_tftc,
                    replicate_id=replicate,
                ))
            except (ValueError, KeyError) as exc:
                print(f"Warning: Skipping row {row_num}: {exc}", file=sys.stderr)
    return plates


def parse_json_input(filepath: str) -> List[DilutionPlate]:
    """Parse plate count data from a JSON file."""
    with open(filepath, "r") as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("JSON input must be a list of plate objects.")

    plates: List[DilutionPlate] = []
    for idx, item in enumerate(data):
        try:
            dilution = float(item.get("dilution_factor", item.get("dilution", 1.0)))
            volume = float(item.get("volume_plated_ml", item.get("volume", 0.1)))
            raw_count = str(item.get("colony_count", item.get("count", "")))
            count, is_tntc, is_tftc = _parse_count_value(raw_count)
            replicate = item.get("replicate_id", item.get("replicate", None))
            if replicate is not None:
                replicate = str(replicate)

            plates.append(DilutionPlate(
                dilution_factor=dilution,
                volume_plated_ml=volume,
                colony_count=count,
                is_tntc=is_tntc,
                is_tftc=is_tftc,
                replicate_id=replicate,
            ))
        except (ValueError, KeyError, TypeError) as exc:
            print(f"Warning: Skipping item {idx}: {exc}", file=sys.stderr)
    return plates


def generate_example_csv(filepath: str) -> None:
    """Write an example CSV template for users to fill in."""
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["dilution_factor", "volume_plated_ml", "colony_count", "replicate_id"])
        writer.writerow(["0.1", "0.1", "TNTC", "A"])
        writer.writerow(["0.1", "0.1", "TNTC", "B"])
        writer.writerow(["0.01", "0.1", "245", "A"])
        writer.writerow(["0.01", "0.1", "231", "B"])
        writer.writerow(["0.001", "0.1", "42", "A"])
        writer.writerow(["0.001", "0.1", "38", "B"])
        writer.writerow(["0.0001", "0.1", "3", "A"])
        writer.writerow(["0.0001", "0.1", "5", "B"])
    print(f"Example CSV template written to: {filepath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dilution-cfu",
        description="Calculate CFU concentrations from serial dilution plate counts.",
        epilog="GMP-compliant tool for microbiological enumeration assays.",
    )
    sub = parser.add_subparsers(dest="command", help="Available commands")

    # --- calculate ---
    calc = sub.add_parser("calculate", help="Calculate CFU from a CSV or JSON input file")
    calc.add_argument("input", help="Path to input file (.csv or .json)")
    calc.add_argument("--sample-mass-g", type=float, default=None,
                      help="Sample mass in grams for CFU/g calculation")
    calc.add_argument("--diluent-ml", type=float, default=9.0,
                      help="Volume of diluent in primary blank (default: 9.0 mL)")
    calc.add_argument("--method", choices=["spread_plate", "pour_plate"],
                      default="spread_plate", help="Plating method used")
    calc.add_argument("--countable-min", type=int, default=None,
                      help="Override minimum countable colonies (default: 30 for bacteria)")
    calc.add_argument("--countable-max", type=int, default=None,
                      help="Override maximum countable colonies (default: 300 for bacteria)")
    calc.add_argument("--organism", choices=["bacteria", "yeast_mold"],
                      default="bacteria", help="Organism type for default countable range")
    calc.add_argument("--json-output", action="store_true",
                      help="Emit machine-readable JSON instead of formatted text")

    # --- quick ---
    quick = sub.add_parser("quick", help="Quick calculation from command-line arguments")
    quick.add_argument("--counts", type=str, required=True,
                       help="Comma-separated colony counts (use TNTC/TFTC for special flags)")
    quick.add_argument("--dilutions", type=str, required=True,
                       help="Comma-separated dilution factors (e.g., 0.1,0.01,0.001)")
    quick.add_argument("--volume", type=float, default=0.1,
                       help="Volume plated per plate in mL (default: 0.1)")
    quick.add_argument("--replicates", type=int, default=2,
                       help="Number of replicate plates per dilution (default: 2)")
    quick.add_argument("--sample-mass-g", type=float, default=None,
                       help="Sample mass in grams for CFU/g calculation")

    # --- example ---
    example = sub.add_parser("example", help="Generate an example CSV template file")
    example.add_argument("output", nargs="?", default="plate_data_example.csv",
                         help="Output file path (default: plate_data_example.csv)")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "example":
        generate_example_csv(args.output)
        return 0

    if args.command == "calculate":
        # Determine countable range
        if args.organism == "yeast_mold":
            default_range = YEAST_MOLD_RANGE
        else:
            default_range = DEFAULT_COUNTABLE_RANGE

        cmin = args.countable_min if args.countable_min is not None else default_range[0]
        cmax = args.countable_max if args.countable_max is not None else default_range[1]

        if cmin >= cmax:
            print("Error: countable-min must be less than countable-max.", file=sys.stderr)
            return 1

        # Parse input file
        fpath: str = args.input
        try:
            if fpath.endswith(".json"):
                plates = parse_json_input(fpath)
            elif fpath.endswith(".csv"):
                plates = parse_csv_input(fpath)
            else:
                print("Error: Input file must have .csv or .json extension.", file=sys.stderr)
                return 1
        except FileNotFoundError:
            print(f"Error: File not found: {fpath}", file=sys.stderr)
            return 1
        except (ValueError, json.JSONDecodeError) as exc:
            print(f"Error parsing input: {exc}", file=sys.stderr)
            return 1

        if not plates:
            print("Error: No valid plate records found in input file.", file=sys.stderr)
            return 1

        calculator = DilutionCFUCalculator(countable_range=(cmin, cmax))

        try:
            result = calculator.calculate_from_plates(
                plates,
                sample_mass_g=args.sample_mass_g,
                dilution_blank_ml=args.diluent_ml,
                method=args.method,
            )
        except ValueError as exc:
            print(f"Calculation error: {exc}", file=sys.stderr)
            return 1

        if args.json_output:
            print(json.dumps({
                "cfu_per_ml": result.cfu_per_ml,
                "cfu_per_g": result.cfu_per_g,
                "uncertainty_percent": result.uncertainty_percent,
                "plates_used": result.plates_used,
                "dilutions_used": result.dilutions_used,
                "method": result.method,
                "countable_range": list(result.countable_range),
                "notes": result.notes,
            }, indent=2))
        else:
            print(calculator.format_result(result))

        return 0

    if args.command == "quick":
        count_strs = [c.strip() for c in args.counts.split(",")]
        dil_strs = [d.strip() for d in args.dilutions.split(",")]

        expected = len(dil_strs) * args.replicates
        if len(count_strs) != expected:
            print(
                f"Error: Expected {expected} colony counts "
                f"({len(dil_strs)} dilutions × {args.replicates} replicates), "
                f"but got {len(count_strs)}.",
                file=sys.stderr,
            )
            return 1

        plates: List[DilutionPlate] = []
        for i, dil_str in enumerate(dil_strs):
            try:
                dilution = float(dil_str)
            except ValueError:
                print(f"Error: Invalid dilution factor: '{dil_str}'", file=sys.stderr)
                return 1
            if dilution <= 0 or dilution > 1:
                print(
                    f"Warning: Dilution factor {dilution} is outside typical range (0 < x ≤ 1).",
                    file=sys.stderr,
                )
            for r in range(args.replicates):
                idx = i * args.replicates + r
                try:
                    count, is_tntc, is_tftc = _parse_count_value(count_strs[idx])
                except ValueError as exc:
                    print(f"Error: {exc}", file=sys.stderr)
                    return 1
                plates.append(DilutionPlate(
                    dilution_factor=dilution,
                    volume_plated_ml=args.volume,
                    colony_count=count,
                    is_tntc=is_tntc,
                    is_tftc=is_tftc,
                    replicate_id=chr(65 + r),
                ))

        calculator = DilutionCFUCalculator()
        try:
            result = calculator.calculate_from_plates(plates, sample_mass_g=args.sample_mass_g)
        except ValueError as exc:
            print(f"Calculation error: {exc}", file=sys.stderr)
            return 1

        print(calculator.format_result(result))
        return 0

    # No command or unrecognized
    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
