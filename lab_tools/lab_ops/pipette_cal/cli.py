#!/usr/bin/env python3
"""Command-line interface for pipette calibration analysis."""

import argparse
import sys
from lab_tools.lab_ops.pipette_cal.io import read_calibration_data
from lab_tools.lab_ops.pipette_cal.core import analyze_calibration


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Analyze gravimetric pipette calibration data against ISO 8655 tolerances."
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Path to input CSV file"
    )
    parser.add_argument(
        "--output", "-o", help="Path to output CSV file (default: stdout)"
    )
    parser.add_argument(
        "--uncertainty",
        action="store_true",
        help="Include combined and expanded (k=2) uncertainty columns.",
    )
    args = parser.parse_args(argv)

    try:
        data = read_calibration_data(args.input)
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        return 1

    if not data:
        print("No valid measurements found.", file=sys.stderr)
        return 1

    results = analyze_calibration(data)

    header = "PipetteID,NominalVolume,TargetVolume,N,MeanVolume_ul,Accuracy_pct,CV_pct,PassFail"
    if args.uncertainty:
        header += ",U_combined_ul,U_expanded_ul"
    out_lines = [header]
    for row in results:
        line = (
            f"{row['PipetteID']},{row['NominalVolume']},{row['TargetVolume']},"
            f"{row['N']},{row['MeanVolume_ul']:.3f},"
            f"{row['Accuracy_pct']:.2f},{row['CV_pct']:.2f},{row['PassFail']}"
        )
        if args.uncertainty:
            uc = row.get("U_combined_ul")
            ue = row.get("U_expanded_ul")
            if isinstance(uc, (int, float)) and isinstance(ue, (int, float)):
                line += f",{uc:.4f},{ue:.4f}"
            else:
                line += ",N/A,N/A"
        out_lines.append(line)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write("\n".join(out_lines) + "\n")
        print(f"Report written to {args.output}")
    else:
        print("\n".join(out_lines))

    return 0


if __name__ == "__main__":
    sys.exit(main())
