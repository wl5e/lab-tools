#!/usr/bin/env python3
"""HPLC System Suitability Calculator

GMP-compliant CLI tool to compute theoretical plates, tailing factor,
resolution, and %RSD from chromatographic peak data in CSV format.
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple


def read_csv(filepath: str, delimiter: str = ',') -> List[Dict[str, str]]:
    """Read CSV file and return list of rows as dicts."""
    try:
        with open(filepath, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            if reader.fieldnames is None:
                raise ValueError("CSV file appears to be empty or malformed.")
            rows = list(reader)
            if not rows:
                raise ValueError("CSV file contains no data rows.")
        return rows
    except FileNotFoundError:
        sys.exit(f"Error: File '{filepath}' not found.")
    except csv.Error as e:
        sys.exit(f"CSV parsing error: {e}")


def validate_columns(rows: List[Dict[str, str]], required_cols: List[str], label: str) -> None:
    """Ensure all required columns exist in the first row."""
    if not rows:
        sys.exit("No data to validate.")
    available = set(rows[0].keys())
    missing = [col for col in required_cols if col not in available]
    if missing:
        sys.exit(f"{label}: Missing required column(s): {', '.join(missing)}")


def to_float(value: str, col_name: str, row_idx: int) -> float:
    """Convert string to float with error context."""
    try:
        return float(value)
    except (ValueError, TypeError):
        sys.exit(f"Row {row_idx + 1}, column '{col_name}': expected a number, got '{value}'")


# ---------- Chromatography calculations ----------
def calc_plates_half_height(t_r: float, wh: float) -> float:
    """Theoretical plates using half-height width."""
    return 5.54 * (t_r / wh) ** 2


def calc_plates_base_width(t_r: float, wb: float) -> float:
    """Theoretical plates using base width."""
    return 16.0 * (t_r / wb) ** 2


def calc_tailing_from_asymmetry(asym: float) -> float:
    """Tailing factor approximated from peak asymmetry (As = b/a at 10 % height)."""
    return (1.0 + asym) / 2.0


def calc_resolution_base(t_r1: float, t_r2: float, wb1: float, wb2: float) -> float:
    """USP resolution using base widths."""
    return 2.0 * (t_r2 - t_r1) / (wb1 + wb2)


def calc_resolution_half(t_r1: float, t_r2: float, wh1: float, wh2: float) -> float:
    """Approximation of resolution using half-height widths (factor 1.18)."""
    return 1.18 * (t_r2 - t_r1) / (wh1 + wh2)


# ---------- Core processing ----------
def process_data(
    rows: List[Dict[str, str]],
    method: str,
    tailing_method: str,
    resolution_method: str,
    injection_col: Optional[str]
) -> List[dict]:
    """Compute SST metrics for each peak and return list of result dicts."""
    # Determine required columns based on methods
    common_cols = ['peak_name', 'retention_time']
    validate_columns(rows, common_cols, "Base columns")

    width_col = 'width_half' if method == 'half' else 'width_base'
    validate_columns(rows, [width_col], f"For method '{method}'")

    if tailing_method == '5pct':
        validate_columns(rows, ['a', 'width_5pct'], "Tailing 5% method")
    else:
        validate_columns(rows, ['asymmetry'], "Tailing asymmetry method")

    if resolution_method == 'base':
        res_width_col = 'width_base'
        validate_columns(rows, ['width_base'], "Resolution base width")
    else:
        res_width_col = 'width_half'
        validate_columns(rows, ['width_half'], "Resolution half-height")

    # Parse numeric fields, group by peak_name for replicate analysis
    peaks: Dict[str, List[dict]] = defaultdict(list)
    for idx, row in enumerate(rows):
        try:
            name = row['peak_name'].strip()
            tr = to_float(row['retention_time'], 'retention_time', idx)
            width_val = to_float(row[width_col], width_col, idx)
            asym_value = None
            if tailing_method == 'asymmetry':
                asym_value = to_float(row['asymmetry'], 'asymmetry', idx)
            a_val = None
            w5_val = None
            if tailing_method == '5pct':
                a_val = to_float(row['a'], 'a', idx)
                w5_val = to_float(row['width_5pct'], 'width_5pct', idx)
            res_width = to_float(row[res_width_col], res_width_col, idx)

            peaks[name].append({
                'tr': tr,
                'width': width_val,
                'asym': asym_value,
                'a': a_val,
                'width_5pct': w5_val,
                'res_width': res_width
            })
        except SystemExit:
            raise   # re-raise sys.exit from to_float
        except Exception as e:
            sys.exit(f"Row {idx + 1}: unexpected error: {e}")

    if not peaks:
        sys.exit("No valid peak data found.")

    # Sort peaks by mean retention time for resolution order
    peak_names = list(peaks.keys())
    mean_tr = {name: statistics.mean([p['tr'] for p in peaks[name]]) for name in peak_names}
    sorted_names = sorted(peak_names, key=lambda n: mean_tr[n])

    # Calculate metrics
    results = []
    for name in sorted_names:
        replicates = peaks[name]
        tr_vals = [p['tr'] for p in replicates]
        tr_mean = statistics.mean(tr_vals)
        tr_rsd = (statistics.stdev(tr_vals) / tr_mean * 100.0) if len(tr_vals) > 1 and tr_mean != 0 else float('nan')

        # Take first replicate values for width, asymmetry, etc. (they should be the same for all replicates)
        rep0 = replicates[0]
        # Plate count
        if method == 'half':
            plates = calc_plates_half_height(tr_mean, rep0['width'])
        else:
            plates = calc_plates_base_width(tr_mean, rep0['width'])

        # Tailing factor
        if tailing_method == 'asymmetry':
            tailing = calc_tailing_from_asymmetry(rep0['asym'])
        else:
            # T = W_0.05 / (2 * a)
            tailing = rep0['width_5pct'] / (2.0 * rep0['a']) if rep0['a'] != 0 else float('inf')

        results.append({
            'name': name,
            'tr_mean': tr_mean,
            'plates': plates,
            'tailing': tailing,
            'resolution': None,
            'tr_rsd': tr_rsd,
            'tr_count': len(tr_vals)
        })

    # Compute resolution between adjacent peaks
    for i in range(len(results) - 1):
        p1 = peaks[sorted_names[i]][0]  # first replicate of peak i
        p2 = peaks[sorted_names[i+1]][0]
        if resolution_method == 'base':
            if p1['res_width'] and p2['res_width']:
                res = calc_resolution_base(mean_tr[sorted_names[i]], mean_tr[sorted_names[i+1]],
                                           p1['res_width'], p2['res_width'])
            else:
                res = float('nan')
        else:
            if p1['res_width'] and p2['res_width']:
                res = calc_resolution_half(mean_tr[sorted_names[i]], mean_tr[sorted_names[i+1]],
                                           p1['res_width'], p2['res_width'])
            else:
                res = float('nan')
        results[i]['resolution'] = res

    return results


def format_table(results: List[dict]) -> str:
    """Format results as an aligned ASCII table."""
    if not results:
        return "No results."
    headers = ["Peak", "Rt (mean)", "Plates", "Tailing", "Resolution", "%RSD", "#Inj"]
    col_widths = [len(h) for h in headers]
    rows = []
    for r in results:
        row = [
            r['name'],
            f"{r['tr_mean']:.4f}",
            f"{r['plates']:.1f}",
            f"{r['tailing']:.3f}" if math.isfinite(r['tailing']) else "Err",
            f"{r['resolution']:.2f}" if r['resolution'] is not None and math.isfinite(r['resolution']) else "N/A",
            f"{r['tr_rsd']:.2f}" if not math.isnan(r['tr_rsd']) else "N/A",
            str(r['tr_count'])
        ]
        rows.append(row)
        for idx, val in enumerate(row):
            col_widths[idx] = max(col_widths[idx], len(val))

    # Build separator and format
    sep = "-".join("-" * w for w in col_widths)
    header_line = "  ".join(h.ljust(w) for h, w in zip(headers, col_widths))
    lines = [header_line, sep]
    for row in rows:
        lines.append("  ".join(v.ljust(w) for v, w in zip(row, col_widths)))
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="HPLC System Suitability Calculator")
    parser.add_argument("input", help="CSV file with peak data")
    parser.add_argument("-m", "--method", choices=["half", "base"], default="base",
                        help="Plate count method (half-height or base width)")
    parser.add_argument("-t", "--tailing-method", choices=["asymmetry", "5pct"], default="asymmetry",
                        help="Tailing factor calculation method")
    parser.add_argument("-r", "--resolution-method", choices=["base", "half"], default="base",
                        help="Resolution calculation method")
    parser.add_argument("-i", "--injection-column", type=str, default=None,
                        help="Column name for injection identifier (optional)")
    parser.add_argument("-d", "--delimiter", type=str, default=",",
                        help="CSV delimiter (default: comma)")
    args = parser.parse_args(argv)

    # Read input
    rows = read_csv(args.input, args.delimiter)

    # Process
    results = process_data(rows, args.method, args.tailing_method, args.resolution_method, args.injection_column)

    # Print table
    table = format_table(results)
    print(table)
    print()
    print("Calculation completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
