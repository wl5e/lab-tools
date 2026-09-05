#!/usr/bin/env python3
"""
LOD/LOQ Calculator
Calculate Limit of Detection (LOD) and Limit of Quantification (LOQ)
from calibration data using linear regression, compliant with ICH Q2(R1).
"""

import argparse
import csv
import math
import sys
from typing import List, Tuple, Dict


def parse_csv(file_path: str, delimiter: str = ',', has_header: bool = True) -> Tuple[List[float], List[float]]:
    """
    Read concentration and signal columns from a CSV file.

    Args:
        file_path: path to the CSV file.
        delimiter: field delimiter.
        has_header: if True, skip the first line as header.

    Returns:
        Tuple of (concentrations, signals) as lists of floats.

    Raises:
        FileNotFoundError: if file does not exist.
        ValueError: if file is empty, has too few columns, or contains non-numeric data.
    """
    try:
        with open(file_path, newline='') as f:
            reader = csv.reader(f, delimiter=delimiter)
            if has_header:
                try:
                    next(reader)  # skip header
                except StopIteration:
                    raise ValueError("CSV file is empty (no header found).")
            conc, signal = [], []
            for row_num, row in enumerate(reader, start=1):
                if not row:
                    continue  # skip blank lines
                if len(row) < 2:
                    raise ValueError(f"Row {row_num}: expected 2 columns, got {len(row)}")
                try:
                    c = float(row[0].strip())
                    s = float(row[1].strip())
                except ValueError:
                    raise ValueError(f"Row {row_num}: non-numeric value encountered ('{row[0]}', '{row[1]}')")
                conc.append(c)
                signal.append(s)
            if len(conc) < 3:
                raise ValueError(f"Insufficient calibration points ({len(conc)}). At least 3 are required.")
            return conc, signal
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")


def linear_regression(x: List[float], y: List[float]) -> Dict[str, float]:
    """
    Perform ordinary least squares regression: y = a + b*x.

    Args:
        x: independent variable (concentration).
        y: dependent variable (signal).

    Returns:
        Dictionary with keys: 'slope', 'intercept', 'r_squared', 'sy', 'n'.
    """
    n = len(x)
    if n < 3:
        raise ValueError(f"Need >=3 points for regression, got {n}.")
    if n != len(y):
        raise ValueError("x and y must have the same length.")

    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi * xi for xi in x)
    sum_y2 = sum(yi * yi for yi in y)

    denominator = n * sum_x2 - sum_x * sum_x
    if abs(denominator) < 1e-12:
        raise ValueError("Cannot compute slope (near-zero variance in x).")

    slope = (n * sum_xy - sum_x * sum_y) / denominator
    intercept = (sum_y - slope * sum_x) / n

    y_mean = sum_y / n
    ss_total = sum((yi - y_mean) ** 2 for yi in y)
    if ss_total == 0:
        r_squared = 1.0  # perfect fit if all y are identical, though slope will be 0
    else:
        y_pred = [intercept + slope * xi for xi in x]
        ss_residual = sum((yi - yp) ** 2 for yi, yp in zip(y, y_pred))
        r_squared = 1 - (ss_residual / ss_total)
        # guard against floating-point noise
        if r_squared < 0:
            r_squared = 0.0
        elif r_squared > 1.0:
            r_squared = 1.0

    # Residual standard error (Sy) with n-2 degrees of freedom
    residuals = [yi - (intercept + slope * xi) for xi, yi in zip(x, y)]
    ss_res = sum(r * r for r in residuals)
    if n <= 2:
        raise ValueError("Cannot compute Sy with n <= 2 (no degrees of freedom for error).")
    sy = math.sqrt(ss_res / (n - 2))

    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'sy': sy,
        'n': n
    }


def calculate_lod_loq(slope: float, sy: float) -> Tuple[float, float]:
    """
    Compute LOD and LOQ using ICH Q2(R1) recommendations:
    LOD = 3.3 * sy / slope
    LOQ = 10  * sy / slope

    Args:
        slope: slope of the calibration line.
        sy: residual standard error.

    Returns:
        (lod, loq) tuple.

    Raises:
        ValueError if slope is zero (or negative).
    """
    if slope <= 0:
        raise ValueError(f"Slope must be positive to calculate meaningful LOD/LOQ. Got {slope:.4f}")
    lod = 3.3 * sy / slope
    loq = 10.0 * sy / slope
    return lod, loq


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='Calculate LOD and LOQ from calibration CSV data (ICH Q2(R1)).'
    )
    parser.add_argument('file', help='Path to CSV file')
    parser.add_argument('--delimiter', '-d', default=',',
                        help="Field delimiter (default ',')")
    parser.add_argument('--no-header', action='store_true',
                        help='CSV file has no header row')
    args = parser.parse_args(argv)

    try:
        conc, signal = parse_csv(args.file, delimiter=args.delimiter,
                                 has_header=not args.no_header)
        reg = linear_regression(conc, signal)
        lod, loq = calculate_lod_loq(reg['slope'], reg['sy'])

        # Pretty print results
        print("\nSignal vs Concentration")
        print("─────────────────────────")
        print(f"  Slope:         {reg['slope']:.4f}")
        print(f"  Intercept:     {reg['intercept']:.4f}")
        print(f"  R²:            {reg['r_squared']:.4f}")
        print(f"  Residual Sy:   {reg['sy']:.4f}")
        print(f"  n:             {reg['n']}")
        print()
        print(f"Limit of Detection  (LOD) = {lod:.3f}")
        print(f"Limit of Quantitation (LOQ) = {loq:.3f}")

    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
