#!/usr/bin/env python3
"""Bioburden I-MR control chart generator with Western Electric rules."""
import csv
import sys
from statistics import mean
import math

def load_data(filepath: str, date_col: str = 'date', value_col: str = 'bioburden'):
    """Load dates and bioburden values from CSV file.

    Args:
        filepath: path to CSV file
        date_col: column name for date (default 'date')
        value_col: column name for bioburden value (default 'bioburden')

    Returns:
        Tuple of (list of dates, list of float values)

    Raises:
        ValueError: if columns missing or data cannot be parsed.
    """
    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("CSV file is empty")
        if date_col not in reader.fieldnames or value_col not in reader.fieldnames:
            raise ValueError(f"CSV must contain '{date_col}' and '{value_col}' columns. Found: {reader.fieldnames}")
        dates = []
        values = []
        for row_num, row in enumerate(reader, start=2):
            try:
                val = float(row[value_col])
            except (ValueError, TypeError) as e:
                raise ValueError(f"Invalid bioburden value in row {row_num}: {row[value_col]}") from e
            dates.append(row[date_col])
            values.append(val)
        if not values:
            raise ValueError("No data rows found.")
    return dates, values

def compute_individuals_limits(values: list) -> dict:
    """Compute Shewhart I-chart control limits using moving ranges.

    Args:
        values: list of bioburden measurements in time order.

    Returns:
        dict with keys: xbar, avg_mr, sigma, ucl_x, lcl_x, ucl_mr, lcl_mr,
        moving_ranges, and zone limits ucl1, lcl1, ucl2, lcl2.
    """
    n = len(values)
    if n < 2:
        return {
            'xbar': mean(values) if values else 0,
            'avg_mr': 0.0,
            'sigma': 0.0,
            'ucl_x': 0.0,
            'lcl_x': 0.0,
            'ucl_mr': 0.0,
            'lcl_mr': 0.0,
            'moving_ranges': [],
            'ucl1': 0.0, 'lcl1': 0.0,
            'ucl2': 0.0, 'lcl2': 0.0,
        }
    xbar = mean(values)
    moving_ranges = [abs(values[i] - values[i-1]) for i in range(1, n)]
    avg_mr = mean(moving_ranges) if moving_ranges else 0.0
    # d2 for n=2 is 1.128, used to estimate sigma
    d2 = 1.128
    sigma = avg_mr / d2 if avg_mr > 0 else 0.0

    # Constants for individual chart (using E2 for n=2, D4 for MR chart)
    E2 = 2.66
    D4 = 3.267
    D3 = 0.0
    ucl_x = xbar + E2 * avg_mr
    lcl_x = max(0.0, xbar - E2 * avg_mr)
    ucl_mr = D4 * avg_mr
    lcl_mr = D3

    # Zone boundaries (±1 sigma, ±2 sigma) for Western Electric rules
    ucl1 = xbar + sigma
    lcl1 = max(0.0, xbar - sigma)
    ucl2 = xbar + 2 * sigma
    lcl2 = max(0.0, xbar - 2 * sigma)

    return {
        'xbar': xbar,
        'avg_mr': avg_mr,
        'sigma': sigma,
        'ucl_x': ucl_x,
        'lcl_x': lcl_x,
        'ucl_mr': ucl_mr,
        'lcl_mr': lcl_mr,
        'moving_ranges': moving_ranges,
        'ucl1': ucl1,
        'lcl1': lcl1,
        'ucl2': ucl2,
        'lcl2': lcl2,
    }

def detect_violations(dates: list, values: list, limits: dict) -> list:
    """Detect control chart violations using Western Electric rules.

    Rules:
        1. One point beyond UCL or LCL (A zone).
        2. 7 consecutive points on same side of centerline.
        3. Two of three consecutive points beyond ±2 sigma (same side).
        4. Four of five consecutive points beyond ±1 sigma (same side).

    Args:
        dates: list of date strings (for labeling).
        values: bioburden values.
        limits: dict from compute_individuals_limits().

    Returns:
        List of tuples (point_index, violation_message).
    """
    violations = []
    n = len(values)
    if n == 0:
        return violations

    xbar = limits['xbar']
    ucl = limits['ucl_x']
    lcl = limits['lcl_x']
    ucl1 = limits['ucl1']
    lcl1 = limits['lcl1']
    ucl2 = limits['ucl2']
    lcl2 = limits['lcl2']

    # Rule 1: point outside control limits
    for i, v in enumerate(values):
        if v > ucl or v < lcl:
            violations.append((i, f"Rule 1: Point {i+1} ({dates[i] if i < len(dates) else '?'}) = {v} beyond control limits (UCL={ucl:.2f}, LCL={lcl:.2f})"))

    # Rule 2: 7 points consecutive on same side of mean
    if n >= 7:
        side = [1 if v > xbar else (-1 if v < xbar else 0) for v in values]
        for i in range(n - 6):
            if all(s == 1 for s in side[i:i+7]):
                violations.append((i, f"Rule 2: Seven points in a row above centerline starting at point {i+1} ({dates[i]})"))
            if all(s == -1 for s in side[i:i+7]):
                violations.append((i, f"Rule 2: Seven points in a row below centerline starting at point {i+1} ({dates[i]})"))

    # Rule 3: two out of three consecutive points beyond ±2 sigma (same side)
    if n >= 3:
        for i in range(n - 2):
            seg = values[i:i+3]
            above = [v > ucl2 for v in seg]
            below = [v < lcl2 for v in seg]
            if sum(above) >= 2 and all(v > xbar for v in seg if v > ucl2):
                violations.append((i, f"Rule 3: Two of three consecutive points beyond {ucl2:.2f} (2σ) above mean, at points {i+1}-{i+3}"))
            if sum(below) >= 2 and all(v < xbar for v in seg if v < lcl2):
                violations.append((i, f"Rule 3: Two of three consecutive points below {lcl2:.2f} (2σ) below mean, at points {i+1}-{i+3}"))
    # Rule 4: four out of five consecutive points beyond ±1 sigma (same side)
    if n >= 5:
        for i in range(n - 4):
            seg = values[i:i+5]
            above = [v > ucl1 for v in seg]
            below = [v < lcl1 for v in seg]
            if sum(above) >= 4 and all(v > xbar for v in seg if v > ucl1):
                violations.append((i, f"Rule 4: Four of five consecutive points beyond {ucl1:.2f} (1σ) above mean, at points {i+1}-{i+5}"))
            if sum(below) >= 4 and all(v < xbar for v in seg if v < lcl1):
                violations.append((i, f"Rule 4: Four of five consecutive points below {lcl1:.2f} (1σ) below mean, at points {i+1}-{i+5}"))

    return violations

def _val_to_row(val, max_val, min_val, height):
    """Map value to terminal row index (0=top)."""
    if max_val == min_val:
        return int(height // 2)
    return int((max_val - val) / (max_val - min_val) * (height - 1))

def _draw_horizontal_line(chart, row, width, char, max_val, min_val, height):
    for col in range(width):
        chart[row][col] = char if chart[row][col] == ' ' else '+'

def format_individual_chart(dates, values, limits, violations):
    """Create ASCII I-chart with mean, UCL/LCL, and points."""
    if not values:
        return "No data"
    max_val = max(values + [limits['ucl_x'], limits['lcl_x'], limits['xbar']]) * 1.05
    min_val = min(values + [limits['ucl_x'], limits['lcl_x'], limits['xbar']]) * 0.95
    if min_val < 0:
        min_val = 0
    if max_val == min_val:
        max_val += 1
    height = 20
    width = len(values)
    chart = [[' ' for _ in range(width)] for _ in range(height)]
    # Draw mean, UCL, LCL lines
    cl_row = _val_to_row(limits['xbar'], max_val, min_val, height)
    ucl_row = _val_to_row(limits['ucl_x'], max_val, min_val, height)
    lcl_row = _val_to_row(limits['lcl_x'], max_val, min_val, height)
    _draw_horizontal_line(chart, cl_row, width, '-', max_val, min_val, height)
    _draw_horizontal_line(chart, ucl_row, width, '·', max_val, min_val, height)
    _draw_horizontal_line(chart, lcl_row, width, '·', max_val, min_val, height)
    # plot data points
    viol_idx = {idx for idx, _ in violations}
    for i, v in enumerate(values):
        row = _val_to_row(v, max_val, min_val, height)
        if 0 <= row < height:
            if i in viol_idx:
                chart[row][i] = 'X'
            else:
                chart[row][i] = '*'
    lines = []
    for r in range(height):
        val_at_row = max_val - r * (max_val - min_val) / (height - 1) if height > 1 else max_val
        lines.append(f"{''.join(chart[r])} {val_at_row:7.2f}")
    # X axis labels
    label_line = ' '.join(d[:4] for d in dates)
    bottom = '-' * (width * 2 - 1)
    footer = f"Mean: {limits['xbar']:.2f}  UCL: {limits['ucl_x']:.2f}  LCL: {limits['lcl_x']:.2f}  σ: {limits['sigma']:.2f}"
    return '\n'.join(lines) + '\n' + bottom + '\n' + label_line + '\n' + footer

def format_moving_range_chart(values, limits):
    """Create ASCII MR chart for moving ranges."""
    mrs = limits['moving_ranges']
    if not mrs:
        return "No moving ranges."
    max_val = max(mrs + [limits['ucl_mr']]) * 1.05
    min_val = 0
    if max_val == 0:
        max_val = 1
    height = 10
    width = len(mrs)
    chart = [[' ' for _ in range(width)] for _ in range(height)]
    ucl_row = _val_to_row(limits['ucl_mr'], max_val, min_val, height)
    _draw_horizontal_line(chart, ucl_row, width, '·', max_val, min_val, height)
    for i, v in enumerate(mrs):
        row = _val_to_row(v, max_val, min_val, height)
        if 0 <= row < height:
            chart[row][i] = '*'
    lines = []
    for r in range(height):
        val_at_row = max_val - r * (max_val - min_val) / (height - 1) if height > 1 else max_val
        lines.append(f"{''.join(chart[r])} {val_at_row:7.2f}")
    label_line = ' '.join(['R' + str(i+1) for i in range(width)])
    footer = f"MR UCL: {limits['ucl_mr']:.2f}  Average MR: {limits['avg_mr']:.2f}"
    return '\n'.join(lines) + '\n' + '-' * (width * 2) + '\n' + label_line + '\n' + footer

def main(argv=None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description='Bioburden I-MR control chart generator with Western Electric rules.',
        epilog='Example: python -m lab_tools.microbiology.bioburden_spc control_data.csv --date-col date --value-col cfu')
    parser.add_argument('csvfile', help='CSV file with date and bioburden columns')
    parser.add_argument('--date-col', default='date', help='Column name for dates (default: date)')
    parser.add_argument('--value-col', default='bioburden', help='Column name for bioburden values (default: bioburden)')
    parser.add_argument('--no-chart', action='store_true', help='Output statistics only, no ASCII charts')
    args = parser.parse_args(argv)

    try:
        dates, values = load_data(args.csvfile, args.date_col, args.value_col)
    except Exception as e:
        print(f"Error reading data: {e}", file=sys.stderr)
        return 1

    if len(values) < 2:
        print("At least 2 data points are required.", file=sys.stderr)
        return 1

    limits = compute_individuals_limits(values)
    violations = detect_violations(dates, values, limits)

    if not args.no_chart:
        print("--- Individual Values (I) Chart ---")
        print(format_individual_chart(dates, values, limits, violations))
        print("\n--- Moving Range (MR) Chart ---")
        print(format_moving_range_chart(values, limits))
    else:
        print(f"Mean: {limits['xbar']:.2f}")
        print(f"UCL: {limits['ucl_x']:.2f}  LCL: {limits['lcl_x']:.2f}")
        print(f"Sigma (σ): {limits['sigma']:.2f}")
        print(f"Avg Moving Range: {limits['avg_mr']:.2f}")
    if violations:
        print("\nViolations detected:")
        for _, msg in violations:
            print(f"  • {msg}")
    else:
        print("\nNo rule violations detected.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
