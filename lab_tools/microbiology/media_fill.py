#!/usr/bin/env python3
"""
Media Fill Contamination Analyzer
=================================
Core analysis logic and CLI for media fill contamination data per USP <1116>.
"""
import argparse
import csv
import math
import sys
from typing import List, Optional, Tuple

class BatchResult:
    def __init__(self, batch_id: str, n: int, x: int, rate: float, ci_lower: float, ci_upper: float):
        self.batch_id = batch_id
        self.n = n
        self.x = x
        self.rate = rate
        self.ci_lower = ci_lower
        self.ci_upper = ci_upper

    def __str__(self):
        return (f"{self.batch_id}: n={self.n}, contaminated={self.x}, "
                f"rate={self.rate:.5f} (CI: {self.ci_lower:.5f} - {self.ci_upper:.5f})")

# Standard normal cumulative distribution function
norm_cdf = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))

# Mapping: significance level α -> two-sided Z-score
ALPHA_TO_Z = {
    0.001: 3.2905,
    0.002: 2.8782,
    0.005: 2.8070,
    0.01: 2.5758,
    0.02: 2.3263,
    0.05: 1.9600,
    0.1: 1.6449,
    0.2: 1.2816
}

def get_z(alpha: float) -> float:
    """Return the two-sided Z-score for a given alpha."""
    # Find exact match or nearest key; raise if not found
    if alpha in ALPHA_TO_Z:
        return ALPHA_TO_Z[alpha]
    # Interpolation? Not needed; raise error
    raise ValueError(
        f"Alpha={alpha} not in precomputed levels. Supported: {sorted(ALPHA_TO_Z.keys())}"
    )

def wilson_ci(n: int, x: int, z: float) -> Tuple[float, float, float]:
    """
    Compute Wilson score confidence interval for a binomial proportion.

    Args:
        n: total trials
        x: number of successes (contaminated units)
        z: two-sided Z-value for desired confidence level

    Returns:
        (rate, lower_bound, upper_bound)
    """
    if n <= 0:
        raise ValueError("n must be positive")
    if x < 0 or x > n:
        raise ValueError("x must be between 0 and n")
    p_hat = x / n
    z2 = z * z
    denominator = 1 + z2 / n
    center = (p_hat + z2 / (2 * n)) / denominator
    margin = z * math.sqrt((p_hat * (1 - p_hat) + z2 / (4 * n)) / n) / denominator
    lower = max(0, center - margin)
    upper = min(1, center + margin)
    return p_hat, lower, upper

def parse_csv(file_path: str) -> List[dict]:
    """
    Read media fill CSV and return a list of rows (order preserved).
    Required columns: batch, n, contaminated
    Optional: date (used for sorting if present)
    """
    rows = []
    with open(file_path, newline='') as f:
        reader = csv.DictReader(f)
        required = {'batch', 'n', 'contaminated'}
        if reader.fieldnames is None:
            raise ValueError("CSV file does not contain a header")
        if not required.issubset(set(reader.fieldnames)):
            missing = required - set(reader.fieldnames)
            raise ValueError(f"Missing required columns: {missing}")
        for row in reader:
            try:
                n = int(row['n'])
                x = int(row['contaminated'])
            except (ValueError, KeyError) as e:
                raise ValueError(f"Invalid numeric field in row {row}: {e}")
            if n <= 0:
                raise ValueError(f"n must be positive (row: {row})")
            if x < 0 or x > n:
                raise ValueError(f"contaminated must be 0 <= x <= n (row: {row})")
            rows.append({
                'batch': row['batch'],
                'n': n,
                'x': x,
                'date': row.get('date', '')
            })
    # Sort by date if date column present and contains values
    if rows and 'date' in rows[0] and any(r['date'] for r in rows):
        try:
            rows.sort(key=lambda r: r['date'])
        except Exception as e:
            raise ValueError(f"Cannot sort by date column: {e}")
    return rows

def analyze_batches(rows: List[dict], alpha: float) -> List[BatchResult]:
    """Compute per-batch rates and confidence intervals."""
    z = get_z(alpha)
    results = []
    for row in rows:
        rate, lo, hi = wilson_ci(row['n'], row['x'], z)
        results.append(BatchResult(row['batch'], row['n'], row['x'], rate, lo, hi))
    return results

def aggregate(results: List[BatchResult], alpha: float) -> BatchResult:
    """Compute overall contamination rate and CI from all batches."""
    total_n = sum(r.n for r in results)
    total_x = sum(r.x for r in results)
    z = get_z(alpha)
    rate, lo, hi = wilson_ci(total_n, total_x, z)
    return BatchResult("TOTAL", total_n, total_x, rate, lo, hi)

def cochran_armitage_trend(results: List[BatchResult]) -> Tuple[float, float]:
    """
    Perform Cochran-Armitage trend test for binomial proportions.
    Scores are the ordered indices (0,1,2,...).
    Returns (Z_statistic, one_sided_p_value).
    """
    if len(results) < 2:
        return 0.0, 0.5  # not enough data
    n_vals = [r.n for r in results]
    x_vals = [r.x for r in results]
    scores = list(range(len(results)))
    total_n = sum(n_vals)
    total_x = sum(x_vals)
    if total_x == 0 or total_x == total_n:
        # No variation; trend undefined – return neutral
        return 0.0, 0.5
    p_hat = total_x / total_n
    # Expected counts under null
    exp = [n * p_hat for n in n_vals]
    # Numerator
    num = sum(s * (x - e) for s, x, e in zip(scores, x_vals, exp))
    # Variance
    sum_ns = sum(n * s for n, s in zip(n_vals, scores))
    sum_ns2 = sum(n * s * s for n, s in zip(n_vals, scores))
    var = p_hat * (1 - p_hat) * (sum_ns2 - (sum_ns ** 2) / total_n)
    if var <= 0:
        return 0.0, 0.5
    Z = num / math.sqrt(var)
    # One-sided p-value for increasing trend (Z > 0)
    p_val = 1 - norm_cdf(Z)
    return Z, p_val

def compute_analysis(input_path: str, output_path: Optional[str], alpha: float) -> str:
    """
    Run full analysis, return string summary, and optionally write CSV report.
    """
    rows = parse_csv(input_path)
    results = analyze_batches(rows, alpha)
    overall = aggregate(results, alpha)
    Z, p_val = cochran_armitage_trend(results)

    lines = ["=== Media Fill Contamination Analysis ===", ""]
    for res in results:
        lines.append(str(res))
    lines.append("")
    lines.append(f"Aggregate: {overall}")
    lines.append(f"Cochran-Armitage trend test (one-sided): Z = {Z:.4f}, p = {p_val:.4f}")
    if p_val < 0.05:
        lines.append("  => Significant increasing trend (p < 0.05).")
    else:
        lines.append("  => No significant increasing trend.")
    summary = "\n".join(lines)

    if output_path:
        _write_report(results, overall, Z, p_val, output_path)
    return summary

def _write_report(results: List[BatchResult], overall: BatchResult,
                  Z: float, p_val: float, path: str):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['batch', 'n', 'contaminated', 'rate', 'ci_lower', 'ci_upper'])
        for r in results:
            writer.writerow([r.batch_id, r.n, r.x, f"{r.rate:.6f}",
                             f"{r.ci_lower:.6f}", f"{r.ci_upper:.6f}"])
        writer.writerow(['TOTAL', overall.n, overall.x, f"{overall.rate:.6f}",
                         f"{overall.ci_lower:.6f}", f"{overall.ci_upper:.6f}"])
        writer.writerow([])
        writer.writerow(['trend_test_Z', Z, 'one_sided_p', p_val])


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description='Analyze media fill contamination data per USP <1116>.'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Path to input CSV file (batch, n, contaminated [,date]).')
    parser.add_argument('--output', '-o', default=None,
                        help='Optional path for result CSV report.')
    parser.add_argument('--alpha', type=float, default=0.05,
                        help='Significance level for confidence intervals (default 0.05).')
    args = parser.parse_args(argv)

    try:
        summary = compute_analysis(args.input, args.output, args.alpha)
        print(summary)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
