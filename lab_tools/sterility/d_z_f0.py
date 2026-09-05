"""D-value and Z-value calculator for sterilization validation.

Compute thermal death kinetics parameters (D-value and Z-value) from
survival-curve data, with linear regression, standard errors and confidence
intervals.

Author: Collins Amatu Gorgerat
License: MIT
"""

import argparse
import csv
import math
import statistics
import sys
from typing import List, Optional, Tuple

from scipy import stats as scipy_stats


def linear_regression(
    x: List[float], y: List[float]
) -> Tuple[float, float, float, float, float, float]:
    """Perform simple linear regression of y on x.

    Returns slope, intercept, r_squared, stderr_slope, stderr_intercept,
    resid_std_err. Raises ValueError if fewer than 2 points or x has zero
    variance.
    """
    n = len(x)
    if n < 2:
        raise ValueError("At least two data points required for regression.")

    mean_x = statistics.mean(x)
    mean_y = statistics.mean(y)

    ssxx = sum((xi - mean_x) ** 2 for xi in x)
    if ssxx == 0:
        raise ValueError("x values must have variance for regression.")

    ssxy = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    slope = ssxy / ssxx
    intercept = mean_y - slope * mean_x

    # residuals
    residuals = [yi - (intercept + slope * xi) for xi, yi in zip(x, y)]
    ss_res = sum(r ** 2 for r in residuals)
    ss_tot = sum((yi - mean_y) ** 2 for yi in y)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 1.0

    # standard errors (require at least 3 points)
    dof = n - 2
    if dof < 1:
        return slope, intercept, r_squared, None, None, None
    sigma2 = ss_res / dof
    se_slope = math.sqrt(sigma2 / ssxx)
    se_intercept = math.sqrt(sigma2 * (1 / n + mean_x ** 2 / ssxx))

    return slope, intercept, r_squared, se_slope, se_intercept, math.sqrt(sigma2)


def compute_d_value(times: List[float], logN: List[float], confidence: float = 0.95) -> dict:
    """Calculate D-value from survival curve data (time in minutes, log10(N)).

    D-value = -1 / slope.
    """
    if len(times) < 2:
        raise ValueError("Need at least two time points to compute D-value.")
    if len(times) != len(logN):
        raise ValueError("times and logN must have same length.")

    slope, intercept, r_sq, se_slope, se_intercept, _ = linear_regression(times, logN)

    if slope >= 0:
        raise ValueError(
            f"Slope is non-negative ({slope:.4f}). Data does not show a death curve; cannot compute D-value."
        )

    D = -1.0 / slope

    se_D = (1.0 / (slope ** 2)) * se_slope if se_slope is not None else None

    dof = len(times) - 2
    t_crit = scipy_stats.t.ppf(1 - (1 - confidence) / 2, dof) if dof >= 1 else None
    if se_D is not None and t_crit is not None:
        lower = D - t_crit * se_D
        upper = D + t_crit * se_D
        ci = (max(0.0, lower), max(0.0, upper))
    else:
        ci = (None, None)

    result = {
        "D_value_min": D,
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_sq,
        "se_D": se_D,
        "confidence_interval_95": ci,
        "unit": "minutes",
    }
    return result


def compute_z_value(
    temperatures: List[float], d_values: List[float], confidence: float = 0.95
) -> dict:
    """Calculate Z-value from D-values at different temperatures.

    Model: log10(D) = a - (1/Z)*T  => Z = -1 / slope.
    """
    if len(temperatures) < 2:
        raise ValueError("Need at least two temperature-Dvalue pairs to compute Z-value.")
    if len(temperatures) != len(d_values):
        raise ValueError("temperatures and d_values must have same length.")

    logD = [math.log10(d) for d in d_values]
    slope, intercept, r_sq, se_slope, se_intercept, _ = linear_regression(temperatures, logD)

    if slope >= 0:
        raise ValueError(
            f"Slope of logD vs T is non-negative ({slope:.4f}). D-value should decrease with temperature; check data."
        )

    Z = -1.0 / slope

    se_Z = (1.0 / (slope ** 2)) * se_slope if se_slope is not None else None

    dof = len(temperatures) - 2
    t_crit = scipy_stats.t.ppf(1 - (1 - confidence) / 2, dof) if dof >= 1 else None
    if se_Z is not None and t_crit is not None:
        lower = Z - t_crit * se_Z
        upper = Z + t_crit * se_Z
        ci = (max(0.0, lower), max(0.0, upper))
    else:
        ci = (None, None)

    result = {
        "Z_value_C": Z,
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_sq,
        "se_Z": se_Z,
        "confidence_interval_95": ci,
        "unit": "°C",
    }
    return result


def read_csv(filename: str) -> List[dict]:
    """Read CSV with columns: temperature (optional), time, logN."""
    rows = []
    with open(filename, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def parse_data(rows: List[dict], require_temp: bool = False) -> Tuple[List[float], List[float], List[float]]:
    """Extract lists of times, logN, and optionally temperatures.

    Returns temps, times, logN. temps may be empty if not present.
    """
    times = []
    logN = []
    temps = []
    for i, row in enumerate(rows, start=1):
        try:
            t = float(row["time"])
            n = float(row["logN"])
            times.append(t)
            logN.append(n)
        except (KeyError, ValueError) as e:
            raise ValueError(f"Row {i}: invalid time or logN: {e}")
        if "temperature" in row and row["temperature"].strip() != "":
            try:
                temp_val = float(row["temperature"])
                temps.append(temp_val)
            except ValueError:
                raise ValueError(f"Row {i}: invalid temperature value.")
        elif require_temp:
            raise ValueError(f"Row {i}: missing temperature, required for Z-value calculation.")
    if require_temp and len(temps) != len(times):
        raise ValueError("Temperature column must have same number of entries as data rows.")

    return temps, times, logN


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="D-value and Z-value calculator for sterilization validation"
    )
    subparsers = parser.add_subparsers(dest="command", help="Subcommand: dvalue or zvalue")
    subparsers.required = True

    d_parser = subparsers.add_parser("dvalue", help="Compute D-value from time-logN data")
    d_parser.add_argument("input", help="CSV file with columns: time, logN")
    d_parser.add_argument(
        "--confidence", type=float, default=0.95,
        help="Confidence level for interval (default 0.95)",
    )

    z_parser = subparsers.add_parser(
        "zvalue",
        help="Compute Z-value from multiple temperature-D-value pairs",
    )
    z_parser.add_argument(
        "input",
        help="CSV file with columns: temperature, time, logN (multiple runs at different temperatures); "
        "the tool will compute D-value per temperature then Z-value.",
    )
    z_parser.add_argument(
        "--confidence", type=float, default=0.95,
        help="Confidence level (default 0.95)",
    )
    z_parser.add_argument(
        "--dvalue-only", action="store_true",
        help="Only output per-temperature D-values (no Z-value)",
    )

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        rows = read_csv(args.input)
    except FileNotFoundError:
        print(f"Error: File not found: {args.input}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        return 1

    if args.command == "dvalue":
        try:
            temps, times, logN = parse_data(rows, require_temp=False)
            result = compute_d_value(times, logN, args.confidence)
            print("D-value Calculation Results")
            print("============================")
            print(f"D-value: {result['D_value_min']:.3f} minutes")
            print(f"Slope: {result['slope']:.6f} log10/min")
            print(f"Intercept: {result['intercept']:.4f}")
            print(f"R²: {result['r_squared']:.6f}")
            if result["se_D"] is not None:
                print(f"Standard Error of D-value: {result['se_D']:.4f} min")
                ci_low, ci_high = result["confidence_interval_95"]
                print(f"95% Confidence Interval: [{ci_low:.4f}, {ci_high:.4f}] minutes")
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    elif args.command == "zvalue":
        try:
            temps, times, logN = parse_data(rows, require_temp=True)
            temp_to_data = {}
            for ti, ni, tempi in zip(times, logN, temps):
                temp_to_data.setdefault(tempi, []).append((ti, ni))

            d_results = {}
            for temp, data_list in sorted(temp_to_data.items()):
                if len(data_list) < 2:
                    print(
                        f"Warning: Temperature {temp} has fewer than 2 time points, "
                        "skipping D-value calculation.",
                        file=sys.stderr,
                    )
                    continue
                t_vals = [d[0] for d in data_list]
                n_vals = [d[1] for d in data_list]
                try:
                    d_val = compute_d_value(t_vals, n_vals, args.confidence)
                    d_results[temp] = d_val
                except ValueError as e:
                    print(f"Warning: Temperature {temp}: {e}", file=sys.stderr)

            if not d_results:
                print("Error: No valid D-values could be computed.", file=sys.stderr)
                return 1

            print("Per-temperature D-values")
            print("--------------------------")
            for temp in sorted(d_results):
                d = d_results[temp]
                D_val = d["D_value_min"]
                print(f"Temperature {temp}°C: D = {D_val:.3f} min (R²={d['r_squared']:.4f})")

            if args.dvalue_only:
                return 0

            temps_list = sorted(d_results.keys())
            d_vals_list = [d_results[t]["D_value_min"] for t in temps_list]
            if len(temps_list) < 2:
                print(
                    "Error: Need at least two temperatures with valid D-values to compute Z-value.",
                    file=sys.stderr,
                )
                return 1

            print("\nZ-value Calculation")
            print("-------------------")
            try:
                z_result = compute_z_value(temps_list, d_vals_list, args.confidence)
                print(f"Z-value: {z_result['Z_value_C']:.2f} °C")
                print(f"Slope of log10(D) vs T: {z_result['slope']:.6f}")
                print(f"R²: {z_result['r_squared']:.4f}")
                if z_result["se_Z"] is not None:
                    ci_low, ci_high = z_result["confidence_interval_95"]
                    print(f"95% CI: [{ci_low:.2f}, {ci_high:.2f}] °C")
            except ValueError as e:
                print(f"Error computing Z-value: {e}", file=sys.stderr)
                return 1

        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
