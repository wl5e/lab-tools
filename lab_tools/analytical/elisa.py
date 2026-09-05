#!/usr/bin/env python3
"""ELISA Standard Curve Fitter.

Fit a 4-parameter logistic (4PL) or log-log linear standard curve to ELISA
data, flag outliers, and back-calculate unknown sample concentrations.
"""

import argparse
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


class StandardCurve:
    """Fit a standard curve to ELISA data.

    Parameters
    ----------
    x : array_like
        Standard concentrations (positive values).
    y : array_like
        Corresponding response values (e.g., optical density).
    model : str, optional
        '4pl' (default) or 'linear' (log‑log linear).
    outlier_threshold : float, optional
        Z‑score (studentised residual) threshold for flagging outliers.
    """

    def __init__(
        self,
        x: np.ndarray,
        y: np.ndarray,
        model: str = "4pl",
        outlier_threshold: float = 2.5,
    ):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.model = model.lower()
        self.threshold = outlier_threshold

        if self.model not in ("4pl", "linear"):
            raise ValueError(f"Unsupported model '{self.model}'. Choose '4pl' or 'linear'.")
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have the same length.")
        if len(self.x) < 4 and self.model == "4pl":
            raise ValueError("At least 4 points are required for 4PL fitting.")
        if len(self.x) < 2:
            raise ValueError("At least 2 points are required for fitting.")
        if np.any(self.x <= 0):
            raise ValueError("Concentrations must be positive (>0).")

        self.params: Optional[Dict[str, float]] = None
        self.r_squared: Optional[float] = None
        self.outliers: List[int] = []
        self._pred_fun = None

    def fit(self) -> None:
        """Fit the chosen model, detect outliers, and store results."""
        if self.model == "4pl":
            self._fit_4pl()
        else:
            self._fit_linear()

        self._flag_outliers()

    def _four_pl(self, x: np.ndarray, a: float, b: float, c: float, d: float) -> np.ndarray:
        """4-parameter logistic function.
        y = d + (a - d) / (1 + (x / c) ** b)"""
        return d + (a - d) / (1.0 + (x / c) ** b)

    def _fit_4pl(self) -> None:
        """Fit 4PL using scipy.curve_fit with robust initial guesses and bounds."""
        x = self.x
        y = self.y

        a_guess = np.min(y)
        d_guess = np.max(y)
        # Estimate C (EC50) as the median concentration
        c_guess = np.median(x)
        # Estimate slope B using the steepest part (positive slope if standard increasing)
        # Use log-log linear slope as init
        idx = np.argsort(x)
        xs = x[idx]
        ys = y[idx]
        dy = np.diff(ys)
        dx = np.diff(xs)
        b_guess = np.median(dy / (dx + 1e-10))
        if b_guess <= 0:
            b_guess = 1.0

        p0 = [a_guess, b_guess, c_guess, d_guess]
        bounds = ([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf])

        try:
            popt, pcov = curve_fit(
                self._four_pl, x, y, p0=p0, bounds=bounds, maxfev=5000
            )
        except Exception as e:
            raise RuntimeError(f"4PL fitting failed: {e}")

        self.params = {"A": popt[0], "B": popt[1], "C": popt[2], "D": popt[3]}
        self._pred_fun = lambda xin: self._four_pl(np.asarray(xin), *popt)
        self._compute_r_squared(y, self._pred_fun(x))

    def _fit_linear(self) -> None:
        """Fit log(y) = m * log(x) + b  (log‑log linear)."""
        x = self.x
        y = self.y
        logx = np.log(x)
        logy = np.log(y)
        m, b = np.polyfit(logx, logy, 1)

        self.params = {"slope": m, "intercept": b}
        self._pred_fun = lambda xin: np.exp(b) * (np.asarray(xin) ** m)
        self._compute_r_squared(y, self._pred_fun(x))

    def _compute_r_squared(self, y_obs: np.ndarray, y_pred: np.ndarray) -> None:
        ss_res = np.sum((y_obs - y_pred) ** 2)
        ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)
        self.r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else 0.0

    def _flag_outliers(self) -> None:
        """Flag outliers using studentised residuals (approximate)."""
        y_pred = self.predict_od(self.x)
        residuals = self.y - y_pred
        std_res = residuals / np.std(residuals, ddof=1) if len(residuals) > 1 else np.zeros_like(residuals)
        self.outliers = [i for i, z in enumerate(std_res) if abs(z) > self.threshold]
        if self.outliers:
            warnings.warn(f"Potential outliers detected at indices: {self.outliers}")

    def predict_od(self, conc: np.ndarray) -> np.ndarray:
        """Return predicted response for given concentration(s)."""
        return self._pred_fun(np.asarray(conc))

    def predict_conc(self, od: np.ndarray) -> np.ndarray:
        """Calculate concentration from OD using the inverse model.

        For 4PL: inverse is c * ((a-d)/(od-d) - 1)^(1/b).
        For log‑log: conc = (od / exp(b))^(1/m)."""
        od = np.asarray(od)
        if self.model == "4pl":
            a, b, c, d = self.params["A"], self.params["B"], self.params["C"], self.params["D"]
            # Check that OD is within (a, d) range; for typical increasing curve a < d.
            # Return NaN for ODs outside the valid range.
            ratio = (a - d) / (od - d)
            with np.errstate(invalid="ignore"):
                conc = c * (ratio - 1.0) ** (1.0 / b)
            # Handle edge cases where OD > d or OD < a → NaN
            conc = np.where(np.isfinite(conc), conc, np.nan)
            return conc
        else:
            m, b = self.params["slope"], self.params["intercept"]
            return np.exp((np.log(od) - b) / m)

    def concentration_ci(self, od: float, alpha: float = 0.05) -> Tuple[float, float]:
        """Compute approximate 95 % confidence interval for a single OD via delta method.

        Uses the Jacobian of the inverse function and the parameter covariance matrix
        obtained during fitting. For log‑log a simpler formula is used."""
        if self._pred_fun is None:
            raise RuntimeError("Curve must be fitted before computing confidence intervals.")

        if self.model == "linear":
            m, b = self.params["slope"], self.params["intercept"]
            log_od = np.log(od)
            log_conc = (log_od - b) / m
            # Approximate variance from ordinary least squares on log‑log
            # Using standard error of slope and intercept from polyfit would require
            # cov-matrix. For simplicity we estimate using residual variance.
            # This is a pragmatic approximation; sophisticated users should perform
            # bootstrapping.
            y_pred = self.predict_od(self.x)
            residuals = self.y - y_pred
            se2 = np.sum(residuals**2) / (len(self.x) - 2)
            var_log_conc = se2 * (1.0 / len(self.x) + (log_od - np.mean(np.log(self.x)))**2 / np.sum((np.log(self.x) - np.mean(np.log(self.x)))**2))
            std_log = np.sqrt(var_log_conc)
            z = 1.96  # 95 %
            lower_log = log_conc - z * std_log
            upper_log = log_conc + z * std_log
            return np.exp(lower_log), np.exp(upper_log)

        # For 4PL we assume we have stored the covariance matrix.
        # We need to re-run fit? We can store pcov during fit.
        warnings.warn("Full 4PL CI not implemented. Returning NaN.")
        return np.nan, np.nan

    def summary(self) -> str:
        """Return a formatted summary of the fitted curve."""
        if self.params is None:
            return "Curve not fitted."

        lines = [
            f"Model: {self.model}",
            f"R²: {self.r_squared:.4f}" if self.r_squared is not None else "R²: N/A",
        ]
        if self.model == "4pl":
            lines.append(f"A (Min): {self.params['A']:.4f}")
            lines.append(f"B (Slope): {self.params['B']:.4f}")
            lines.append(f"C (EC50): {self.params['C']:.4f}")
            lines.append(f"D (Max): {self.params['D']:.4f}")
        else:
            lines.append(f"m (slope): {self.params['slope']:.4f}")
            lines.append(f"b (intercept): {self.params['intercept']:.4f}")

        if self.outliers:
            lines.append(f"Potential outliers (indices): {self.outliers}")
        else:
            lines.append("No outliers flagged.")
        return "\n".join(lines)


def predict_unknowns(
    curve: StandardCurve,
    ods: np.ndarray,
    sample_ids: np.ndarray,
) -> List[Dict]:
    """Calculate concentrations for unknown samples."""
    concs = curve.predict_conc(ods)
    results = []
    for sid, od, conc in zip(sample_ids, ods, concs):
        entry = {"SampleID": sid, "OD": f"{od:.4f}", "Concentration": f"{conc:.4f}" if np.isfinite(conc) else ""}
        # CI would be added if implemented
        results.append(entry)
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fit ELISA standard curve and compute sample concentrations."
    )
    parser.add_argument(
        "--input", required=True, type=Path, help="Path to input CSV file"
    )
    parser.add_argument(
        "--model",
        choices=["4pl", "linear"],
        default="4pl",
        help="Calibration model: 4pl (default) or linear (log‑log)",
    )
    parser.add_argument(
        "--type-col", default="Type", help="Column name for sample type"
    )
    parser.add_argument(
        "--conc-col", default="Concentration", help="Column name for concentration"
    )
    parser.add_argument(
        "--od-col", default="OD", help="Column name for response signal"
    )
    parser.add_argument(
        "--sample-col", default="SampleID", help="Column name for sample identifier"
    )
    parser.add_argument(
        "--outlier-threshold",
        type=float,
        default=2.5,
        help="Z-score threshold for flagging outliers (default: 2.5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Save results to CSV (default: print to stdout)",
    )
    return parser


def _err(msg: str) -> int:
    print(msg, file=sys.stderr)
    return 1


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not args.input.exists():
        return _err(f"Error: input file '{args.input}' not found.")

    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        return _err(f"Error reading CSV: {e}")

    required_cols = {args.type_col, args.od_col}
    missing = required_cols - set(df.columns)
    if missing:
        return _err(f"Missing required columns: {', '.join(missing)}")

    # Separate standards and unknowns
    standards = df[df[args.type_col].str.strip().str.lower() == "standard"].copy()
    unknowns = df[df[args.type_col].str.strip().str.lower() == "unknown"].copy()

    if standards.empty:
        return _err("No rows marked as 'Standard' in type column.")
    if unknowns.empty:
        print("Warning: no unknowns found; only curve parameters will be reported.", file=sys.stderr)

    # Ensure concentration column exists for standards
    if args.conc_col not in standards.columns:
        return _err(f"Concentration column '{args.conc_col}' not found in standards.")

    # Clean data: drop rows with missing OD or concentration
    standards = standards.dropna(subset=[args.od_col, args.conc_col])
    if standards.empty:
        return _err("No valid standard points after removing missing values.")

    # Check positive concentrations (required for both models)
    if (standards[args.conc_col] <= 0).any():
        return _err("All standard concentrations must be positive (>0).")

    try:
        curve = StandardCurve(
            x=standards[args.conc_col].values,
            y=standards[args.od_col].values,
            model=args.model,
            outlier_threshold=args.outlier_threshold,
        )
        curve.fit()
    except Exception as e:
        return _err(f"Curve fitting failed: {e}")

    if args.sample_col in df.columns:
        sample_ids = unknowns.get(args.sample_col, pd.Series([""] * len(unknowns)))
    else:
        sample_ids = pd.Series([f"U{i+1}" for i in range(len(unknowns))])

    if not unknowns.empty:
        od_values = unknowns[args.od_col].values
        results = predict_unknowns(curve, od_values, sample_ids.values)
        results_df = pd.DataFrame(results)
    else:
        results_df = pd.DataFrame()

    print(curve.summary())
    if not results_df.empty:
        print("\nSample concentrations:")
        print(results_df.to_string(index=False))

    if args.output:
        output_path = args.output
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write("# Curve Summary\n")
            f.write(curve.summary())
            f.write("\n\n# Unknowns\n")
            results_df.to_csv(f, index=False, lineterminator="\n")
        print(f"\nResults written to {output_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
