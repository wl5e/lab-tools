#!/usr/bin/env python3
"""
GMP-compliant MPN (Most Probable Number) calculator using maximum likelihood estimation.

Computes the concentration of microorganisms (per mL, per 100 mL, etc.) from serial dilution
tube count data with equal replicates, and provides 95% confidence intervals based on the
Poisson model.
"""

import argparse
import json
import math
import sys
from typing import Dict, List, Optional, Union


def compute_mpn(
    tubes: int,
    volumes: List[float],
    positive_counts: List[int],
    per_unit: float = 100.0
) -> Dict[str, Optional[Union[float, str]]]:
    """
    Calculate MPN and 95% confidence interval using MLE.

    Parameters
    ----------
    tubes : int
        Number of replicate tubes per dilution.
    volumes : list of float
        Sample volume per tube for each dilution (same unit).
    positive_counts : list of int
        Number of positive tubes for each dilution, same order as volumes.
    per_unit : float, optional
        Factor to convert concentration to desired reporting unit.
        Default 100 gives MPN per 100 mL.

    Returns
    -------
    dict
        Keys: 'mpn', 'lower_ci', 'upper_ci', 'method', 'unit_str', 'warning'.
        Values are floats or None for boundary cases. 'method' is always 'mle'.

    Raises
    ------
    ValueError
        If inputs are inconsistent or out of allowed range.
    """
    n = tubes
    v = list(volumes)
    p = list(positive_counts)

    if len(v) != len(p):
        raise ValueError("Volumes and positive counts must have the same length.")
    if n <= 0:
        raise ValueError("Number of tubes must be a positive integer.")
    for i, pi in enumerate(p):
        if pi < 0 or pi > n:
            raise ValueError(f"Positive count for dilution {i} ({pi}) out of range 0..{n}.")
    for vi in v:
        if vi <= 0.0:
            raise ValueError(f"Volume values must be positive, got {vi}.")

    unit_str = f"per {per_unit} mL"

    # --------------------------------------------------------------------
    # All tubes negative -> MPN < detection limit, one‑sided upper CI.
    # --------------------------------------------------------------------
    if all(pi == 0 for pi in p):
        total_vol = n * sum(v)
        if total_vol == 0:
            raise ValueError("Total volume examined is zero.")
        # Exact one‑sided 95% upper limit: P(all neg) = 0.05
        lam_upper = -math.log(0.05) / total_vol  # lambda/unit vol
        return {
            "mpn": 0.0,
            "lower_ci": 0.0,
            "upper_ci": lam_upper * per_unit,
            "method": "mle",
            "unit_str": unit_str,
            "warning": "All tubes negative. Upper confidence limit computed using one-sided exact probability.",
        }

    # --------------------------------------------------------------------
    # All tubes positive -> MPN > max measurable, one‑sided lower CI.
    # --------------------------------------------------------------------
    if all(pi == n for pi in p):
        def all_pos_prob(lam: float) -> float:
            """Probability all tubes positive at concentration lambda."""
            prob = 1.0
            for vi in v:
                prob *= (1.0 - math.exp(-lam * vi)) ** n
            return prob

        # Bisection to solve all_pos_prob(lam) = 0.05
        lo = 0.0
        hi = 1000.0
        # Increase hi until prob >= 0.05
        while all_pos_prob(hi) < 0.05:
            hi *= 2
            if hi > 1e12:
                raise RuntimeError("Cannot bracket root for all‑positive lower CI.")
        for _ in range(80):
            mid = (lo + hi) / 2
            if all_pos_prob(mid) >= 0.05:
                ho = mid  # keep the bracket
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-10:
                break
        lam_lower = (lo + hi) / 2
        return {
            "mpn": None,
            "lower_ci": lam_lower * per_unit,
            "upper_ci": None,
            "method": "mle",
            "unit_str": unit_str,
            "warning": "All tubes positive. MPN > max measurable. Lower confidence limit provided.",
        }

    # --------------------------------------------------------------------
    # General case: solve score equation via bisection.
    # --------------------------------------------------------------------
    def score(lam: float) -> float:
        """Score function S(λ) = Σ v_i * (p_i/(1 - e^{-λ v_i}) - n)."""
        total = 0.0
        for vi, pi in zip(v, p):
            if pi == 0:
                total += vi * (-n)  # contribution = -n*vi
            else:
                # Handle extreme lambda*vi values to avoid under/overflow
                lv = lam * vi
                if lv > 700.0:   # exp(-700) ~ 9.8e-305 ≈ 0
                    denom = 1.0
                else:
                    expn = math.exp(-lv)
                    denom = 1.0 - expn
                    if denom == 0.0:  # lv near zero, approximate denominator
                        denom = lv
                total += vi * (pi / denom - n)
        return total

    # Find a lower bound lam_low where score is large positive.
    lam_low = 1e-12
    # Increase lam_low if needed until score > 0.
    for _ in range(100):
        s_low = score(lam_low)
        if s_low > 0:
            break
        lam_low /= 10.0
        if lam_low < 1e-300:
            lam_low = 1e-300
            break
    else:
        raise RuntimeError("Cannot find lower bound with positive score.")

    # Find an upper bound lam_high where score <= 0.
    lam_high = 1.0
    for _ in range(200):
        s_high = score(lam_high)
        if s_high <= 0:
            break
        lam_high *= 2
        if lam_high > 1e10:
            raise RuntimeError("Cannot find upper bound with non‑positive score.")
    else:
        raise RuntimeError("Failed to bracket root.")

    # Bisection until interval width < tolerance.
    tol = 1e-10
    for _ in range(120):
        mid = (lam_low + lam_high) / 2
        s_mid = score(mid)
        if s_mid == 0.0:
            lam_hat = mid
            break
        if s_mid > 0:
            lam_low = mid
        else:
            lam_high = mid
        if lam_high - lam_low < tol:
            break
    else:
        raise RuntimeError("MPN root finding did not converge.")
    lam_hat = (lam_low + lam_high) / 2

    # --------------------------------------------------------------------
    # Fisher information & asymptotic confidence interval.
    # --------------------------------------------------------------------
    I = 0.0
    for vi in v:
        lv = lam_hat * vi
        if lv > 700.0:
            exp_neg = 0.0
            one_minus = 1.0
        else:
            exp_neg = math.exp(-lv)
            one_minus = 1.0 - exp_neg
        if one_minus == 0.0:  # degenerate (lam_hat=0, but handled already)
            continue
        I += n * vi * vi * exp_neg / one_minus

    if I <= 0.0:
        # Should not happen for general case, but guard.
        se = 0.0
    else:
        se = 1.0 / math.sqrt(I)

    z = 1.96  # 95% confidence
    # Asymptotic CI on the log scale (lambda > 0): keeps the lower bound
    # positive and matches standard MPN confidence limits better than a
    # Wald interval on the raw scale (which clamps the lower bound to 0).
    se_log = se / lam_hat if lam_hat > 0 else 0.0
    lam_lower = lam_hat * math.exp(-z * se_log)
    lam_upper = lam_hat * math.exp(z * se_log)

    mpn_val = lam_hat * per_unit
    lower_ci = lam_lower * per_unit
    upper_ci = lam_upper * per_unit

    return {
        "mpn": mpn_val,
        "lower_ci": lower_ci,
        "upper_ci": upper_ci,
        "method": "mle",
        "unit_str": unit_str,
        "warning": None,
    }


# ---------------------------------------------------------------------------
# Output formatting helpers
# ---------------------------------------------------------------------------
def format_result(result: dict) -> str:
    """Format the result dict as a human‑readable string."""
    unit = result["unit_str"]
    mpn_fmt = f"{result['mpn']:.2f}" if result["mpn"] is not None else "> max measurable"
    low_fmt = f"{result['lower_ci']:.2f}" if result["lower_ci"] is not None else "N/A"
    high_fmt = f"{result['upper_ci']:.2f}" if result["upper_ci"] is not None else "N/A"
    lines = [
        f"MPN ({unit}): {mpn_fmt}",
        f"95% CI: [{low_fmt}, {high_fmt}]",
        f"Method: {result['method'].upper()}",
    ]
    if result.get("warning"):
        lines.append(f"Warning: {result['warning']}")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="GMP-compliant MPN calculator using maximum likelihood estimation."
    )
    parser.add_argument(
        "-n", "--tubes", type=int, required=True,
        help="Number of replicate tubes per dilution.",
    )
    parser.add_argument(
        "-v", "--volumes", nargs="+", type=float, required=True,
        help="Sample volumes per tube (space separated). e.g. 10 1 0.1",
    )
    parser.add_argument(
        "-p", "--positive", nargs="+", type=int, required=True,
        help="Positive tube counts, same order as volumes.",
    )
    parser.add_argument(
        "-u", "--unit", type=float, default=100.0,
        help="Unit factor for reporting (default: 100 for per 100 mL).",
    )
    parser.add_argument(
        "-o", "--output", choices=["text", "json"], default="text",
        help="Output format.",
    )
    args = parser.parse_args(argv)

    try:
        result = compute_mpn(args.tubes, args.volumes, args.positive, args.unit)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except RuntimeError as e:
        print(f"Runtime error: {e}", file=sys.stderr)
        return 2

    if args.output == "json":
        json_result = {}
        for k, v in result.items():
            if v is None:
                json_result[k] = None
            else:
                json_result[k] = v
        print(json.dumps(json_result, indent=2))
    else:
        print(format_result(result))

    return 0


if __name__ == "__main__":
    sys.exit(main())