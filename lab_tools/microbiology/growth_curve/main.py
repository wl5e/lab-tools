"""Command-line interface for growth curve analysis."""
import argparse
import sys
import os
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from .analyzer import fit_growth_curve, MODEL_MAP


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Fit growth models to OD measurements and export parameters."
    )
    parser.add_argument(
        "input",
        help="Path to input CSV file.",
    )
    parser.add_argument(
        "--well-col",
        default="Well",
        help="Column name that identifies the well/replicate (default: 'Well').",
    )
    parser.add_argument(
        "--time-col",
        default="Time",
        help="Column name for time (default: 'Time').",
    )
    parser.add_argument(
        "--od-col",
        default="OD",
        help="Column name for optical density (default: 'OD').",
    )
    parser.add_argument(
        "--model",
        choices=list(MODEL_MAP.keys()) + ["all"],
        default="all",
        help="Model to fit: 'logistic', 'gompertz', or 'all' (both, default).",
    )
    parser.add_argument(
        "--output",
        default="growth_parameters.csv",
        help="Output CSV file path (default: growth_parameters.csv).",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate plots of fitted curves (requires matplotlib).",
    )
    parser.add_argument(
        "--plot-dir",
        default="./plots",
        help="Directory to save plot images (default: ./plots).",
    )
    parser.add_argument(
        "--delimiter",
        default=",",
        help="CSV delimiter (default: ',').",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    # Load data
    try:
        df = pd.read_csv(args.input, sep=args.delimiter)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        return 1

    required_cols = [args.well_col, args.time_col, args.od_col]
    for col in required_cols:
        if col not in df.columns:
            print(f"Missing required column '{col}' in input CSV.", file=sys.stderr)
            return 1

    # Drop rows with missing essential values
    df = df.dropna(subset=required_cols).copy()
    df[args.time_col] = pd.to_numeric(df[args.time_col], errors="coerce")
    df[args.od_col] = pd.to_numeric(df[args.od_col], errors="coerce")
    df = df.dropna(subset=[args.time_col, args.od_col])

    models_to_fit = list(MODEL_MAP.keys()) if args.model == "all" else [args.model]

    # Group by well
    grouped = df.groupby(args.well_col)
    results = []

    for well, group in grouped:
        t = group[args.time_col].values
        y = group[args.od_col].values
        if len(t) < 4:
            print(f"Well '{well}' has fewer than 4 points; skipping.", file=sys.stderr)
            continue

        for mdl in models_to_fit:
            try:
                params = fit_growth_curve(t, y, model=mdl)
            except Exception as e:
                print(f"Fit error for well '{well}', model {mdl}: {e}", file=sys.stderr)
                params = {
                    "A": float("nan"),
                    "mu": float("nan"),
                    "lag": float("nan"),
                    "r_squared": float("nan"),
                    "model": mdl,
                    "doubling_time": float("nan"),
                }
            results.append({
                "well": well,
                "model": params["model"],
                "lag_time": params["lag"],
                "max_growth_rate": params["mu"],
                "max_OD": params["A"],
                "doubling_time": params["doubling_time"],
                "r_squared": params["r_squared"],
            })

    if not results:
        print("No valid results to write.", file=sys.stderr)
        return 1

    # Write output CSV
    fieldnames = ["well", "model", "lag_time", "max_growth_rate",
                  "max_OD", "doubling_time", "r_squared"]
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"Parameters written to {args.output}")

    # Optional plotting
    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib is required for plotting. Install it with: pip install matplotlib", file=sys.stderr)
            return 1

        plot_dir = Path(args.plot_dir)
        plot_dir.mkdir(parents=True, exist_ok=True)
        for well, group in grouped:
            t = group[args.time_col].values
            y = group[args.od_col].values
            plt.figure()
            plt.scatter(t, y, label="Data", s=20)
            t_smooth = np.linspace(t.min(), t.max(), 200)
            for mdl in models_to_fit:
                try:
                    params = fit_growth_curve(t, y, model=mdl)
                    func = MODEL_MAP[mdl]
                    y_fit = func(t_smooth, params["A"], params["mu"], params["lag"])
                    plt.plot(t_smooth, y_fit, label=mdl.capitalize())
                except Exception:
                    pass
            plt.xlabel(args.time_col)
            plt.ylabel(args.od_col)
            plt.title(f"Well {well}")
            plt.legend()
            plt.tight_layout()
            out_plot = plot_dir / f"well_{well}.png"
            plt.savefig(out_plot, dpi=150)
            plt.close()
            print(f"Plot saved to {out_plot}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
