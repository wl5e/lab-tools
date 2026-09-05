"""Command-line interface for endotoxin limit and MVD calculations."""

import argparse
import sys
from typing import List, Optional

from lab_tools.sterility.endotoxin.core import (
    calc_endotoxin_limit_product,
    calc_endotoxin_limit_volume,
    calc_mvd,
)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Endotoxin limit and MVD calculator per USP <85>",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py product --dose-mg 10 --concentration 2
  python main.py product --dose-mg 10 --concentration 2 --k 0.2
  python main.py product --dose-mg 10        (returns EU/mg)
  python main.py product --dose-vol 5
  python main.py mvd --endotoxin-limit 0.5 --lysate-sensitivity 0.125
        """,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # product subcommand
    product_parser = subparsers.add_parser(
        "product", help="Calculate endotoxin limit for a drug product"
    )
    dose_group = product_parser.add_mutually_exclusive_group(required=True)
    dose_group.add_argument(
        "--dose-mg", type=float, metavar="MG",
        help="Max mass dose (mg active / kg body weight / hour)",
    )
    dose_group.add_argument(
        "--dose-vol", type=float, metavar="ML",
        help="Max volume dose (mL product / kg body weight / hour)",
    )
    product_parser.add_argument(
        "--concentration", type=float, metavar="MG/ML",
        help="Concentration of active substance (mg/mL). If provided, output in EU/mL.",
    )
    product_parser.add_argument(
        "--k", type=float, default=None,
        help="K factor in EU/kg (overrides route default)",
    )
    product_parser.add_argument(
        "--route", type=str, default="IV",
        choices=["IV", "intrathecal"],
        help="Administration route (default: IV, sets K=5 for IV, 0.2 for intrathecal if --k not given)",
    )

    # mvd subcommand
    mvd_parser = subparsers.add_parser("mvd", help="Calculate Maximum Valid Dilution")
    mvd_parser.add_argument(
        "--endotoxin-limit", type=float, required=True, metavar="EU/ML",
        help="Endotoxin limit of product in EU/mL",
    )
    mvd_parser.add_argument(
        "--lysate-sensitivity", type=float, required=True, metavar="EU/ML",
        help="Lysate sensitivity (λ) in EU/mL",
    )

    args = parser.parse_args(argv)

    try:
        if args.command == "product":
            # Determine K factor
            if args.k is not None:
                k = args.k
            elif args.route == "IV":
                k = 5.0
            elif args.route == "intrathecal":
                k = 0.2
            else:
                raise ValueError(f"Unknown route '{args.route}'")

            if args.dose_mg is not None:
                result = calc_endotoxin_limit_product(
                    dose_mg_per_kg_hour=args.dose_mg,
                    concentration_mg_per_ml=args.concentration,
                    k=k,
                )
            else:  # dose_vol provided
                result = calc_endotoxin_limit_volume(
                    dose_ml_per_kg_hour=args.dose_vol,
                    k=k,
                )

            # Display results
            print(f"Endotoxin Limit: {result['endotoxin_limit']} {result['unit']}")
            if "K" in result:
                print(f"K factor used: {result['K']} EU/kg")
            if "dose_mg_per_kg_hour" in result:
                print(f"Max dose (mass): {result['dose_mg_per_kg_hour']} mg/kg/h")
            if "dose_ml_per_kg_hour" in result:
                print(f"Max dose (volume): {result['dose_ml_per_kg_hour']} mL/kg/h")
            if result.get("concentration_mg_per_ml") is not None:
                print(f"Concentration: {result['concentration_mg_per_ml']} mg/mL")

        elif args.command == "mvd":
            mvd = calc_mvd(args.endotoxin_limit, args.lysate_sensitivity)
            print(f"Maximum Valid Dilution (MVD): {mvd:.2f}")
            print(f"Endotoxin Limit: {args.endotoxin_limit} EU/mL")
            print(f"Lysate Sensitivity: {args.lysate_sensitivity} EU/mL")
        else:
            print("Unknown command", file=sys.stderr)
            return 2
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
