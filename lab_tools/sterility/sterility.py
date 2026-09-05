"""Sterility test sample size calculations per USP <71> and EP 2.6.1.

Author: Collins Amatu Gorgerat
License: MIT
"""

import argparse
import sys
from typing import List, Optional


def calculate_sample_size(batch_size, product_type, pharmacopeia="USP"):
    """Calculate the minimum number of units to test for sterility.

    Args:
        batch_size (int): Total number of container units in the batch.
        product_type (str): One of 'liquid', 'solid', 'ophthalmic'.
        pharmacopeia (str): 'USP' or 'EP'.

    Returns:
        tuple: (actual_sample_size, warning_string_or_None)

    Raises:
        ValueError: If inputs are invalid.
    """
    if not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError(f"Batch size must be a positive integer, got {batch_size}")

    valid_types = {"liquid", "solid", "ophthalmic"}
    if product_type not in valid_types:
        raise ValueError(
            f"Invalid product type '{product_type}', must be one of {valid_types}"
        )

    if pharmacopeia not in ("USP", "EP"):
        raise ValueError("Pharmacopeia must be 'USP' or 'EP'")

    # Common tables (USP & EP share the same harmonized rules for these categories)
    if product_type == "liquid":
        if batch_size < 100:
            required = max(4, int(batch_size * 0.1 + 1e-9))  # 10% or 4, whichever greater
        elif batch_size <= 500:
            required = 10
        else:  # >500
            required = min(20, int(batch_size * 0.02 + 1e-9))  # 2% or 20, whichever less
    elif product_type == "solid":
        if batch_size <= 200:
            required = max(2, int(batch_size * 0.05 + 1e-9))  # 5% or 2, whichever greater
        else:
            required = 10
    else:  # ophthalmic
        # For ophthalmic and other noninjectable preparations, the same as solid rule applies
        if batch_size <= 200:
            required = max(2, int(batch_size * 0.05 + 1e-9))
        else:
            required = 10

    # Never test more units than exist in the batch
    actual = min(required, batch_size)
    warning = None
    if required > batch_size:
        warning = (
            f"Required sample size ({required}) exceeds batch size; "
            f"using all {batch_size} units."
        )

    return actual, warning


def calculate_volume_per_container(product_type, container_volume_ml):
    """Determine the volume to be withdrawn from each container for sterility testing.

    Args:
        product_type (str): 'liquid' or 'ophthalmic'.
        container_volume_ml (float): Nominal fill volume per container (mL).

    Returns:
        tuple: (volume_to_test_ml, description)

    Raises:
        ValueError: If product type is not applicable or container volume invalid.
    """
    if container_volume_ml <= 0:
        raise ValueError("Container volume must be positive.")

    if product_type == "liquid":
        if container_volume_ml <= 100:
            return container_volume_ml, "entire contents"
        else:
            return 10.0, "10 mL"
    elif product_type == "ophthalmic":
        if container_volume_ml <= 5:
            return container_volume_ml, "entire contents"
        else:
            return 5.0, "5 mL"
    else:
        raise ValueError(
            "Volume per container is only applicable for 'liquid' or 'ophthalmic' products."
        )


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Calculate sterility test sample size per USP/EP."
    )
    parser.add_argument(
        "--batch-size", type=int, required=True,
        help="Total number of units in the batch",
    )
    parser.add_argument(
        "--product-type", choices=["liquid", "solid", "ophthalmic"], required=True,
        help="Type of product",
    )
    parser.add_argument(
        "--pharmacopeia", choices=["USP", "EP"], default="USP",
        help="Pharmacopeia (default: USP)",
    )
    parser.add_argument(
        "--container-volume", type=float,
        help="Container fill volume in mL (optional, for liquid/ophthalmic products)",
    )
    args = parser.parse_args(argv)

    # Validate batch size early
    if args.batch_size <= 0:
        print("Error: Batch size must be a positive integer.", file=sys.stderr)
        return 1

    try:
        sample_size, warning = calculate_sample_size(
            args.batch_size, args.product_type, args.pharmacopeia
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print(f"Product type: {args.product_type}")
    print(f"Pharmacopeia: {args.pharmacopeia}")
    print(f"Batch size: {args.batch_size}")
    print(f"Minimum number of units to test: {sample_size}")

    if warning:
        print(f"WARNING: {warning}")

    if args.container_volume is not None:
        if args.product_type not in ("liquid", "ophthalmic"):
            print(
                "Warning: container volume only relevant for liquid/ophthalmic products; ignoring."
            )
        else:
            if args.container_volume <= 0:
                print("Error: Container volume must be positive.", file=sys.stderr)
                return 1
            try:
                vol, desc = calculate_volume_per_container(
                    args.product_type, args.container_volume
                )
            except ValueError as e:
                print(f"Error: {e}", file=sys.stderr)
                return 1

            total_vol = vol * sample_size
            print(f"Volume per container to test: {desc} ({vol} mL)")
            print(f"Total volume needed for each medium: {total_vol} mL")

    return 0


if __name__ == "__main__":
    sys.exit(main())
