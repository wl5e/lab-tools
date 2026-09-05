#!/usr/bin/env python3
"""
PCR annealing temperature calculator.

Supports:
- Wallace rule (2°C per A/T, 4°C per G/C) for oligos ≤14 nt
- SantaLucia unified nearest-neighbour parameters with salt correction
"""

import argparse
import math
import sys

# SantaLucia (1998) unified DNA nearest‑neighbour parameters
# ΔH in kcal/mol, ΔS in cal/(mol·K)
NN_PARAMETERS = {
    "AA": ( -7.9, -22.2 ),
    "TT": ( -7.9, -22.2 ),
    "AT": ( -7.2, -20.4 ),
    "TA": ( -7.2, -21.3 ),
    "AC": ( -8.4, -22.4 ),
    "TG": ( -8.4, -22.4 ),
    "CA": ( -8.5, -22.7 ),
    "GT": ( -8.5, -22.7 ),
    "AG": ( -7.8, -21.0 ),
    "TC": ( -7.8, -21.0 ),
    "GA": ( -8.2, -22.2 ),
    "CT": ( -8.2, -22.2 ),
    "CG": (-10.6, -27.2 ),
    "GC": ( -9.8, -24.4 ),
    "GG": ( -8.0, -19.9 ),
    "CC": ( -8.0, -19.9 ),
}

# Gas constant in cal/(mol·K)
R = 1.987


def validate_sequence(seq: str) -> str:
    """Raise ValueError if sequence contains invalid characters."""
    if not seq:
        raise ValueError("Sequence must be non‑empty")
    seq = seq.upper()
    allowed = {'A', 'C', 'G', 'T'}
    if any(ch not in allowed for ch in seq):
        invalid = [ch for ch in seq if ch not in allowed]
        raise ValueError(f"Invalid nucleotides found: {invalid}")
    return seq


def wallace_tm(sequence: str) -> float:
    """
    Calculate Ta using the Wallace rule: Tm = 2*(A+T) + 4*(G+C).
    Valid for oligos up to 14 nt; longer sequences produce a warning.
    """
    seq = validate_sequence(sequence)
    if len(seq) > 14:
        print("Warning: Wallace rule is inaccurate for sequences longer than 14 nt",
              file=sys.stderr)
    a_t = seq.count('A') + seq.count('T')
    g_c = seq.count('G') + seq.count('C')
    return 2 * a_t + 4 * g_c


def santalucia_tm(sequence: str, primer_conc_nM: float = 500,
                  salt_conc_mM: float = 50) -> float:
    """
    Calculate Tm using SantaLucia unified NN parameters.

    Parameters
    ----------
    sequence : str
        DNA sequence (5'→3', case‑insensitive).
    primer_conc_nM : float
        Total primer (strand) concentration in nM.
    salt_conc_mM : float
        Monovalent cation (Na⁺) concentration in mM.

    Returns
    -------
    float : Tm in °C, rounded to 1 decimal place.
    """
    seq = validate_sequence(sequence)
    if len(seq) < 3:
        raise ValueError("SantaLucia method requires at least 3 nucleotides")

    # Sum nearest‑neighbour pairs
    delta_h = 0.0
    delta_s = 0.0
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        try:
            h, s = NN_PARAMETERS[pair]
        except KeyError:
            # Should never happen after validation
            raise ValueError(f"Invalid dinucleotide pair: {pair}")
        delta_h += h
        delta_s += s

    # Convert concentration to M
    conc_M = primer_conc_nM * 1e-9
    salt_M = salt_conc_mM * 1e-3

    # Base Tm in Kelvin at 1 M NaCl
    if delta_s <= 0 and conc_M > 0:
        # Standard equation
        denom = delta_s + R * math.log(conc_M / 4)
        if denom == 0:
            raise ValueError("Denominator is zero – cannot calculate Tm")
        tm_kelvin = (delta_h * 1000) / denom
    else:
        raise ValueError("Invalid thermodynamic parameters or concentration")

    # Salt correction (SantaLucia 1998)
    if salt_M > 0:
        tm_celsius = tm_kelvin - 273.15 + 16.6 * math.log10(salt_M)
    else:
        tm_celsius = tm_kelvin - 273.15

    return round(tm_celsius, 1)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Calculate PCR primer annealing temperature using Wallace or SantaLucia methods."
    )
    parser.add_argument(
        "sequence",
        type=str,
        help="DNA sequence (5' → 3'), e.g. ATGCATGCATGACGT"
    )
    parser.add_argument(
        "-m", "--method",
        choices=["wallace", "santalucia", "both"],
        default="santalucia",
        help="Thermodynamic model (default: santalucia)"
    )
    parser.add_argument(
        "-c", "--conc",
        type=float,
        default=500.0,
        help="Primer concentration in nM (default: 500)"
    )
    parser.add_argument(
        "-s", "--salt",
        type=float,
        default=50.0,
        help="Monovalent cation concentration in mM (default: 50)"
    )

    args = parser.parse_args(argv)

    # Validate early to give a clean error
    try:
        validate_sequence(args.sequence)
    except ValueError as e:
        parser.error(str(e))

    if args.method in ("wallace", "both"):
        try:
            tm_w = wallace_tm(args.sequence)
            print(f"Ta (Wallace)     : {tm_w:.1f} °C")
        except ValueError as e:
            print(f"Wallace error: {e}", file=sys.stderr)

    if args.method in ("santalucia", "both"):
        try:
            tm_s = santalucia_tm(args.sequence, args.conc, args.salt)
            print(f"Ta (SantaLucia)  : {tm_s:.1f} °C")
        except ValueError as e:
            print(f"SantaLucia error: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
