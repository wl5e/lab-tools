"""Command-line interface."""

import argparse
import sys
from lab_tools.molecular.restriction.core import digest_fragments
from lab_tools.molecular.restriction.enzymes import ENZYMES


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='Simulate restriction enzyme digest of DNA sequences.'
    )
    parser.add_argument(
        '-s', '--sequence',
        help='DNA sequence string (only ATCG)',
        default=None
    )
    parser.add_argument(
        '-f', '--fasta',
        help='FASTA file (first entry will be used)',
        default=None
    )
    parser.add_argument(
        '-e', '--enzymes',
        help='Comma-separated list of enzyme names (e.g., EcoRI,BamHI)'
    )
    parser.add_argument(
        '-c', '--circular',
        action='store_true',
        help='Treat sequence as circular'
    )
    parser.add_argument(
        '--list-enzymes',
        action='store_true',
        help='List all available enzymes and exit'
    )
    args = parser.parse_args(argv)

    if args.list_enzymes:
        print("Available enzymes:")
        for name, (site, tc, bc) in ENZYMES.items():
            print(f"  {name:10} {site:10} cut: {tc}/{bc}")
        return 0

    if not args.sequence and not args.fasta:
        parser.error("Either --sequence or --fasta must be provided.")
    if not args.enzymes:
        parser.error("Enzyme list is required (use --enzymes).")

    # Obtain sequence
    if args.sequence:
        seq = args.sequence.strip().upper()
        if not seq:
            parser.error("Sequence cannot be empty.")
    else:
        seq = _read_fasta(args.fasta)
        if not seq:
            parser.error("FASTA file is empty or invalid.")

    enzymes = [e.strip() for e in args.enzymes.split(',') if e.strip()]
    if not enzymes:
        parser.error("At least one enzyme required.")

    try:
        fragments = digest_fragments(seq, enzymes, circular=args.circular)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if not fragments:
        print("No fragments produced (unexpected).")
        return 0

    # Print table
    header = (f"{'#':>4} {'Start':>6} {'End':>6} {'Length':>7} "
              f"{'Left Enzyme':<12} {'Left Overhang':<14} "
              f"{'Right Enzyme':<12} {'Right Overhang':<14}")
    print(header)
    print("-" * len(header))
    for i, frag in enumerate(fragments, 1):
        print(f"{i:4} {frag['start']:6} {frag['end']:6} {frag['length']:7} "
              f"{frag['left_enzyme'] or 'None':<12} {frag['left_overhang'] or '-':<14} "
              f"{frag['right_enzyme'] or 'None':<12} {frag['right_overhang'] or '-':<14}")

    return 0


def _read_fasta(filepath: str) -> str:
    with open(filepath) as f:
        lines = f.readlines()
    if not lines or not lines[0].startswith('>'):
        raise ValueError("File does not appear to be FASTA.")
    seq_parts = []
    for line in lines[1:]:
        line = line.strip()
        if line.startswith('>'):
            break
        seq_parts.append(line.upper())
    return ''.join(seq_parts)


if __name__ == '__main__':
    sys.exit(main())
