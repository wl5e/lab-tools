#!/usr/bin/env python3
"""qPCR ΔΔCq Analyzer.

Core logic for relative gene expression analysis (ΔΔCq method) with error
propagation, plus a command-line interface.
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from typing import Dict, List, NamedTuple, Optional


class Result(NamedTuple):
    """One row of output data."""
    sample: str
    gene: str
    mean_cq_gene: float
    sd_cq_gene: float
    mean_cq_ref: float
    sd_cq_ref: float
    dCq: float
    sd_dCq: float
    ddCq: float
    sd_ddCq: float
    fold_change: float
    fold_change_low: float
    fold_change_high: float


def _read_cq_data(csv_path: str, delimiter: str) -> Dict[str, Dict[str, List[float]]]:
    """Parse CSV and return raw Cq values organised as data[sample][gene] = list of floats."""
    data: Dict[str, Dict[str, List[float]]] = defaultdict(lambda: defaultdict(list))

    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError('CSV file appears to be empty or lacks a header row.')
        fieldnames = [name.strip() for name in reader.fieldnames]
        missing = {'Sample', 'Gene', 'Cq'} - set(fieldnames)
        if missing:
            raise ValueError(f'Missing required columns: {",".join(missing)}. Found columns: {",".join(fieldnames)}')

        for row_num, row in enumerate(reader, start=2):
            sample = row.get('Sample', '').strip()
            gene = row.get('Gene', '').strip()
            cq_str = row.get('Cq', '').strip()
            if not sample or not gene:
                print(f'Warning: row {row_num} has empty Sample or Gene, skipping.', file=sys.stderr)
                continue
            try:
                cq_val = float(cq_str)
            except (ValueError, TypeError):
                print(f'Warning: row {row_num} Cq value "{cq_str}" is not numeric, skipping.', file=sys.stderr)
                continue
            data[sample][gene].append(cq_val)

    if not data:
        raise ValueError('No valid (Sample, Gene, Cq) rows found in the input file.')
    return data


def _compute_means_and_sds(data: Dict[str, Dict[str, List[float]]]) -> Dict[str, Dict[str, tuple]]:
    """For each (sample, gene) pair compute (mean, sample_sd)."""
    summary: Dict[str, Dict[str, tuple]] = defaultdict(dict)
    for sample, genes in data.items():
        for gene, values in genes.items():
            mean = statistics.mean(values)
            sd = statistics.stdev(values) if len(values) > 1 else 0.0
            summary[sample][gene] = (mean, sd)
    return summary


def _find_control_sample(samples: List[str], control_id: str) -> str:
    """Return the sample name that contains the control_id substring (case-insensitive).
    Raise if none or multiple match."""
    matches = [s for s in samples if control_id.lower() in s.lower()]
    if not matches:
        raise ValueError(f'No sample matches control pattern "{control_id}". Available samples: {", ".join(sorted(samples))}')
    if len(matches) > 1:
        raise ValueError(f'Multiple samples match control pattern "{control_id}": {", ".join(matches)}. Provide a more specific substring.')
    return matches[0]


def analyze_qpcr(
    csv_path: str,
    ref_gene: str,
    control_id: str,
    delimiter: str = ','
) -> List[Result]:
    """Run ΔΔCq analysis and return list of Result rows.

    Parameters
    ----------
    csv_path : path to CSV with columns Sample, Gene, Cq.
    ref_gene : name of the reference (housekeeping) gene.
    control_id : substring that uniquely identifies the control sample.
    delimiter : field delimiter in CSV.

    Returns
    -------
    list of Result namedtuples.
    """
    # 1. Parse
    raw = _read_cq_data(csv_path, delimiter)

    # 2. Summary statistics
    means = _compute_means_and_sds(raw)

    # 3. Identify control sample
    all_samples = sorted(means.keys())
    control_sample = _find_control_sample(all_samples, control_id)

    if ref_gene not in means.get(control_sample, {}):
        available = list(means[control_sample].keys()) if control_sample in means else []
        raise ValueError(f'Reference gene "{ref_gene}" not found in control sample "{control_sample}". Genes available: {available}')

    # 4. Compute ΔCq for every sample/gene (excluding ref_gene)
    dCq_data: Dict[str, Dict[str, tuple]] = {}  # sample -> gene -> (dCq, sd_dCq)
    ref_mean_control, ref_sd_control = means[control_sample][ref_gene]

    for sample, genes in means.items():
        if sample in dCq_data:
            continue
        dCq_data[sample] = {}
        if ref_gene not in genes:
            print(f'Warning: reference gene "{ref_gene}" missing in sample "{sample}", skipping all genes for this sample.', file=sys.stderr)
            continue
        ref_mean, ref_sd = genes[ref_gene]
        for gene, (mean_g, sd_g) in genes.items():
            if gene == ref_gene:
                continue
            dCq = mean_g - ref_mean
            # Error propagation: SD_dCq = sqrt(SD_g^2 + SD_ref^2)
            sd_dCq = math.sqrt(sd_g**2 + ref_sd**2)
            dCq_data[sample][gene] = (dCq, sd_dCq)

    # 5. Control ΔCq values for each gene
    if control_sample not in dCq_data:
        raise ValueError(f'Control sample "{control_sample}" has no analysable genes (may lack reference gene or GOIs).')
    control_dCq = dCq_data[control_sample]

    # 6. Calculate ΔΔCq and fold change for all non-control samples
    results = []
    for sample, genes in dCq_data.items():
        if sample == control_sample:
            continue
        for gene, (dCq, sd_dCq) in genes.items():
            if gene not in control_dCq:
                print(f'Warning: gene "{gene}" not found in control sample, skipping for sample "{sample}".', file=sys.stderr)
                continue
            dCq_ctl, sd_ctl = control_dCq[gene]
            ddCq = dCq - dCq_ctl
            # Error propagation: SD_ddCq = sqrt(SD_dCq^2 + SD_dCq_control^2)
            sd_ddCq = math.sqrt(sd_dCq**2 + sd_ctl**2)

            # Fold change 2^(-ddCq)
            fc = math.pow(2, -ddCq)
            fc_high = math.pow(2, -(ddCq - sd_ddCq))
            fc_low = math.pow(2, -(ddCq + sd_ddCq))

            # Retrieve raw mean/SD of the gene and ref for completeness
            mean_gene, sd_gene = means[sample][gene]
            mean_ref, sd_ref = means[sample][ref_gene]

            results.append(Result(
                sample=sample,
                gene=gene,
                mean_cq_gene=round(mean_gene, 4),
                sd_cq_gene=round(sd_gene, 4),
                mean_cq_ref=round(mean_ref, 4),
                sd_cq_ref=round(sd_ref, 4),
                dCq=round(dCq, 4),
                sd_dCq=round(sd_dCq, 4),
                ddCq=round(ddCq, 4),
                sd_ddCq=round(sd_ddCq, 4),
                fold_change=round(fc, 4),
                fold_change_low=round(fc_low, 4),
                fold_change_high=round(fc_high, 4),
            ))

    return results


def export_csv(results: List[Result], output, delimiter: str = ',') -> None:
    """Write results as CSV to a file path or a file-like object."""
    fieldnames = [
        'Sample', 'Gene',
        'mean_Cq_gene', 'SD_Cq_gene',
        'mean_Cq_ref', 'SD_Cq_ref',
        'dCq', 'SD_dCq',
        'ddCq', 'SD_ddCq',
        'fold_change', 'fold_change_low', 'fold_change_high'
    ]
    if isinstance(output, str):
        f = open(output, 'w', newline='')
        own_handle = True
    else:
        f = output
        own_handle = False

    try:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for r in results:
            writer.writerow({
                'Sample': r.sample,
                'Gene': r.gene,
                'mean_Cq_gene': r.mean_cq_gene,
                'SD_Cq_gene': r.sd_cq_gene,
                'mean_Cq_ref': r.mean_cq_ref,
                'SD_Cq_ref': r.sd_cq_ref,
                'dCq': r.dCq,
                'SD_dCq': r.sd_dCq,
                'ddCq': r.ddCq,
                'SD_ddCq': r.sd_ddCq,
                'fold_change': r.fold_change,
                'fold_change_low': r.fold_change_low,
                'fold_change_high': r.fold_change_high,
            })
    finally:
        if own_handle:
            f.close()


def main(argv: Optional[List[str]] = None) -> int:
    """Command-line interface for qPCR ΔΔCq analysis."""
    parser = argparse.ArgumentParser(
        description='Relative gene expression analysis (ΔΔCq method) with error propagation.'
    )
    parser.add_argument('-i', '--input', required=True, help='Input CSV file with columns: Sample,Gene,Cq')
    parser.add_argument('-o', '--output', default=None, help='Output CSV file (default: stdout)')
    parser.add_argument('--ref-gene', default='GAPDH', help='Reference gene name (default: GAPDH)')
    parser.add_argument('--control', default='control',
                        help='Substring that uniquely identifies the control sample (case-insensitive, default: control)')
    parser.add_argument('--sep', default=',', help='CSV delimiter (default: comma)')
    args = parser.parse_args(argv)

    try:
        results = analyze_qpcr(
            csv_path=args.input,
            ref_gene=args.ref_gene,
            control_id=args.control,
            delimiter=args.sep
        )
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        return 1

    if not results:
        print('No results generated.', file=sys.stderr)
        return 0

    if args.output:
        export_csv(results, args.output, delimiter=args.sep)
        print(f'Results written to {args.output}', file=sys.stderr)
    else:
        export_csv(results, sys.stdout, delimiter='\t')

    return 0


if __name__ == '__main__':
    sys.exit(main())
