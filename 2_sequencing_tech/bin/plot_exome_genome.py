#!/usr/bin/env python3.6

import argparse
from collections import defaultdict
import sys
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
from matplotlib import lines


def millions(x, pos):
    millions = x / 1_000_000
    return f'{millions:.0f}Mb'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('-m', '--metadata', required=True)
    ap.add_argument('-r', '--run-dir', required=True)
    args = ap.parse_args()

    metadata = {}
    runs = []
    with open(args.metadata) as metadata_h:
        header = metadata_h.readline().rstrip(os.linesep).split(',')
        for line in metadata_h:
            sp_line = line.rstrip(os.linesep).split(',')
            run = sp_line[1]
            runs.append(run)
            metadata[run] = {header[idx]: ent for idx, ent in enumerate(sp_line)}

    subsets = [1, 2, 3, 4, 5, 6, 8, 10]
    run_results = defaultdict(list)

    for run in runs:
        for idx, subset in enumerate(subsets):
            subset_dir = f'{args.run_dir}/{run}/subsets/{subset}'
            if os.path.exists(subset_dir):
                if len(run_results[run]) != idx:
                    sys.stderr.write(f'Missing subset for {run}\n')
                    sys.exit(1)
                with open(f'{subset_dir}/exome_results.txt') as exome_h:
                    # Skip header
                    exome_h.readline()
                    sp_line = exome_h.readline().rstrip(os.linesep).split('\t')
                    exome_total = int(sp_line[1])
                    genome_total = int(sp_line[2])
                    run_results[run].append((exome_total, genome_total))

    pyplot.rcParams["font.family"] = "DejaVu Serif"
    pyplot.rcParams["axes.spines.top"] = False
    pyplot.rcParams["axes.spines.right"] = False
    fig, ax1 = pyplot.subplots(figsize=(6,4))

    for run, subset_results in sorted(run_results.items(), key=lambda x: x[0]):
        exome_results = [exome_result for exome_result, genome_result in subset_results]
        genome_results = [genome_result for exome_result, genome_result in subset_results]
        pyplot.plot(subsets[:len(exome_results)], exome_results, ls='--', color='#' + metadata[run]['color'])
        pyplot.plot(subsets[:len(genome_results)], genome_results, ls='-', color='#' + metadata[run]['color'])

    ax1.set_xticks(subsets)
    ax1.set_xticklabels([f'{subset}Gb' for subset in subsets])
    ax1.yaxis.set_major_formatter(FuncFormatter(millions))

    legend_handles = [lines.Line2D([], [], color='black', label='Genome'), lines.Line2D([], [], color='black', ls='--', label='Exome')]
    ax1.legend(handles=legend_handles)
    ax1.set_ylabel('Total Bases Covered by Assemblies', fontsize=12)
    ax1.set_xlabel('Subset Size', fontsize=12)
    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    pyplot.savefig(args.output)
