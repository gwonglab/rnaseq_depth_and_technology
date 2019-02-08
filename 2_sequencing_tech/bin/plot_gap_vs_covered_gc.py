#!/usr/bin/env python3

import sys
import os
import argparse
import numpy
from itertools import cycle
from collections import Counter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib import lines


def cumulative_values(arr):
    new_arr = []
    cur_sum = 0
    for val in arr:
        cur_sum += val
        new_arr.append(cur_sum)
    return new_arr


def cumulant_x_y(arr):
    gc_counts = Counter(arr)
    total = sum(gc_counts.values())
    x_vals = [x for x in sorted(gc_counts.keys())]
    y_vals = []
    cur_cum = 0
    for x in x_vals:
        cur_cum += gc_counts[x]
        y_vals.append(cur_cum / total)
    return (x_vals, y_vals)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-t', '--title', default="")
    p.add_argument('-o', '--out', required=True)
    p.add_argument('files', nargs='+')
    p.add_argument('-m', '--metadata', required=True)
    args = p.parse_args()
    if args.title.startswith('BGISEQ'):
        args.title = args.title.replace('BGISEQ', 'DNBSeq$^\mathrm{TM}$', 1)

    args.files.sort()

    suffixes = ['covered_info.txt', 'gap_info.txt']
    for f, suffix in zip(args.files, cycle(suffixes)):
        if not f.endswith(suffix):
            sys.stderr.write(f'Error: {f} should end with {suffix}\n')
            sys.exit(1)
    if len(args.files) % 2 != 0:
        sys.stderr.write('Wrong number of files.\n')
        sys.exit(1)

    metadata = {}
    with open(args.metadata) as metadata_h:
        header = metadata_h.readline().rstrip(os.linesep).split(',')
        for line in metadata_h:
            sp_line = line.rstrip(os.linesep).split(',')
            library = sp_line[1]
            metadata[library] = {header[idx]: ent for idx, ent in enumerate(sp_line)}

    fig, ax1 = pyplot.subplots()

    legend_handles = []
    assembled_leg = lines.Line2D([], [], color='black', label='Assembled')
    gap_leg = lines.Line2D([], [], color='black', ls='--', label='Gap')
    legend_handles.append(assembled_leg)
    legend_handles.append(gap_leg)

    for covered_file, gap_file in zip(args.files[::2], args.files[1::2]):
        gc_gap_proportions = numpy.loadtxt(gap_file, dtype=float, usecols=(0))
        gc_gap_proportions.sort()
        gc_gap_proportions *= 100
        gc_gap_proportions = numpy.around(gc_gap_proportions)

        gc_covered_proportions = numpy.fromfile(covered_file, dtype=float, sep=' ')
        gc_covered_proportions.sort()
        gc_covered_proportions *= 100
        gc_covered_proportions = numpy.around(gc_covered_proportions)

        library = os.path.dirname(covered_file).split('/')[1]

        gap_x, gap_y = cumulant_x_y(gc_gap_proportions)
        cov_x, cov_y = cumulant_x_y(gc_covered_proportions)
        pyplot.plot(gap_x, gap_y, color='#' + metadata[library]['color'], linestyle='--')
        pyplot.plot(cov_x, cov_y, color='#' + metadata[library]['color'], linestyle='-')

    pyplot.title(args.title)
    pyplot.ylabel("Transcript Proportion")
    pyplot.xlabel("GC %")
    pyplot.ylim(0, 1.1)
    pyplot.xlim(0, 100)

    first_legend = pyplot.legend(handles=legend_handles, loc=4)
    pyplot.gca().add_artist(first_legend)

    pyplot.savefig(args.out, dpi=600)


if __name__ == '__main__':
    main()
