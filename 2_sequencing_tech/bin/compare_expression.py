#!/usr/bin/env python3

import sys
import os
import argparse
from collections import defaultdict
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('transcript_gc')
    ap.add_argument('illumina_quant')
    ap.add_argument('zebra_quant')
    ap.add_argument('-c', '--tpm-column', default=5, type=int)
    ap.add_argument('-t', '--title', default="")
    ap.add_argument('-o', '--output-file', default='results.png')
    ap.add_argument('-x1', '--xstart', type=float)
    ap.add_argument('-x2', '--xend', type=float)
    ap.add_argument('-fx')
    ap.add_argument('-fy')
    args = ap.parse_args()
    args.tpm_column -= 1

    transcript_gc_proportion = {}
    with open(args.transcript_gc) as f:
        for line in f:
            sp_line = line.rstrip(os.linesep).split('\t')
            transcript = sp_line[0]
            gc_proportion = float(sp_line[1])
            transcript_gc_proportion[transcript] = gc_proportion

    transcript_tpms = defaultdict(lambda: {'zebra': 0, 'illumina': 0})

    with open(args.illumina_quant) as f:
        f.readline()
        ill_name = args.illumina_quant.split('/')[1]
        for line in f:
            sp_line = line.rstrip(os.linesep).split('\t')
            name = sp_line[0]
            transcript_tpms[name]['illumina'] = float(sp_line[args.tpm_column])

    with open(args.zebra_quant) as f:
        f.readline()
        zeb_name = args.zebra_quant.split('/')[1]
        for line in f:
            sp_line = line.rstrip(os.linesep).split('\t')
            name = sp_line[0]
            transcript_tpms[name]['zebra'] = float(sp_line[args.tpm_column])

    x = []
    y = []
    for transcript, tpms in transcript_tpms.items():
        illumina_tpm = tpms['illumina']
        zebra_tpm = tpms['zebra']
        if zebra_tpm > 10 and illumina_tpm > 10:
            x.append(transcript_gc_proportion[transcript])
            y.append(illumina_tpm / zebra_tpm)

    reorder = sorted(range(len(x)), key=lambda i: x[i])
    x = [x[i] for i in reorder]
    y = numpy.log([y[i] for i in reorder])
    fit_line = numpy.polyfit(x, y, 1)

    p = numpy.poly1d(fit_line)
    xl = [0.0, 1.0]
    yl = [p(xx) for xx in xl]

    fig, ax = pyplot.subplots()

    ax.plot(x, y, 'r.')
    ax.plot(xl, yl, 'b-')
    # sys.stdout.write("{}\n".format(args.title))
    ax.set_title(bytes(args.title, 'utf-8').decode("unicode_escape"))
    ax.set_ylabel("log({} expression/{} expression)".format(ill_name, zeb_name))
    ax.set_xlabel("Transcript GC Proportion")
    ax.set_ylim(-5, 5)
    if args.xstart is not None and args.xend is not None:
        ax.set_xlim(args.xstart, args.xend)
    ax.text(0.31, -3.0, "{:.2f}x + {:.2f}".format(*fit_line))
    sys.stderr.write("Writing to {}...".format(args.output_file))
    pyplot.savefig(args.output_file, dpi=600)
