#!/usr/bin/env python3

import sys
import os
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter, defaultdict
from operator import itemgetter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


def kmer_gc_proportions(seq, depths, avg_depth):
    if len(seq) < 100:
        return
    kmer_gc_proportions.counter += 1
    seq = seq.upper()
    base_count = Counter(seq[:100])
    avg_kmer_depth = sum(depths[:100]) / 100
    yield(((base_count['G'] + base_count['C']) / 100, avg_kmer_depth / avg_depth))
    for idx, base in enumerate(seq[100:], 100):
        base_count[seq[idx - 100]] -= 1
        base_count[base] += 1
        avg_kmer_depth = sum(depths[idx - 100:idx + 1]) / 100
        yield(((base_count['G'] + base_count['C']) / 100, avg_kmer_depth / avg_depth))


kmer_gc_proportions.counter = 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--depths-file', required=True)
    ap.add_argument('-t', '--transcript-names-file', required=True)
    ap.add_argument('-f', '--transcript-fasta-file', required=True)
    ap.add_argument('--title', required=True)
    ap.add_argument('-l', '--library', required=True)
    ap.add_argument('-m', '--metadata', required=True)
    ap.add_argument('-o', '--output-file', required=True)
    args = ap.parse_args()

    metadata = {}
    with open(args.metadata) as metadata_h:
        header = metadata_h.readline().rstrip(os.linesep).split(',')
        for line in metadata_h:
            sp_line = line.rstrip(os.linesep).split(',')
            library = sp_line[1]
            metadata[library] = {header[idx]: ent for idx, ent in enumerate(sp_line)}

    color = metadata[args.library]['color']
    wanted_transcripts = set()

    with open(args.transcript_names_file) as transcript_names_h:
        wanted_transcripts.update((line.strip() for line in transcript_names_h))

    transcript_seqs = {}
    with open(args.transcript_fasta_file) as transcript_fasta_h:
        for title, seq in SimpleFastaParser(transcript_fasta_h):
            transcript_name = title.split(' ')[0]
            if transcript_name in wanted_transcripts:
                transcript_seqs[transcript_name] = seq

    transcript_depths = {}
    with open(args.depths_file) as depths_h:
        for line in depths_h:
            sp_line = line.rstrip(os.linesep).split('\t')
            transcript = sp_line[0]
            if transcript in wanted_transcripts:
                transcript_depths[transcript] = [int(entry) for entry in sp_line[1:]]

    all_kmer_gc_proportions = []
    for transcript in wanted_transcripts:
        depths = transcript_depths[transcript]
        seq = transcript_seqs[transcript]
        avg_depth = sum(depths) / len(depths)
        all_kmer_gc_proportions.extend(kmer_gc_proportions(seq, depths, avg_depth))

    all_kmer_gc_proportions.sort(key=itemgetter(0))

    kmer_gc_proportions_hsh = defaultdict(list)
    for gc_prop, depth_prop in all_kmer_gc_proportions:
        kmer_gc_proportions_hsh[gc_prop].append(depth_prop)
    for gc_prop, depth_prop_list in kmer_gc_proportions_hsh.items():
        kmer_gc_proportions_hsh[gc_prop] = sum(depth_prop_list) / len(depth_prop_list)

    fig, ax1 = pyplot.subplots()
    ax1.set_ylim((0, 1.5))

    x_vals = [x for x, y in kmer_gc_proportions_hsh.items()]
    y_vals = [y for x, y in kmer_gc_proportions_hsh.items()]

    ax1.scatter(x_vals, y_vals, color='#{}'.format(color))
    ax1.set_xlabel('GC Proportion (100-base windows)')
    ax1.set_ylabel('Relative coverage')
    ax1.set_title(args.title)

    pyplot.savefig(args.output_file, dpi=600)
    sys.stderr.write(f'transcripts: {kmer_gc_proportions.counter}\n')
