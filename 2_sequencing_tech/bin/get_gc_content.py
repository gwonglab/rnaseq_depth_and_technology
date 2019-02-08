#!/usr/bin/env python3

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter


def max_kmer_gc(seq, kmer_size):
    current_count = Counter(seq[0:kmer_size])
    max_gc_freq = (current_count['G'] + current_count['C']) / float(sum(current_count.values()))
    for i in range(kmer_size, len(seq)):
        current_count[seq[i - kmer_size]] -= 1
        current_count[seq[i]] += 1
        current_gc_freq = (current_count['G'] + current_count['C']) / float(sum(current_count.values()))
        max_gc_freq = max(current_gc_freq, max_gc_freq)
    return max_gc_freq


def seq_gc(seq):
        seq = seq.upper()
        count = Counter(seq)
        return (count['G'] + count['C']) / float(sum(count.values()))


with open(sys.argv[1]) as fasta_f:
    for title, seq in SimpleFastaParser(fasta_f):
        sp_title = title.split()
        sys.stdout.write("{}\t{}\n".format(sp_title[0], max_kmer_gc(seq, 100)))
