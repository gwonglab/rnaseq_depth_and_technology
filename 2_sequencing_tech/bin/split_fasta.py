#!/usr/bin/env python3

import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fasta_in')
    ap.add_argument('splits', type=int)
    ap.add_argument('fasta_out_prefix')
    args = ap.parse_args()

    digits = len(str(args.splits))
    output_files = []
    output_lengths = [0] * args.splits
    for i in range(0, args.splits):
        padded_num = str(i).zfill(digits)
        output_files.append(open('{0}_part_{1}.fa'.format(args.fasta_out_prefix, padded_num), 'w'))

    with open(args.fasta_in) as fasta_h:
        for title, seq in SimpleFastaParser(fasta_h):
            min_length = min(output_lengths)
            min_idx = output_lengths.index(min_length)
            output_files[min_idx].write(">{}\n{}\n".format(title, seq))
            output_lengths[min_idx] += len(seq)

    for f in output_files:
        f.close()


if __name__ == '__main__':
    main()
