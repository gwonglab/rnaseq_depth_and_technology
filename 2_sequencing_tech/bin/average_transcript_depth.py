#!/usr/bin/env python3

import os
import argparse

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-t', '--transcript-names-file')
    ap.add_argument('-o', '--output-file', required=True)
    ap.add_argument('-c', '--coverage-file', required=True)
    args = ap.parse_args()

    transcript_names = None
    if args.transcript_names_file is not None:
        transcript_names = set()
        with open(args.transcript_names_file) as transcript_names_h:
            transcript_names.update((line.strip() for line in transcript_names_h))

    with open(args.coverage_file) as coverage_h, \
            open(args.output_file, 'w') as output_h:
        for line in coverage_h:
            sp_line = line.rstrip(os.linesep).split('\t')
            transcript = sp_line[0]
            base_coverage = [int(entry) for entry in sp_line[1:]]
            average_depth = sum(base_coverage) / len(base_coverage)
            if transcript_names is None or transcript in transcript_names:
                output_h.write(f'{transcript}\t{average_depth}\n')
