#!/usr/bin/env python3

import sys
import os
import argparse
from collections import defaultdict, Counter


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-g', '--gtf-file', required=True)
    ap.add_argument('-o', '--output-file', required=True)
    args = ap.parse_args()

    chrom_exon_cov = defaultdict(Counter)
    with open(args.gtf_file) as gtf_h:
        for line in gtf_h:
            if line.startswith('#'):
                continue
            sp_line = line.rstrip(os.linesep).split('\t')
            chrom = sp_line[0]
            entry_type = sp_line[2]
            if entry_type != 'exon':
                continue

            start = int(sp_line[3])
            end = int(sp_line[4])

            chrom_exon_cov[chrom].update(range(start, end + 1))

    all_transcripts = set()
    duplicate_transcripts = set()
    with open(args.gtf_file) as gtf_h:
        for line in gtf_h:
            if line.startswith('#'):
                continue
            sp_line = line.rstrip(os.linesep).split('\t')
            chrom = sp_line[0]
            entry_type = sp_line[2]
            if entry_type != 'exon':
                continue

            try:
                additional_info = {key: value[1:-1] for key, value in (entry.strip().split(' ') for entry in sp_line[8].split(';')[:-1])}
            except:
                sys.stderr.write(f'{line}\n')
                sys.exit(1)
            transcript = additional_info['transcript_id']
            start = int(sp_line[3])
            end = int(sp_line[4])
            exon_cov = (chrom_exon_cov[chrom][idx] for idx in range(start, end + 1))
            duplicate = any((chrom_exon_cov[chrom][idx] > 1 for idx in range(start, end + 1)))
            if duplicate:
                duplicate_transcripts.add(transcript)
            all_transcripts.add(transcript)

    wanted_transcripts = all_transcripts - duplicate_transcripts
    with open(args.output_file, 'w') as output_h:
        for transcript in wanted_transcripts:
            output_h.write(f'{transcript}\n')
