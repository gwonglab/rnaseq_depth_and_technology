#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-m', '--metadata', required=True)
    ap.add_argument('-r', '--run-dir', required=True)
    ap.add_argument('-o', '--output-dir', required=True)
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

    os.makedirs(args.output_dir, exist_ok=True)

    subsets = [1, 2, 3, 4, 5, 6, 8, 10]

    for subset in subsets:
        subset_runs = []
        platform_transcripts = defaultdict(set)
        platform_transcripts['BGISEQ'] = set()
        platform_transcripts['HiSeq'] = set()

        for run in runs:
            platform = metadata[run]["platform"]
            complete_transcripts_fn = f'{args.run_dir}/{run}/subsets/{subset}/unambiguous_complete_transcript_names.txt'
            if os.path.exists(complete_transcripts_fn):
                subset_runs.append(run)
                with open(complete_transcripts_fn) as complete_transcripts_h:
                    for line in complete_transcripts_h:
                        platform_transcripts[platform].add(line.strip())

        for platform, complete_transcripts in platform_transcripts.items():
            with open(f'{args.output_dir}/{platform}_{subset}_complete_transcripts.txt', 'w') as out_h:
                for complete_transcript in complete_transcripts:
                    out_h.write(f'{complete_transcript}\n')


if __name__ == '__main__':
    main()
