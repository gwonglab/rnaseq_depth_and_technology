#!/usr/bin/env python3

import argparse
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-r', '--run-dir', required=True)
    ap.add_argument('-m', '--metadata', required=True)
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

    subsets = [1, 2, 3, 4, 5, 6, 8, 10]

    os.makedirs(args.output_dir, exist_ok=True)

    for subset in subsets:
        subset_runs = []
        subset_runs_complete = {}

        for run in runs:
            complete_transcripts_fn = f'{args.run_dir}/{run}/subsets/{subset}/unambiguous_complete_transcript_names.txt'
            if os.path.exists(complete_transcripts_fn):
                subset_runs.append(run)
                subset_runs_complete[run] = set()
                with open(complete_transcripts_fn) as complete_transcripts_h:
                    for line in complete_transcripts_h:
                        subset_runs_complete[run].add(line.strip())

        with open(f'{args.output_dir}/{subset}gb.txt', 'w') as output_h, \
                open(f'{args.output_dir}/{subset}gb_percent.txt', 'w') as perc_out_h:
            output_h.write("\t{}\n".format("\t".join((str(run) for run in subset_runs))))
            perc_out_h.write("\t{}\n".format("\t".join((str(run) for run in subset_runs))))

            start_idx = 0
            for run in subset_runs:
                output_h.write(f'{run}')
                perc_out_h.write(f'{run}')
                # Overlapped counts
                for idx, compared_run in enumerate(subset_runs):
                    if idx < start_idx:
                        output_h.write('\t')
                        continue
                    shared_count = len(subset_runs_complete[run] & subset_runs_complete[compared_run])
                    output_h.write(f'\t{shared_count}')
                output_h.write('\n')

                for idx, compared_run in enumerate(subset_runs):
                    shared_count = len(subset_runs_complete[run] & subset_runs_complete[compared_run])
                    fraction = shared_count / len(subset_runs_complete[run])
                    perc_out_h.write(f'\t{fraction}')
                perc_out_h.write('\n')
                start_idx += 1


if __name__ == '__main__':
    main()
