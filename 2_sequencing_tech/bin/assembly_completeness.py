#!/usr/bin/env python3

import last
import argparse
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import matplotlib.patches as patches


def get_gap_ranges(depths):
    ranges = []
    in_gap = False
    for idx, depth in enumerate(depths):
        if depth == 0:
            if not in_gap:
                ranges.append([idx, idx])
                in_gap = True
            else:
                ranges[-1][1] = idx
        else:
            in_gap = False
    return ranges


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-a', '--alignment', required=True)
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('-c', '--complete-transcripts', required=True)
    args = ap.parse_args()

    covered_transcript_alignments, complete_transcript_alignments, transcript_coverage, query_counts = last.analyse_transcript_alignment(args.alignment, complete_scaffold_threshold=0.0)

    full_transcript_count = len(complete_transcript_alignments.keys())
    unambiguous_transcript_alignments = defaultdict(list)

    for alignments in complete_transcript_alignments.values():
        for alignment in alignments:
            if query_counts[alignment.query] == 1:
                unambiguous_transcript_alignments[alignment.subject].append(alignment)
    unambiguous_transcript_count = len(unambiguous_transcript_alignments.keys())

    unambiguous_transcripts = set(unambiguous_transcript_alignments.keys())
    complete_transcripts = set(complete_transcript_alignments.keys())
    with open(args.output, 'w') as output_h, open(args.complete_transcripts, 'w') as complete_h:
        output_h.write("{}\t{}\n".format(full_transcript_count, unambiguous_transcript_count))
        for transcript in unambiguous_transcripts:
            complete_h.write("{}\n".format(transcript))
