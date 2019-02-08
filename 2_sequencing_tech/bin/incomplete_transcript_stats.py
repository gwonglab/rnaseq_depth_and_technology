#!/usr/bin/env python3

import last
import sys
import os
import argparse
import array
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import matplotlib.patches as patches


def get_low_cov_seq(transcript_seq, covered_transcript_alignments, transcript_read_coverage=Counter()):
    transcript_seq = transcript_seq.upper()
    scaffold_coverage = array.array('I', [0]) * len(transcript_seq)
    for last_entry in covered_transcript_alignments:
        last.add_last_coverage(scaffold_coverage, last_entry)

    seq = ""
    for idx, (sc1, base) in enumerate(zip(scaffold_coverage, transcript_seq)):
        if sc1 == 0 and transcript_read_coverage[idx] <= 1:
            seq += base
    return seq


def get_low_cov_gaps(transcript_seq, covered_transcript_alignments, transcript_read_coverage=Counter()):
    scaffold_coverage = array.array('I', [0]) * len(transcript_seq)
    for last_entry in covered_transcript_alignments:
        last.add_last_coverage(scaffold_coverage, last_entry)

    gaps = []
    cur_region = None
    for idx, (sc1, base) in enumerate(zip(scaffold_coverage, transcript_seq)):
        covered = sc1 > 0 or transcript_read_coverage[idx] > 0
        if covered:
            if cur_region is not None:
                gaps.append(transcript_seq[cur_region[0]:cur_region[1] + 1])
            cur_region = None
        elif cur_region is None:
            cur_region = [idx, idx]
        else:
            cur_region[1] = idx
    if cur_region is not None:
        gaps.append(transcript_seq[cur_region[0]:cur_region[1] + 1])
    return gaps


def kmer_gc_frequencies(seq, kmer_size, min_size):
    seq = seq.upper()
    if kmer_size > len(seq) and len(seq) < min_size:
        return []
    current_count = Counter(seq[0:kmer_size])
    gc_freqs = [(current_count['G'] + current_count['C']) / float(sum(current_count.values()))]
    if kmer_size > len(seq):
        return gc_freqs

    for i in range(kmer_size, len(seq)):
        current_count[seq[i - kmer_size]] -= 1
        current_count[seq[i]] += 1
        gc_freqs.append((current_count['G'] + current_count['C']) / float(sum(current_count.values())))
    return gc_freqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-a', '--transcript-alignment', required=True)
    ap.add_argument('-g', '--read-genome-coverage', required=True)
    ap.add_argument('-f', '--fasta', required=True)
    ap.add_argument('-o', '--output-dir', required=True)
    ap.add_argument('--other-complete-transcripts', required=True)
    ap.add_argument('--same-complete-transcripts', required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    other_complete_transcripts = set()
    with open(args.other_complete_transcripts) as complete_h:
        for line in complete_h:
            other_complete_transcripts.add(line.strip())

    same_complete_transcripts = set()
    with open(args.same_complete_transcripts) as complete_h:
        for line in complete_h:
            same_complete_transcripts.add(line.strip())

    transcript_seqs = {}
    for title, seq in SimpleFastaParser(open(args.fasta)):
        transcript_seqs[title.split()[0]] = seq

    covered_transcript_alignments, complete_transcript_alignments, transcript_coverage, query_counts = last.analyse_transcript_alignment(args.transcript_alignment, complete_scaffold_threshold=0.0)

    # Look at transcripts that were assembled from any library on other platform, but were not assembled in any library on same platform
    missing_transcripts = other_complete_transcripts - same_complete_transcripts

    read_genome_depth = last.get_read_coverage(args.read_genome_coverage, missing_transcripts)

    with open('{}/missing_read_depth.txt'.format(args.output_dir), 'w') as missing_h:
        for transcript in missing_transcripts:
            missing_count = 0
            low_read_cov_count = 0
            transcript_scaffold_coverage = array.array('I', [0]) * len(transcript_seqs[transcript])
            transcript_length = len(transcript_scaffold_coverage)
            for last_entry in covered_transcript_alignments[transcript]:
                last.add_last_coverage(transcript_scaffold_coverage, last_entry)
            for idx, coverage in enumerate(transcript_scaffold_coverage):
                if coverage == 0:
                    missing_count += 1
                    # Is the read exon coverage low?
                    if transcript in read_genome_depth:
                        if read_genome_depth[transcript][idx] <= 1:
                            low_read_cov_count += 1
                    else:
                        low_read_cov_count += 1
            missing_h.write('{}\t{}\t{}\n'.format(transcript_length, missing_count, low_read_cov_count))

    with open('{}/k100_covered_info.txt'.format(args.output_dir), 'w') as covered_h, \
            open('{}/wholeseq_read_gap_info.txt'.format(args.output_dir), 'w') as low_cov_full_h:

        for transcript in missing_transcripts:
            covered = last.get_covered(transcript_seqs[transcript], covered_transcript_alignments[transcript])
            for region in covered:
                gc_percents = kmer_gc_frequencies(region, 100, 10)
                for gc_percent in gc_percents:
                    covered_h.write("{}\n".format(gc_percent))

            if transcript in read_genome_depth:
                read_gaps = get_low_cov_gaps(transcript_seqs[transcript], covered_transcript_alignments[transcript], read_genome_depth[transcript])
            else:
                read_gaps = get_low_cov_gaps(transcript_seqs[transcript], covered_transcript_alignments[transcript])

            for gap in read_gaps:
                if len(gap) >= 10:
                    gc_percent = last.seq_gc_content(gap)
                    low_cov_full_h.write('{}\t{}\n'.format(gc_percent, len(gap)))


if __name__ == "__main__":
    main()
