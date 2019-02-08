import attr
import os
import array
from collections import defaultdict, Counter


@attr.s
class LastEntry(object):
    score = attr.ib()
    eg2 = attr.ib()
    e = attr.ib()
    subject = attr.ib()
    subject_start = attr.ib()
    subject_align_length = attr.ib()
    subject_strand = attr.ib()
    subject_length = attr.ib()
    subject_align_string = attr.ib()
    query = attr.ib()
    query_start = attr.ib()
    query_align_length = attr.ib()
    query_strand = attr.ib()
    query_length = attr.ib()
    query_align_string = attr.ib()
    query_alignment_proportion = attr.ib()
    subject_alignment_proportion = attr.ib()
    gapless_alignment_size = attr.ib()


def last_entries(f):
    for line in f:
        if line[0] != 'a':
            continue
        sp_lines = [line.rstrip(os.linesep).split()]
        sp_lines.append(next(f).rstrip(os.linesep).split())
        sp_lines.append(next(f).rstrip(os.linesep).split())
        score = int(sp_lines[0][1].split('=')[1])
        eg2 = float(sp_lines[0][2].split('=')[1])
        e = float(sp_lines[0][3].split('=')[1])
        subject = sp_lines[1][1]
        subject_start = int(sp_lines[1][2])
        subject_align_length = int(sp_lines[1][3])
        subject_strand = sp_lines[1][4]
        subject_length = int(sp_lines[1][5])
        subject_align_string = sp_lines[1][6]
        query = sp_lines[2][1]
        query_start = int(sp_lines[2][2])
        query_align_length = int(sp_lines[2][3])
        query_strand = sp_lines[2][4]
        query_length = int(sp_lines[2][5])
        query_align_string = sp_lines[2][6]
        alignment_length = 0
        for qbase, sbase in zip(query_align_string, subject_align_string):
            if qbase != '-' and qbase != 'N' and qbase != 'n' and sbase != '-':
                alignment_length += 1
        yield(LastEntry(score, eg2, e, subject, subject_start, subject_align_length, subject_strand, subject_length, subject_align_string, query, query_start, query_align_length, query_strand, query_length, query_align_string, alignment_length / float(query_length), alignment_length / float(subject_length), alignment_length))


def best_last_entries(f):
    group = []
    prev_query = ""
    multi_count = 0
    for last_entry in last_entries(f):
        query_name = last_entry.query
        if query_name != prev_query:
            if len(group) > 0:
                if len(group) > 1:
                    multi_count += 1
                group.sort(key=lambda x: x.subject_alignment_proportion, reverse=True)
                group.sort(key=lambda x: x.query_alignment_proportion, reverse=True)
                group.sort(key=lambda x: x.score, reverse=True)
                best_score = group[0].score
                for le in filter(lambda x: x.score == best_score, group):
                    yield(le)
            group = [last_entry]
            prev_query = query_name
        else:
            group.append(last_entry)

    if len(group) > 0:
        if len(group) > 1:
            multi_count += 1
        group.sort(key=lambda x: x.subject_alignment_proportion, reverse=True)
        group.sort(key=lambda x: x.query_alignment_proportion, reverse=True)
        group.sort(key=lambda x: x.score)
        best_score = group[0].score
        for le in filter(lambda x: x.score == best_score, group):
            yield(le)


def add_last_coverage(subject_coverage, last_entry):
    idx = last_entry.subject_start
    for qbase, sbase in zip(last_entry.query_align_string, last_entry.subject_align_string):
        if qbase == '-':
            idx += 1
        elif sbase == '-':
            continue
        else:
            subject_coverage[idx] += 1
            idx += 1


def analyse_transcript_alignment(fn, complete_scaffold_threshold=0.95, complete_transcript_threshold=0.95):
    covered_transcript_alignments = defaultdict(list)
    transcript_coverage = {}
    complete_transcript_alignments = defaultdict(list)
    query_counts = defaultdict(int)
    with open(fn) as f:
        for last_entry in best_last_entries(f):
            if last_entry.query_alignment_proportion >= complete_scaffold_threshold and last_entry.subject_alignment_proportion >= complete_transcript_threshold:
                complete_transcript_alignments[last_entry.subject].append(last_entry)
                query_counts[last_entry.query] += 1
            covered_transcript_alignments[last_entry.subject].append(last_entry)
            if last_entry.subject not in transcript_coverage:
                transcript_coverage[last_entry.subject] = array.array('i', [0]) * last_entry.subject_length

            add_last_coverage(transcript_coverage[last_entry.subject], last_entry)
    return (covered_transcript_alignments, complete_transcript_alignments, transcript_coverage, query_counts)


def get_gaps(transcript_seq, covered_transcript_alignments):
    scaffold_coverage = array.array('I', [0]) * len(transcript_seq)
    for last_entry in covered_transcript_alignments:
        add_last_coverage(scaffold_coverage, last_entry)

    gaps = []
    cur_region = None
    for idx, (sc1, base) in enumerate(zip(scaffold_coverage, transcript_seq)):
        if sc1 > 0:
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


def get_covered(transcript_seq, covered_transcript_alignments):
    scaffold_coverage = array.array('I', [0]) * len(transcript_seq)
    for last_entry in covered_transcript_alignments:
        add_last_coverage(scaffold_coverage, last_entry)

    covered = []
    cur_region = None
    for idx, (sc1, base) in enumerate(zip(scaffold_coverage, transcript_seq)):
        if sc1 == 0:
            if cur_region is not None:
                covered.append(transcript_seq[cur_region[0]:cur_region[1] + 1])
            cur_region = None
        elif cur_region is None:
            cur_region = [idx, idx]
        else:
            cur_region[1] = idx
    if cur_region is not None:
        covered.append(transcript_seq[cur_region[0]:cur_region[1] + 1])
    return covered


def get_read_coverage(fn, wanted_transcripts):
    read_coverage = {}
    with open(fn) as f:
        for line in f:
            sp_line = line.rstrip(os.linesep).split("\t")
            transcript = sp_line[0]
            if transcript not in wanted_transcripts:
                continue
            read_coverage[transcript] = array.array('I')
            for idx, depth in enumerate((int(depth) for depth in sp_line[1:])):
                read_coverage[transcript].append(depth)
    return read_coverage


def seq_gc_content(seq):
    seq = seq.upper()
    gc_count = 0
    for base in seq:
        if base == 'G' or base == 'C':
            gc_count += 1
    return gc_count / float(len(seq))


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


def get_transcript_read_coverage(fn, missing_complete_transcripts):
    transcript_read_coverage = {}
    with open(fn) as f:
        for last_entry in best_last_entries(f):
            if last_entry.subject not in missing_complete_transcripts:
                continue
            elif last_entry.subject not in transcript_read_coverage:
                transcript_read_coverage[last_entry.subject] = array.array('I', [0]) * last_entry.subject_length

            add_last_coverage(transcript_read_coverage[last_entry.subject], last_entry)
    return transcript_read_coverage


def count_connectedness(transcript_coverage, transcript_read_coverage):
    # low_depth=0, covered_by_transcript=1, disconnected=2, connected=3
    connectedness = array.array('I')
    min_depth = 5

    for depth in transcript_coverage:
        connectedness.append(1 if depth > 0 else 0)

    prev = 0
    for idx, depth in enumerate(transcript_read_coverage):
        if (prev == 3 or prev == 1) and depth >= min_depth and connectedness[idx] != 1:
            connectedness[idx] = 3
        elif prev != 3 and prev != 1 and depth >= min_depth and connectedness[idx] != 1:
            connectedness[idx] = 2
        prev = connectedness[idx]
    prev = 0
    for idx, depth in reversed(list(enumerate(transcript_read_coverage))):
        if (prev == 3 or prev == 1) and depth >= min_depth and connectedness[idx] != 1:
            connectedness[idx] = 3
        prev = connectedness[idx]
    return (connectedness.count(0), connectedness.count(1), connectedness.count(2), connectedness.count(3))
