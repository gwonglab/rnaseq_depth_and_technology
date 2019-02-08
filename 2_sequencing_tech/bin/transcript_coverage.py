#!/usr/bin/env python3

import sys
import os
import argparse
import attr
from collections import defaultdict


@attr.s
class Transcript(object):
    name = attr.ib()
    chromosome = attr.ib()
    exons = attr.ib(default=attr.Factory(list))
    direction = attr.ib(default="+")
    depths = attr.ib(default=attr.Factory(list))

    def add_exon(self, exon):
        if self.direction == "+":
            self.exons.append(exon)
        else:
            self.exons.insert(0, exon)

    def set_depths(self, chromosome_depths):
        for exon in self.exons:
            chromosome_range = range(exon.start, exon.end + 1) if self.direction == "+" else range(exon.end, exon.start - 1, -1)
            for chromosome_position in chromosome_range:
                if chromosome_position in chromosome_depths:
                    self.depths.append(chromosome_depths[chromosome_position])
                else:
                    self.depths.append(0)

    def print_depths(self, f=sys.stdout):
        f.write("{}\t{}\n".format(self.name, "\t".join((str(depth) for depth in self.depths))))


@attr.s
class Exon(object):
    start = attr.ib()
    end = attr.ib()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-g', '--gtf', required=True)
    ap.add_argument('-d', '--depth', required=True)
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    # Tracks wanted bases and their depths
    wanted_bases = defaultdict(dict)

    # Find bases that are contained in exons.  We only need to track them.
    transcripts = {}
    with open(args.gtf) as gtf_f:
        for line in gtf_f:
            if line[0] == "#":
                continue
            sp_line = line.rstrip(os.linesep).split("\t")
            feature_type = sp_line[2]
            if feature_type != "exon":
                continue
            chromosome = sp_line[0]
            # GTF uses 1-based numbering.  We switch to 0-based.
            start = int(sp_line[3]) - 1
            end = int(sp_line[4]) - 1
            direction = sp_line[6]
            info = {}
            for info_name, value in (info_line.split(" ") for info_line in sp_line[-1].rstrip(";").split("; ")):
                value = value[1:-1]
                if info_name != 'tag':
                    info[info_name] = value
            transcript_id = info['transcript_id']
            if transcript_id not in transcripts:
                transcripts[transcript_id] = Transcript(transcript_id, chromosome, direction=direction)
            transcripts[transcript_id].add_exon(Exon(start, end))
            for idx in range(start, end + 1):
                wanted_bases[chromosome][idx] = 0

    # Read depths
    with(open(args.depth)) as depth_f:
        for line in depth_f:
            sp_line = line.rstrip(os.linesep).split("\t")
            chromosome = sp_line[0]
            position = int(sp_line[1])
            if position in wanted_bases[chromosome]:
                depth = int(sp_line[2])
                wanted_bases[chromosome][position] = depth

    with open(args.output, "w") as f:
        for transcript in transcripts.values():
            transcript.set_depths(wanted_bases[transcript.chromosome])
            transcript.print_depths(f)
