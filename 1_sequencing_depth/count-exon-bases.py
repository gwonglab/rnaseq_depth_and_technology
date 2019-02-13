#!/usr/bin/python3

# usage ./count-exon-bases.py exons.gff align1.psl align2.psl ...

# reports how many bases are aligned in each of the psl files inside and
# outside exons listing in exons.gff
# outputs bases_in_exons bases_in_genome for each input PSL file
# bases with alignments to more than one scaffold are only counted once

import sys
import collections
import gzip

# Exon coverage in the genome is recorded by storing a list of "deltas". +1 is
# recorded at the start of an exon, -1 at the end, and 0 elsewhere.  The
# cumulative sum is then the number of exons at each position.  This collapses
# multiple exons sharing start/end positions, and is much simpler than finding
# exon overlaps directly.

# table of start and stop locations
# 3 indices - table[dataset_number][chromosome_name][position]
deltas = [collections.defaultdict(lambda: collections.defaultdict(int))]

# record exon start and stop locations from GFF/GTF file in deltas[0][...][...]
for line in open(sys.argv[1], 'r'):
   if line[0] != '#':
      field = line.rstrip().split('\t')
      if field[2] == 'exon':
         chrom = field[0]
         start = int(field[3]) - 1   # remove 1 base shift between GFF (1-based)
         end   = int(field[4]) - 1   # and Blat (zero-based) coordinates
         deltas[0][chrom][start] += 1
         deltas[0][chrom][end+1] -= 1

   # if sequence-region comment is present read the chromosome length from it
   # otherwise the max. exon coordinate is used to approximate the end position
   elif line[:18] == '##sequence-region ':
      field = line.rstrip().split(' ')
      chrom = field[1]
      end   = int(field[3]) - 1
      deltas[0][chrom][end+1] = 0  # record zero exon delta at chrom. end

# count how many bases are 1) in the genome and 2) in at least 1 exon 
exome_total = 0
genome_total = 0
for chrom in sorted(deltas[0].keys()):
   old_pos = 0
   exon_depth = 0
   exon_bases = 0  # accumlate count of bases in exons on this chromosome
   
   # exon_depth increases/decreases by 1 at the start/end of each exon
   # exon_depth is constant between old_pos and pos-1
   for pos, step in sorted(deltas[0][chrom].items()):
      # if run of bases old_pos .. pos-1 has at least one exon add to count
      if exon_depth > 0:
         exon_bases += (pos - old_pos)
      exon_depth += step
      old_pos = pos

   # add this chromosome's specific numbers to totals
   exome_total += exon_bases
   genome_total += pos  

print('Bases in {}\tGenome {}\tExome {}'.format(
                           sys.argv[1], genome_total, exome_total))


# routine to process Blat PSL table lines - record +1/-1 deltas at start/end of
# echo aligned block provided the score is within 98% of the assembly length
def recordDeltas(best, deltas):
   if int(best[17]) >= 1:     # skip unless there is at least 1 aligned block
      length = int(best[10])
      score = int(best[0])
      if (length - score) * 50 <= length:   # check score >= 98% * length
         target = best[13]
         block_sz = map(int, best[18].rstrip(',').split(','))
         t_starts = map(int, best[20].rstrip(',').split(','))
         for bl_start, bl_length in zip(t_starts, block_sz):
            deltas[target][bl_start] += 1
            deltas[target][bl_start + bl_length] -= 1

# loop over the Blat PSL files remaining in the command line 
for filename in sys.argv[2:]:
   deltas.append(collections.defaultdict(lambda: collections.defaultdict(int)))

   # read through BLAT psLayout version 3 tables and populate deltas with the
   # deltas for the aligned blocks of the best alignment for each assembly
   # this code requires that all lines for each query are adjacent in the file
   query_name = ''
   best = ['0' for i in range(21)]
   best_score = -1

   for linenumber, line in enumerate(open(filename, 'r')):
      if linenumber > 4:                        # skip four lines of header
         fields = line.split()
         if query_name != fields[9] and linenumber > 5:
            # if current line is a new query_name, process the best line 
            recordDeltas(best, deltas[-1])      # for the prior query_name
            query_name = fields[9]              # remember the current query
            best_score = -1                     # reset best score
         score = int(fields[0])
         if score > best_score:  # if line has a better score then remember it
            best = fields
            best_score = score

   # at the end-of-file process the last query_name
   recordDeltas(best, deltas[-1])

# build set with all the chromosome names from all the inputs in deltas
chrom_list = set()
for x in deltas:
   chrom_list = chrom_list.union(set(x.keys()))

# iterate over all the chromosome names
base_counts = collections.defaultdict(int)
for chrom in chrom_list:
   # build list of all delta sites (any input) for this chromosome
   sites = set()
   for x in deltas:
      sites = sites.union(set(x[chrom].keys()))

   # walk the site list integrating the depth delta lists
   depths = [0 for x in deltas]
   oldposition = 0
   
   # sites is a list of positions where aligned blocks begin/end in any input
   # so in each input all bases in the span oldposition+1 .. position have
   # the same depth of coverage - these are stored in depths
   
   for position in sorted(sites):
      # generate a 'pattern' like X--XX based on which of the inputs
      # has exonic bases over this span
      pattern = ''.join(['X' if x > 0 else '-'  for x in depths])
      
      # mantain a running total of bases seen with each unique pattern
      base_counts[pattern] += position - oldposition

      # update position and the depths by adding the deltas 
      depths = list(map(sum, zip(depths, [x[chrom][position] for x in deltas])))
      oldposition = position

# loop over the observed patterns and sum the base counts for each case 
input_count = len(sys.argv) - 2
exome_total    = [0 for x in range(input_count)]
genome_total   = [0 for x in range(input_count)]
for pattern, counts in base_counts.items():
   for i in range(input_count):
      if pattern[i+1] == 'X':
         genome_total[i] += counts
         if pattern[0] == 'X':
            exome_total[i] += counts

print('{}\t{}\t{}'.format('Filename', 'Genome (bp)', 'Exome (bp)'))
for i in range(input_count):
   print('{}\t{}\t{}'.format(sys.argv[i+2], genome_total[i], exome_total[i]))

