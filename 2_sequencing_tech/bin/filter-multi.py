#!/usr/bin/python2

# usage ./filter-multi.py exons.gff align1.psl align2.psl ... special.psl  

# reports how many bases are aligned in each of the psl files inside and outside exons
# outputs four columns for each input PSL file:
# bases_in_exons
# bases_in_genome
# bases_in_exon_shared_with_special
# bases_in_genome_shared_with_special


import sys
import collections
import gzip

# populate tables of start and stop locations
deltas = [collections.defaultdict(lambda: collections.defaultdict(int))]

# populate tables of exon start and stop locations from GFF/GTF file
for line in open(sys.argv[1], 'r'):
   if line[0] != '#':
      field = line.rstrip().split('\t')
      if field[2] == 'exon':
         chrom = field[0].split('.', 1)[0]
         start = int(field[3]) - 1   # correct 1 base shift between GFF and blat file
         end   = int(field[4]) - 1
         deltas[0][chrom][start] += 1
         deltas[0][chrom][end+1] -= 1



# record a +1 at the start/-1 at the end of each aligned block
# provided the score is within 98% of the assembly length
def recordDeltas(best, deltas):
   length = int(best[10])
   score = int(best[0])
   if (length - score) * 50 <= length:
      target = best[13]
      block_sz = map(int, best[18].rstrip(',').split(','))
      t_starts = map(int, best[20].rstrip(',').split(','))
      for bl_start, bl_length in zip(t_starts, block_sz):
         deltas[target][bl_start] += 1
         deltas[target][bl_start + bl_length] -= 1

for filename in sys.argv[2:]:
   deltas.append(collections.defaultdict(lambda: collections.defaultdict(int)))
   # read through a BLAT *.PSL table and populate deltas with deltas for the
   # the aligned blocks of the best alignment for each assembly
   assembly = ''
   best = ['0' for i in range(21)]
   best_score = 0
   #deltas = collections.defaultdict(lambda: collections.defaultdict(int))
   for linenumber, line in enumerate(file(filename, 'r')):
         fields = line.split()
         if assembly != fields[9] and linenumber > 0:
            recordDeltas(best, deltas[-1])
            assembly = fields[9]
            best_score = -1
         score = int(fields[0])
         if score > best_score:
            best = fields
            best_score = score

   recordDeltas(best, deltas[-1])

# build list of all the chromosomes - use dict to handle repeats
chrom_list = collections.defaultdict(int)
for d in deltas:
   for chrom in d.iterkeys():
      chrom_list[chrom] +=1

# iterate over all the chromosome names
base_counts = collections.defaultdict(int)
for chrom in chrom_list.iterkeys():
   # build list of all sites for this chromosome - dict to handle repeats
   sites = collections.defaultdict(int)
   for x in deltas:
      for site in x[chrom]:
         sites[site] += 1

   # walk the site list integrating the depth delta lists
   depths = [0 for x in deltas]
   oldposition = 0
   for position in sorted(sites.keys()):

      pattern = ['X' if x > 0 else '-'  for x in depths]
      base_counts[ ''.join(pattern) ] += position - oldposition
      depths = map(sum, zip(depths, [x[chrom][position] for x in deltas]))

      oldposition = position

input_count = len(sys.argv) - 2
exon_total     = [0 for x in range(input_count)]
genome_total   = [0 for x in range(input_count)]
exon_overlap   = [0 for x in range(input_count)]
genome_overlap = [0 for x in range(input_count)]
for pattern in base_counts.iterkeys():
   for i in range(input_count):
      if pattern[i+1] == 'X':
         genome_total[i] += base_counts[pattern]
         if pattern[-1] == 'X':
            genome_overlap[i] += base_counts[pattern]
         if pattern[0] == 'X':
            exon_total[i] += base_counts[pattern]
            if pattern[-1] == 'X':
               exon_overlap[i] += base_counts[pattern]

print "File\tExon Total\tGenome Total\tExon Overlap\tGenome Overlap"
for i in range(input_count):
   print sys.argv[2+i], '\t', exon_total[i], '\t', genome_total[i], '\t', exon_overlap[i], '\t', genome_overlap[i], '\t'

## output the base_counts
#output = {}
#for pattern in sorted(base_counts.iterkeys()):
#   if pattern not in output:
#      pat1 = '-' + pattern[1:]
#      pat2 = 'X' + pattern[1:]
#      print pattern[1:], '\t', base_counts[pat2], '\t', base_counts[pat1] + base_counts[pat2]
#      output[pat1] = 0 
#      output[pat2] = 0
