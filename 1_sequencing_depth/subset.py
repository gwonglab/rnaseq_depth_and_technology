#!/usr/bin/python3

import sys
import random
import os.path

random.seed(0)

# Usage: ./subset.py <subset-size> output-basename statistics-file fastq-files-1.fq 
# ./subset.py 100000 tests/subset tests/171p5-00?/171p5-00?-1.fq
# This routine is generally within 1 read-pair of the desired size.
# A better system could half this to a max difference of +/- 0.5 read-pair.



# Iterator to fetch pair of sequences, one from each of two FastQ files.
# assumes each entry is four lines and both files have the same read orders
# WARNING: This routine ASSUMES the read names are the same. It does not check.
# It will produce garbage if this constraint on the input is violated.
class readFastQ():
   def __init__(self, filename1, filename2):
      self.file1 = open(filename1, 'r') 
      self.file2 = open(filename2, 'r') 
   def __iter__(self):
      return self
   def __next__(self):
      # read next four lines from both files
      lines1 = [self.file1.__next__().rstrip() for x in range(4)]
      lines2 = [self.file2.__next__().rstrip() for x in range(4)]
      name = lines1[0][1:].rsplit('/', 1)[0]
      name = name.rsplit(' ', 1)[0]
      return (name, lines1[1], lines1[3], lines2[1], lines2[3])


# sum the base counts listed in the read statistics file
base_total = 0
for line in open(sys.argv[3]):
   base_total += int(line.split()[4])
print("Total Bases:", base_total)


# now read the FastQ files -- choosing read-pairs for output

# open the output files
output1 = open('{}-1.fq'.format(sys.argv[2]), 'w')
output2 = open('{}-2.fq'.format(sys.argv[2]), 'w')

# target are the number of bases still to be output
# remaining are the number of bases yet to be considered
target = int(sys.argv[1])
remaining = base_total

# now read the FastQ files -- choosing read-pairs for output
for file1name in sys.argv[4:]:
   # generate the second FastQ filename from the first
   file2name = file1name.replace('-1.fq', '-2.fq')
   print('Now reading:', file1name, file2name)

   for cluster in readFastQ(file1name, file2name):
      bases = len(cluster[1]) + len(cluster[3])

      # want to output target bases of the remaining bases
      if random.random() < float(target) / remaining:
         output1.write('@{}/1\n{}\n+\n{}\n'.format(cluster[0],
                                       cluster[1], cluster[2]))
         output2.write('@{}/2\n{}\n+\n{}\n'.format(cluster[0],
                                       cluster[3], cluster[4]))
         target -= bases
      remaining -= bases
output2.close()
output1.close()

# report how close to the target size the output is
print('Final output size:', int(sys.argv[1]) - target)
