#!/usr/bin/python2

import sys
import random

# Usage: ./subset.py <subset-size> output-basename fastq-files-1.fq 
# ./subset.py 100000 tests/subset tests/171p5-00?/171p5-00?-1.fq


# Iterator to fetch pair of sequences, one from each of two FastQ files.
class readFastQ():
   def __init__(self, filename1, filename2):
      self.file1 = open(filename1, 'r') 
      self.file2 = open(filename2, 'r') 
   def __iter__(self):
      return self
   def next(self):
      lines1 = [self.file1.next().rstrip() for x in range(4)]
      lines2 = [self.file2.next().rstrip() for x in range(4)]
      name = lines1[0][1:].rsplit('/', 1)[0]
      name = name.rsplit(' ', 1)[0]
      return (name, lines1[1], lines1[3], lines2[1], lines2[3])

# read through all FastQ files and count bases
base_total = 0
for file1name in sys.argv[3:]:
   # generate the second FastQ filename from the first
   file2name = file1name.replace('1.fastq', '2.fastq')
   print 'Now reading:', file1name, file2name
   for cluster in readFastQ(file1name, file2name):
      bases = len(cluster[1]) + len(cluster[3])
      base_total += bases

print "Total Bases:", base_total

# now repeat-- choosing read-pairs for output
random.seed(0)
output1 = open('{}/1.fq'.format(sys.argv[2]), 'w')
output2 = open('{}/2.fq'.format(sys.argv[2]), 'w')

target = int(sys.argv[1])
remaining = base_total
for file1name in sys.argv[3:]:
   # generate the second FastQ filename from the first
   file2name = file1name.replace('1.fastq', '2.fastq')
   print 'Now reading:', file1name, file2name
   for cluster in readFastQ(file1name, file2name):
      bases = len(cluster[1]) + len(cluster[3])
      # want to output target bases of the remaining bases
      if random.random() < float(target) / remaining:
         #print target, remaining
         output1.write('@{}/1\n{}\n+\n{}\n'.format(cluster[0],
                                       cluster[1], cluster[2]))
         output2.write('@{}/2\n{}\n+\n{}\n'.format(cluster[0],
                                       cluster[3], cluster[4]))
         target -= bases
      remaining -= bases
output2.close()
output1.close()

print 'Final output size:', int(sys.argv[1]) - target
