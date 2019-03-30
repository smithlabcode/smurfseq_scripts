#!/usr/bin/env python3

# readsFromFrags.py: generates simulated SMURF-seq reads by randomly
# combining fragments of actual reads that map uniquely
#
# Copyright (C) 2019 Rish Prabakar and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import sys, argparse
from random import shuffle
import pysam

def complement(base):
  if   base == 'A': return 'T'
  elif base == 'C': return 'G'
  elif base == 'G': return 'C'
  elif base == 'T': return 'A'
  else: return 'N'


def rev_comp(seq):
  return ''.join(list(map(complement, seq))[::-1])


class Frag:
  def __init__(self, chrom, pos, seq):
    self.chrom = chrom
    self.pos = str(pos) # we won't do arithmetic on the pos
    self.seq = seq


def main():

  parser = argparse.ArgumentParser()
  '[sam file] [frag count] [frag len] [read count] [out file]'
  parser.add_argument('-f', '--n-frags', dest='nFragsPerRead', required=True,
                      help= 'number of frags per read')
  parser.add_argument('-l', '--frag-len', dest='fragLen', type=int,
                      required=True, help= 'length of frags to simulate')
  parser.add_argument('-n', '--n-reads', dest='nReads', required=True,
                      help= 'number of reads to simulate')
  parser.add_argument('-r', '--read-len', dest='readLen', required=True,
                      help= 'length of reads to simulate')
  parser.add_argument('-s', '--sam-file', dest='samFileName',
                      required=True, help= 'mapped reads in SAM file')
  parser.add_argument('-o', '--outfile', dest='outFile',
                      required=True, help= 'output filename')
  parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                      required=False, help= 'print progress info to stderr')
  args = parser.parse_args()

  goodFrags = []

  samFile = pysam.AlignmentFile(args.samFileName, 'r', check_sq=False)
  for frag in samFile.fetch():

    # make sure the alignmnet is unique
    if (not frag.has_tag('XA') and
        not frag.has_tag('SA') and
        not frag.is_unmapped and
        frag.mapping_quality > 1):

      fragChrom, fragPos = frag.query_name.split('-')[1].split(':')
      fragPos = int(fragPos)
      if fragChrom == frag.reference_name:
        if abs(fragPos - frag.reference_start) <= 10:
          fragSeq = frag.query_alignment_sequence
          if len(fragSeq) == args.fragLen:
            if frag.is_reverse: fragSeq = revComp(fragSeq)
            goodFrags.append(Frag(frag.reference_name,
                                  frag.reference_start, fragSeq))

  shuffle(goodFrags)
  outFile = open(args.outFile, 'w')
  i = 0
  j = 0
  while j < len(seqList) and i < args.nReads:
    readSeq = ''
    readPos = []
    while len(outSeq) < args.readLen and j < len(seqList):
      readSeq += goodFrags[j].seq
      readPos.append(goodFrags[j].chrom + ':' + goodFrags[j].pos)
      j += 1
    print('>Read%d ' % i + ','.join(outPos), file=outFile)
    print(outSeq, file=outFile)
    i += 1

  out_file.close()


if __name__ == '__main__':
  main()
