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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--n-frags', dest='nFragsPerRead', type=int,
                        required=True, help= 'number of frags per read')
    parser.add_argument('-l', '--frag-len', dest='fragLen', type=int,
                        required=True, help= 'length of frags to simulate')
    parser.add_argument('-r', '--n-reads', dest='nReads', type=int,
                        required=True, help= 'number of reads to simulate')
    parser.add_argument('-b', '--bed-file', dest='fragBedFileName',
                        required=True, help= 'good frags in bed file')
    parser.add_argument('-n', '--name', dest='dataSet',
                        required=True, help= 'name of data set')
    parser.add_argument('-o', '--outfile', dest='outFile',
                        required=True, help= 'output filename')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print progress info to stderr')
    args = parser.parse_args()

    if args.VERBOSE:
        print('loading frags', file=sys.stderr)
    frags = []
    for frag in open(args.fragBedFileName):
        chrom, pos, dummy, seq = frag.split()[:4]
        frags.append((chrom, pos, seq))
    if args.VERBOSE:
        print('n frags: %d' % len(frags), file=sys.stderr)

    if args.VERBOSE:
        print('shuffling frags', file=sys.stderr)
    shuffle(frags)

    if args.VERBOSE:
        print('generating and outputting reads', file=sys.stderr)
    outFile = open(args.outFile, 'w')
    readLen = args.nFragsPerRead*args.fragLen
    i = 0
    j = 0
    while j < len(frags) and i < args.nReads:
        readSeq = ''
        readPos = []
        while len(readSeq) < readLen and j < len(frags):
            readSeq += frags[j][2]
            readPos.append(frags[j][0] + ':' + frags[j][1])
            j += 1
        print('>%s_%d_%d ' % (args.dataSet, args.fragLen, i),
              ','.join(readPos), file=outFile)
        print(readSeq, file=outFile)
        i += 1
    outFile.close()

if __name__ == '__main__':
    main()
