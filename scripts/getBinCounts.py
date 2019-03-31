#!/usr/bin/env python3

# getBinCounts.py: counts the number of reads in each of a given set
# of genomics bins
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

import sys, os, argparse
import pysam, numpy
from bisect import bisect_right


class Bin:
    def __init__(self, line):
        self.chrom, self.start, self.end = line.strip().split()[:3]
        self.start = int(self.start)
        self.end = int(self.end)
        self.count = 0
    def __str__(self):
        return '%s\t%d\t%d\t%d' % (self.chrom, self.start, self.end, self.count)
    def increment(self):
        self.count += 1


def loadChromInfo(chromInfoFile):
    chromSizes = {}
    chromNames = []
    for i in open(chromInfoFile):
        chrom, chromLen = i.strip().split()[:2]
        chromLen = int(chromLen)
        chromNames.append(chrom)
        chromSizes[chrom] = chromLen
    chromNames.sort()

    chromOffsets = {}
    totalLength = 0
    for i in chromNames:
        chromOffsets[i] = totalLength
        totalLength += chromSizes[i]
    return chromOffsets, totalLength


def loadBins(chromOffsets, binsFile):
    bins = [Bin(i) for i in open(binsFile)]
    binSorter = [(chromOffsets[bins[i].chrom] + bins[i].start, i)
                 for i in range(len(bins))]
    binSorter.sort()
    binPosAbs = [i[0] for i in binSorter]
    binLookup = [i[1] for i in binSorter]
    return bins, binPosAbs, binLookup


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile-sam', dest='inSamFile', required=True,
                        help= 'input reads in SAM format')
    parser.add_argument('-c', '--infile-chroms', dest='inChromFile',
                        required=True, help= '2 col file: chrom and chromLen')
    parser.add_argument('-b', '--infile-bins', dest='inBinsFile',
                        required=True, help= 'bins file in 3-col BED format')
    parser.add_argument('-o', '--outfile-bins', dest='outBinsFile',
                        required=True, help= 'output file for bin counts')
    parser.add_argument('-s', '--outfile-stats', dest='outStatsFile',
                        required=True, help= 'output file for stats')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print progress info to stderr')
    args = parser.parse_args()
    ########################################################################


    if args.VERBOSE:
        print('loading chrom sizes', file=sys.stdout)
    chromOffsets, genomeSize = loadChromInfo(args.inChromFile)
    if args.VERBOSE:
        print('total chroms: %d' % len(chromOffsets), file=sys.stdout)
        print('genome size: %d' % genomeSize, file=sys.stdout)


    # the "binLookup" below is used so the bins can later be output in
    # the same order they appear in the input file, which is needed
    # because other parts of the pipeline for CNV analysis are
    # designed (poorly) to depend on information associated with line
    # numbers across files, rather than relative to genomic
    # coordinates
    if args.VERBOSE:
        print('loading bins', file=sys.stdout)
    bins, binPosAbs, binLookup = loadBins(chromOffsets, args.inBinsFile)
    nBins = len(bins)
    # binChroms used to make sure only chroms that are part of the
    # bins are counted
    chromsWithBins = dict([(i.chrom, True) for i in bins])
    if args.VERBOSE:
        print('total bins: %d' % nBins, file=sys.stdout)
        print('chroms with bins: %d' % len(chromsWithBins), file=sys.stdout)


    if args.VERBOSE:
        print('counting reads in bins', file=sys.stdout)
    totalReads = 0
    mapfile = pysam.AlignmentFile(args.inSamFile, 'r', check_sq=False)
    for read in mapfile.fetch():
        if read.reference_name in chromsWithBins:
            readPosAbs = chromOffsets[read.reference_name] + read.reference_start
            binIdx = binLookup[bisect_right(binPosAbs, readPosAbs) - 1]
            bins[binIdx].increment()
            totalReads += 1
    readsPerBin = float(totalReads)/nBins
    if args.VERBOSE:
        print('total used reads: %d' % totalReads, file=sys.stdout)
        print('mean reads per bin: %d' % (readsPerBin))


    if args.VERBOSE:
        print('writing output', file=sys.stdout)
    out = open(args.outBinsFile, 'w')
    for i in range(nBins):
        print('%s\t%f' % (str(bins[i]), bins[i].count/readsPerBin), file=out)

    statsOut = open(args.outStatsFile, 'w')
    print('total reads: %d' % totalReads, file=statsOut)
    print('median bin count: %d' % numpy.median([i.count for i in bins]),
          file=statsOut)


if __name__ == "__main__":
    main()
