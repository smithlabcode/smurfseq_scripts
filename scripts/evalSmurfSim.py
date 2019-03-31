#!/usr/bin/env python3

# evalSmurfSim.py: tabulates the summary statistics for performance of
# a SMURF-seq simulation
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
import pysam
import numpy

class Interval:
    """
    Just stores the chrom name, start and end for a genomic interval
    and only to aggregate data
    """

    def __init__(self, chrom, pos, length):
        self.chrom = chrom
        self.start = int(pos)
        self.length = int(length)
        self.end = self.start + self.length

    def matches(self, other):
        """r = self, m = other"""
        if self.chrom != other.chrom:
            return False
        else:
            overlap = 0
            if other.start < self.start:
                overlap = other.end - self.start
            elif other.end < self.end:
                overlap = other.length
            else:
                overlap = self.end - other.start
            return overlap >= self.length/2


def overlaps(a, b):
    """two list of 2 elements each: a = [x,y] and b = [x',y']"""
    return ((b[0] <= a[0] and a[0] < b[1]) or
            (a[0] <= b[0] and b[0] < a[1]))


def getStartPosFromCigar(samRead):
    """
    Extract the mapping start position from the cigar string
    correcting the 1-offset depending on strand.
    """
    if samRead.is_reverse: # reverse strand
        if (samRead.cigartuples[-1][0] == 4 or
            samRead.cigartuples[-1][0] == 5):
            return samRead.cigartuples[-1][1] - 1 # fix sam format 1-offset
    else: # forward strand
        if (samRead.cigartuples[0][0] == 4 or
            samRead.cigartuples[0][0] == 5):
            return samRead.cigartuples[0][1] - 1 # fix sam format 1-offset
    return 0


def loadPredictedQueryLocations(mappedSamFile):
    """
    Takes a SAM file name and returns a dictionary with keys
    corresponding to the reads and locations corresponding to the
    start and end positions of the mappings within that query read.
    """
    mappedFrags = {}
    mapfile = pysam.AlignmentFile(mappedSamFile, 'r', check_sq=False)
    for read in mapfile.fetch():
        if read.query_alignment_sequence and not read.is_unmapped:
            startPos = getStartPosFromCigar(read)
            if read.query_name not in mappedFrags:
                mappedFrags[read.query_name] = []
            endPos = startPos + read.query_alignment_length
            mappedFrags[read.query_name].append([startPos, endPos])
    mapfile.close()
    return mappedFrags


def loadPredictedRefLocations(mappedSamFile):
    """
    Takes a SAM file name and returns a dictionary with keys
    corresponding to read names. For each read, the value is a list of
    reference genome intervals where the predicted fragments are
    predicted to map.
    """
    mapLoc = {}
    mapfile = pysam.AlignmentFile(mappedSamFile, 'r', check_sq=False)
    for r in mapfile.fetch():
        if ((r.query_alignment_sequence) and
            (not r.is_unmapped) and (r.mapping_quality > 0)):
            if r.query_name not in mapLoc:
                mapLoc[r.query_name] = []
            mapLoc[r.query_name].append(Interval(r.reference_name,
                                                 r.reference_start,
                                                 r.query_alignment_length))
    mapfile.close()
    return mapLoc


def loadActualRefLocations(fastaFile, fragLen):
    """
    Accepts a FASTA file name of simulated SMURF-seq reads with the
    locations of fragments within each read encoded in the name line of
    each read. The locations are added to a dictionary as a list of
    genomic intervals of each fragment the read, and the keys are the
    read names.
    """
    mapLoc = {}
    for line in open(fastaFile):
        if line[0] == '>': # process only name lines
            name = line.strip().split()[0][1:]
            readPosPart = line.strip().split()[1]
            mapLoc[name] = [Interval(x.split(':')[0], x.split(':')[1], fragLen)
                            for x in readPosPart.split(',')]
    return mapLoc


def loadPredictedRefBases(mappedSamFile):
    """
    Accepts that name of a SAM format file for mapped simulated
    SMURF-seq reads. Exctracts the predicted mapping locations within
    the reference genome as a count of nucleotides covering each genomic
    position. These are contained in a two-level dictionary, with the
    first having the chrom name as key, and the second level having the
    chromosome position as key, with value equal to the number of mapped
    fragments overlapping that position in the genome. This is terribly
    inefficient, but convenient.
    """
    refBases = {}
    mapfile = pysam.AlignmentFile(mappedSamFile, 'r', check_sq=False)
    for m in mapfile.fetch():
        if m.query_alignment_sequence and not m.is_unmapped:
            if m.reference_name not in refBases:
                refBases[m.reference_name] = {}
            endPos = m.reference_start + m.reference_length
            for i in range(m.reference_start, endPos):
                if i not in refBases[m.reference_name]:
                    refBases[m.reference_name][i] = 0
                refBases[m.reference_name][i] += 1
    mapfile.close()
    return refBases


def getActualRefBases(actualRefLocations):
    """
    Accepts a dictionary with read names as keys. For each read, the
    associated value is a list of reference genome Interval. These are
    converted into a two-level dictionary with the same organization as
    the return value of the loadPredictedRefBases function. The
    difference here is that for the 'actual' ref bases we extract them
    from the ref locations rather than loading them from a file (as is
    done in the loadPredictedRefBases).
    """
    refBases = {}
    for readName in actualRefLocations:
        for r in actualRefLocations[readName]:
            if r.chrom not in refBases:
                refBases[r.chrom] = {}
            for i in range(r.start, r.end):
                if i not in refBases[r.chrom]:
                    refBases[r.chrom][i] = 0
                refBases[r.chrom][i] += 1
    return refBases


def getNucsAccuracy(actual, predicted):
    """
    tabulate the true positives, false positive and false negatives
    for identifying the mapping locations of each nucleotide in the
    simulated SMURF-seq read within the reference genome
    """
    tpNucs = 0
    fpNucs = 0
    for i in predicted:
        if i in actual:
            for j in predicted[i]:
                pVal = predicted[i][j]
                if j in actual[i]:
                    tpNucs += min(actual[i][j], pVal)
                    if pVal > actual[i][j]:
                        fpNucs += (pVal - actual[i][j])
                else: fpNucs += pVal
        else: fpNucs += sum(predicted[i].values())

    fnNucs = 0
    for i in actual:
        if i in predicted:
            for j in actual[i]:
                aVal = actual[i][j]
                if j in predicted[i]:
                    if aVal > predicted[i][j]:
                        fnNucs += aVal - predicted[i][j]
                else: fnNucs += aVal
        else: fnNucs += sum(actual[i].values())
    return tpNucs, fpNucs, fnNucs


def countQueryBasesMapped(mappedFrags):
    """
    Accepts the dictionary 'mappedFrags' which is keyed on the name of a
    read and containing a list of ordered pairs (start, end) for the
    identified fragment locations within the read. Returns the number of
    mapped bases by collapsing overlapping fragments.
    """
    mappedBases = 0
    for readName in mappedFrags:
        mapped = mappedFrags[readName]
        mapped.sort()
        nonOverl = []
        for i in mapped:
            if len(nonOverl) > 0 and overlaps(i, nonOverl[-1]):
                nonOverl[-1][1] = i[1]
            else: nonOverl.append(i)
        mappedBases += sum([(x[1] - x[0]) for x in nonOverl])
    return mappedBases


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta-file', dest='fastaFile', required=True,
                        help= 'reads in FASTA format')
    parser.add_argument('-s', '--sam-file', dest='samFileName',
                        required=True, help= 'mapped reads in SAM file')
    parser.add_argument('-l', '--frag-len', dest='fragLen', type=int,
                        required=True, help= 'actual fragment length')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print progress info to stderr')
    args = parser.parse_args()


    ########################################################################
    # predictedRefFragLengths = []
    # for i in predictedRefLocations:
    #   predictedRefFragLengths += [j.length for j in predictedRefLocations[i]]


    ########################################################################
    # load the reference genome locations of fragments for the predicted
    # mappings (from a SAM file) and for the actual mappings (from the
    # name line of the simulalted FASTA file).
    if args.VERBOSE:
        print('loading predicted and actual reference locations',
              file=sys.stderr)
    predictedRefLocations = loadPredictedRefLocations(args.samFileName)
    actualRefLocations = loadActualRefLocations(args.fastaFile, args.fragLen)


    ########################################################################
    ## process the predicted and actual reference locations for the
    ## simulated reads, and determine the true positives, false
    ## positives and flase negatives by comparing predicted with actual.
    if args.VERBOSE:
        print('computing fragment mapping accuracy summary stats',
              file=sys.stderr)
    nFalsePositives = 0
    nFalseNegatives = 0
    nTruePositives = 0
    truePositiveFragLengths = []

    for readName in actualRefLocations:
        numMatch = 0
        for r in actualRefLocations[readName]:
            matchFound = False
            minLengthDiscrepancy = sys.maxsize
            bestFragLen = 0
            nPredictedRefLocations = 0
            if readName in predictedRefLocations:
                nPredictedRefLocations = len(predictedRefLocations[readName])
                for m in predictedRefLocations[readName]:
                    if r.matches(m):
                        matchFound = True
                        currLenDiscr = abs(m.length - r.length)
                        if currLenDiscr < minLengthDiscrepancy:
                            minLengthDiscrepancy = currLenDiscr
                            bestFragLen = m.length
            if matchFound:
                truePositiveFragLengths.append(bestFragLen)
                nTruePositives += 1
                numMatch += 1
            else:
                nFalseNegatives += 1
        nFalsePositives += (nPredictedRefLocations - numMatch)


    ########################################################################
    ## summarize the statistics on precition and recall for fragment
    ## reference locations
    predictedPositiveFrag = nTruePositives + nFalsePositives
    conditionPositiveFrag = nTruePositives + nFalseNegatives
    precision = float(nTruePositives)/predictedPositiveFrag
    recall = float(nTruePositives)/conditionPositiveFrag

    meanTruPosFragLen = numpy.mean(truePositiveFragLengths)
    sdTruPosFragLen = numpy.std(truePositiveFragLengths)
    medianTruPosFragLen = numpy.median(truePositiveFragLengths)


    ########################################################################
    ## summarize statistics for individual nucleotide predictions
    if args.VERBOSE:
        print('loading ref bases for predicted mappigns', file=sys.stderr)
    predictedRefBases = loadPredictedRefBases(args.samFileName)
    if args.VERBOSE:
        print('extracting ref bases for actual mappings', file=sys.stderr)
    actualRefBases = getActualRefBases(actualRefLocations)

    if args.VERBOSE:
        print('computing accuracy of mapping per simulated base',
              file=sys.stderr)
    tpNucs, fpNucs, fnNucs = getNucsAccuracy(actualRefBases, predictedRefBases)
    predictedPositiveNucs = tpNucs + fpNucs
    conditionPositiveNucs = tpNucs + fnNucs
    precisionNucs = float(tpNucs)/predictedPositiveNucs
    recallNucs = float(tpNucs)/conditionPositiveNucs


    ########################################################################
    ## determine the fraction of nucleotides in reads that are mapped
    if args.VERBOSE:
        print('computing fraction of query bases mapped', file=sys.stderr)
    nTotalQueryBases = sum([len(i.strip())
                            for i in open(args.fastaFile) if i[0] != '>'])
    predictedQueryFrags = loadPredictedQueryLocations(args.samFileName)
    nQueryBasesMapped = countQueryBasesMapped(predictedQueryFrags)
    fracBasesMapped = float(nQueryBasesMapped)/nTotalQueryBases


    ########################################################################
    ## format and print the results
    if args.VERBOSE:
        print('writing output', file=sys.stderr)
    outputFormat = ['%.4f', # precision
                    '%.4f', # recall
                    '%.4f', # precision nuc
                    '%.4f', # recall nuc
                    '%.3f', # frac bases mapped
                    '%.3f', # mean true pos frag len
                    '%.3f', # sd of true pos frag len
                    '%d']   # median true pos frag len

    outputFormat = '\t'.join(outputFormat)
    print(outputFormat %
          (precision,
           recall,
           precisionNucs,
           recallNucs,
           fracBasesMapped,
           meanTruPosFragLen,
           sdTruPosFragLen,
           medianTruPosFragLen))

if __name__ == '__main__':
    main()
