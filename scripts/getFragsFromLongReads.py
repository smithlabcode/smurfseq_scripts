#!/usr/bin/env python3

# getFragsFromLongReads.py: generates short candidate fragments in a
# 6-col BED format, with the sequene in the 4th column, from mapped
# long reads provided in SAM format.
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
from random import randint
import pysam

def complement(base):
    """complement nucleotide or return N if not a valid nuc"""
    if   base == 'A': return 'T'
    elif base == 'C': return 'G'
    elif base == 'G': return 'C'
    elif base == 'T': return 'A'
    else: return 'N'

def revcomp(seq):
    """get the reverse complement of a DNA sequence"""
    return ''.join(list(map(complement, seq))[::-1])


## functions below use the cigar string as the argument "a" and the
## [0] element indictes m+mm/ins/del as 0/1/2. The [1] element
## indicates the length of the insertion or deletion (or 0 if m/mm).
def notDeletion(a):
    return a[0] == 0 or a[0] == 1

def isInsertion(a):
    return a[0] == 1

def operationSize(a):
    return a[1]

def findFragEndInReference(read, fragStart, fragEnd,
                           refStart, cigarIdx, cigarOffset):

    refEnd = refStart
    fragPos = fragStart # keeps current position in the fragment

    # take care of case where we are inside the match or insertion
    # from the previous fragment; after this step, the cigarOffset
    # will be 0 unless the remaining part of the current
    # match/insertion is longer than the size of a fragment
    if cigarIdx < len(read.cigartuples) and cigarOffset > 0:
        remaining = operationSize(read.cigartuples[cigarIdx]) - cigarOffset
        fragSize = fragEnd - fragStart
        if remaining <= fragSize:
            refEnd += remaining
            fragPos += remaining
            cigarIdx += 1
            cigarOffset = 0
        else:
            refEnd += fragSize
            fragPos += fragSize
            cigarOffset += fragSize

    # Now iterate along the cigar operations, adding to the reference
    # position and the fragment position as needed. If we hit the
    # fragment length (fragPos == fragEnd), then the cigarOffset might
    # be set to non-zero, which means the cigarIndex is not
    # incremented, so the same cigar entry will be used in the next
    # call to this function.
    while cigarIdx < len(read.cigartuples) and fragPos < fragEnd:
        opSize = operationSize(read.cigartuples[cigarIdx])
        if notDeletion(read.cigartuples[cigarIdx]):
            # not deletion in query, advance ref and frag
            fragRemaining = (fragEnd - fragPos)
            if opSize <= fragRemaining:
                if not isInsertion(read.cigartuples[cigarIdx]):
                    refEnd += opSize
                fragPos += opSize
                cigarIdx += 1
                cigarOffset = 0
            else:
                if not isInsertion(read.cigartuples[cigarIdx]):
                    refEnd += fragRemaining
                fragPos += fragRemaining
                cigarOffset += fragRemaining
        else:
            cigarIdx += 1
            refEnd += opSize # deletion in query, advance in reference

    return refEnd, cigarIdx, cigarOffset

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam-file', dest='samFileName',
                        required=True, help= 'mapped reads in SAM file')
    parser.add_argument('-l', '--frag-len', dest='fragLen', type=int,
                        required=True, help= 'actual fragment length')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print progress info to stderr')
    args = parser.parse_args()

    samFile = pysam.AlignmentFile(args.samFileName, 'r', check_sq=False)
    for read in samFile.fetch():
        # make sure the alignment is good
        if (not read.has_tag('XA') and
            not read.has_tag('SA') and
            not read.is_unmapped and
            read.mapping_quality > 0):

            fragStart = randint(0, args.fragLen-1) # first frag start
            fragEnd = fragStart + args.fragLen # first frag end
            refStart = read.reference_start # first ref start
            cigarIdx = 0
            cigarOffset = 0
            while fragEnd < len(read.query_alignment_sequence):
                refEnd, cigarIdx, cigarOffset = \
                    findFragEndInReference(read, fragStart, fragEnd,
                                           refStart, cigarIdx, cigarOffset)
                fragSeq = read.query_alignment_sequence[fragStart:fragEnd]
                if read.is_reverse: fragSeq = revcomp(fragSeq)
                print('%s\t%d\t%d\t%s\t0\t+' % (read.reference_name,
                                                refStart, refEnd, fragSeq))
                fragStart = fragEnd
                fragEnd += args.fragLen
                refStart = refEnd

if __name__ == '__main__':
    main()
