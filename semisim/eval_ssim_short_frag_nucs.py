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
  def __init__(self, chrom, pos, length):
    self.chrom = chrom
    self.pos = int(pos)
    self.length = int(length)
  def matches(self, other):
    """r = self, m = other"""
    if self.chrom != other.chrom:
      return False
    else:
      overlap = 0
      if other.pos < self.pos:
        overlap = other.pos + other.length - self.pos
      elif (other.pos + other.length < self.pos + self.length):
        overlap = other.length
      else:
        overlap = self.pos + self.length - other.pos
      return overlap >= self.length/2

def overlaps(a, b):
  """two list of 2 elements each: a = [x,y] and b = [x',y']"""
  return ((b[0] <= a[0] and a[0] < b[1]) or
          (a[0] <= b[0] and b[0] < a[1]))

def countQueryBasesMapped(mappedFrags):
  """
  mappedFrags is a dictionary keyed on the name of a read and
  containing a list of ordered pairs (start, end) for the identified
  fragment locations within the read
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

def loadPredictedRefBases(mappedSamFile):
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

def getStartPosFromCigar(samRead):
  """get the mapping start position"""
  if samRead.is_reverse: # reverse strand
    if (samRead.cigartuples[-1][0] == 4 or samRead.cigartuples[-1][0] == 5):
      return samRead.cigartuples[-1][1] - 1 # fix sam format 1-offset
  else: # forward strand
    if (samRead.cigartuples[0][0] == 4 or samRead.cigartuples[0][0] == 5):
      return samRead.cigartuples[0][1] - 1 # fix sam format 1-offset
  return 0

def loadPredictedQueryLocations(mappedSamFile):
  """
  takes a SAM file name and returns a dictionary with keys
  corresponding to the reads and locations corresponding to the
  start and end positions of the mappings within that query read
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
  takes a SAM file name and returns a dictionary with keys
  corresponding to reads and for each read, a list of the intervals in
  the reference genome where parts of that read are predicted to map.
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

def loadActualRefLocations(fastaFileName, fragLen):
  mapLoc = {}
  for line in open(fastaFileName):
    if line[0] == '>': # process only name lines
      name = line.strip().split()[0][1:]
      readPosPart = line.strip().split()[1]
      mapLoc[name] = [Interval(x.split(':')[0], x.split(':')[1], fragLen)
                      for x in readPosPart.split(',')]
  return mapLoc

def getActualRefBases(actualRefLocations):
  refBases = {}
  for readName in actualRefLocations:
    for r in actualRefLocations[readName]:
      if r.chrom not in refBases:
        refBases[r.chrom] = {}
      for i in range(r.pos, (r.pos + r.length)):
        if i not in refBases[r.chrom]:
          refBases[r.chrom][i] = 0
        refBases[r.chrom][i] += 1
  return refBases

def getNucsAccuracy(actual, predicted):
  tpNucs = 0
  fpNucs = 0
  for i in predicted:
    if i in actual:
      for j in predicted[i]:
        if j in actual[i]:
          tpNucs += min(actual[i][j], predicted[i][j])
          if predicted[i][j] > actual[i][j]:
            fpNucs += (predicted[i][j] - actual[i][j])
        else: fpNucs += predicted[i][j]
    else: fpNucs += sum(predicted[i].values())

  fnNucs = 0
  for i in actual:
    if i in predicted:
      for j in actual[i]:
        if j in predicted[i]:
          if actual[i][j] > predicted[i][j]:
            fnNucs += actual[i][j] - predicted[i][j]
        else: fnNucs += actual[i][j]
    else: fnNucs += sum(actual[i].values())
  return tpNucs, fpNucs, fnNucs

def main():

  ### Parse the arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', '--fasta-file', dest='fastaFileName', required=True,
                      help= 'reads in FASTA format')
  parser.add_argument('-s', '--sam-file', dest='samFileName',
                      required=True, help= 'mapped reads in SAM file')
  parser.add_argument('-l', '--frag-len', dest='fragLen', type=int,
                      required=True, help= 'actual fragment length')
  args = parser.parse_args()


  ########################################################################
  # load the data from the sam file, both the predicted query
  # fragments and their predicted mapping locations in the reference
  # genome
  predictedQueryFrags = loadPredictedQueryLocations(args.samFileName)
  predictedRefLocations = loadPredictedRefLocations(args.samFileName)

  # predictedRefFragLengths = []
  # for i in predictedRefLocations:
  #   predictedRefFragLengths += [j.length for j in predictedRefLocations[i]]

  ########################################################################
  ## process the true simulated read locations and determine the
  ## counts for those that map correctly, those that are not mapped,
  ## and the mapping locations that are false positives
  actualRefLocations = loadActualRefLocations(args.fastaFileName, args.fragLen)

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
  predictedRefBases = loadPredictedRefBases(args.samFileName)
  actualRefBases = getActualRefBases(actualRefLocations)

  tpNucs, fpNucs, fnNucs = getNucsAccuracy(actualRefBases, predictedRefBases)
  predictedPositiveNucs = tpNucs + fpNucs
  conditionPositiveNucs = tpNucs + fnNucs
  precisionNucs = float(tpNucs)/predictedPositiveNucs
  recallNucs = float(tpNucs)/conditionPositiveNucs

  ########################################################################
  ## determine the fraction of nucleotides in reads that are mapped
  nTotalQueryBases = sum([len(i.strip())
                          for i in open(args.fastaFileName) if i[0] != '>'])
  nQueryBasesMapped = countQueryBasesMapped(predictedQueryFrags)
  unmappedBases = nTotalQueryBases - nQueryBasesMapped
  fracBasesUnmapped = float(unmappedBases)/nTotalQueryBases

  ########################################################################
  ## format and print the results
  outputFormat = ['%.4f', # precision
                  '%.4f', # recall
                  '%.4f', # precision nuc
                  '%.4f', # recall nuc
                  '%.3f', # bases unmapped
                  '%.3f', # mean true pos frag len
                  '%.3f', # sd of true pos frag len
                  '%d']   # median true pos frag len

  outputFormat = '\t'.join(outputFormat)
  print(outputFormat %
        (precision,
         recall,
         precisionNucs,
         recallNucs,
         fracBasesUnmapped,
         meanTruPosFragLen,
         sdTruPosFragLen,
         medianTruPosFragLen))

if __name__ == '__main__':
  main()
