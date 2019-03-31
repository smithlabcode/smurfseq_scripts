#!/usr/bin/env python3

# smurfSimDemux.py: take apart the different data sets and simulated
# fragment lengths from mapped simulated SMURF-seq reads (so they can
# be analyzed separately more easily)
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
import pysam

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in-sam-file', dest='inSamFile', required=True,
                        help= 'input reads in SAM format')
    parser.add_argument('-p', '--out-prefix', dest='outPrefix',
                        required=True, help= 'output file prefix')
    parser.add_argument('-s', '--out-suffix', dest='outSuffix',
                        required=True, help= 'output file suffix (before .sam)')
    parser.add_argument('-d', '--out-dir', dest='outDir',
                        required=True, help= 'output directory')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print progress info to stderr')
    args = parser.parse_args()

    outFiles = {}
    mapfile = pysam.AlignmentFile(args.inSamFile, 'r', check_sq=False)
    for read in mapfile.fetch():
        dataSetFrag = '_'.join(read.query_name.split('_')[:2])
        if dataSetFrag not in outFiles:
            fn = os.path.join(args.outDir,
                              args.outPrefix + '_' + dataSetFrag +
                              args.outSuffix + '.sam')
            if args.VERBOSE:
                print('creating output file: %s' % fn, file=sys.stderr)
            outFiles[dataSetFrag] = \
                pysam.AlignmentFile(fn, 'w', template=mapfile)
        outFiles[dataSetFrag].write(read)

if __name__ == '__main__':
  main()
