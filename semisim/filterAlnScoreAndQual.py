#!/usr/bin/env python3

# filterAlnScoreAndQual.py: filter mapping results in SAM format based
# on alignment score and quality
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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-sam', dest='inputFile',
                        required=True, help= 'mapped reads in SAM file')
    parser.add_argument('-o', '--output-sam', dest='outputFile',
                        required=True, help= 'output file name')
    parser.add_argument('-s', '--min-alignment-score', dest='minAlnScore',
                        type=float, required=True,
                        help= 'minimum alignment score to keep')
    parser.add_argument('-q', '--min-alignment-quality', dest='minAlnQual',
                        type=float, required=True,
                        help= 'minimum alignment quality to keep')
    args = parser.parse_args()

    samInput = pysam.AlignmentFile(args.inputFile, 'r', check_sq=False)
    samOutput = pysam.AlignmentFile(args.outputFile, 'w', template=samInput)

    for read in samInput.fetch():
        if (read.mapping_quality >= args.minAlnQual and
            read.has_tag('AS') and read.get_tag('AS') >= args.minAlnScore):
            samOutput.write(read)

    samInput.close()
    samOutput.close()


if __name__ == '__main__':
    main()
