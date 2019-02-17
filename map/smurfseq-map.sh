#!/bin/bash
###
#   Copyright (C) 2018 University of Southern California
#
#   Authors: Rishvanth Prabakar
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
###

if [ -z "$BWA" ]; then
  echo "Set path to BWA"
  exit 1
fi

if [ -z "$SAMTOOLS" ]; then
  echo "Set path to SAMTOOLS"
  exit 1
fi

if [ -z "$SMURFDIR" ]; then
  echo "Set path to SMURFDIR (smurfseq-scripts directory)"
  exit 1
fi

if [ $# -ne 2 ]; then
  echo "./smurfseq-map.sh [reference] [read]"
  exit 1
fi

$BWA mem -x ont2d -k 1 -W 5 -t 6 $1 $2 > $2.sam
$SMURFDIR/map/extractunmapped.py $2.sam $2 30 > $2.unmap.fastq
$BWA mem -x ont2d -k 1 -W 5 -t 6 $1 $2.unmap.fastq > $2.unmap.sam
$SAMTOOLS merge $2.2x.sam $2.sam $2.unmap.sam
rm $2.sam $2.unmap.sam $2.unmap.fastq
