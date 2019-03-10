#!/bin/bash
###     
#   Copyright (C) 2016 University of Southern California
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

if [ $# -ne 2 ]; then
  echo "./cnv-5k [sam file] [sample name]"
  exit 1
fi

$SMURFDIR/cnv/get_unique_maps.py $1 > $1.unique.sam
$SMURFDIR/cnv/get_bin_counts.py $1.unique.sam $SMURFDIR/cnv/5k/hg19.chrom.sizes.5k.txt $SMURFDIR/cnv/5k/hg19.bin.boundaries.5k.sorted.txt $2.varbin.5k.txt $2.varbin.5k.stats.txt
# R CMD BATCH '--args '$2.varbin.5k.txt' '$2' '$SMURFDIR/cnv/5k/hg19.new.varbin.gc.content.5k.txt' '$SMURFDIR/cnv/5k/hg19.badbins.5k.txt'' $SMURFDIR/cnv/cbs.r cbs.r.out
$SMURFDIR/cnv/cbs_smurf.r $2.varbin.5k.txt $2 $SMURFDIR/cnv/5k/hg19.new.varbin.gc.content.5k.txt $SMURFDIR/cnv/5k/hg19.badbins.5k.txt
