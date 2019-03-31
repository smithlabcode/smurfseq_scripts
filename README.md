# SMURF-seq

SMURF-seq is a protocol for sequencing short reads on a long-read
sequencer by randomly concatenating short fragments. This repo
contains the scripts we used to conduct initial analysis of SMURF-seq
data in the context of copy-number profiling (with sequencing done on
the Oxford MinION instrument), and to benchmark mapping methods for
performance in mapping SMURF-seq reads.

## Simulating SMURF-seq reads for evaluating mappers

The simulation procedure takes mapped long reads and uses them to
generate short fragments with a known mapping location. Since the
mapping location is known, we can assess a mapping strategy by how
well it recovers the known mapping locations. The steps are as follows.

1. Select and map long reads: The data set should be from WGS using
   the sequencing technology of interest (in our case, the Oxford
   MinION). Use a standard long-read mapper, and make sure the output
   is in SAM format.
2. Generate candidate short fragments: From the SAM format of the
   long-read mapping locations, generate candidate short fragments
   with known mapping locations. This is done with a script in the
   `scipts` directory as follows:
   ```
   $ ./getFragsFromLongReads.py -s long_reads_mapped.sam -l 100 > candidates.bed
   ```
   The candidates are encoded in 6-col BED format with the name column
   (the 4th) encoding the DNA sequence of the fragment. The parameter
   `-l` indicates the length of the candidate fragments to generate.
   Only uniquely mapping long reads are used for generating candidate
   fragments. The coordinates for each line in the `candidates.bed` file
   are the reference genome mapping coordinates for the short fragment
   obtained by arithmetic on the long-read reference location accounting
   for indels in the mapping.
3. Filter the candidates to exclude deadzones: Any short fragments that
   would not map uniquely without the rest of the long read are excluded
   by using bedtools and a file of deadzones. We obtain the deadzones using
   the program `deadzones` available from `http://github.com/smithlabcode/utils`
   and it will generate deadzones for a given read length. We used 40 bp, which
   is conservative if the candidate fragments are longer than 40 bp. The
   deadzone k-mer should not be larger than the candidate fragment size. The
   program is run as follows:
   ```
   $ ./deadzones -k 40 -o dz_hg19_40.bed hg19.fa
   ```
   This will use lots of memory and might be slow, but it has parameters
   to make it use less memory and it only needs to be done once.
   To filter the candidate fragments, use `bedtools` as follows:
   ```
   bedtools intersect -v -a candidates.bed -b dz_hg19_40.bed > good_frags.bed
   ```
   The `good_candidates.bed` will be used to generate the simulated
   SMURF-seq reads.
4. Randomly combine fragments into long reads: This step uses another script
   from the `scripts` directory:
   ```
   ./readsFromFrags.py -n FAB42704 -f 10 -l 100 -r 10000 \
       -b good_frags.bed -o simulated.fa
   ```
   Above, the `-n` parameter indicates an identifying for the original
   data set, which goes into the read name so we can map multiple
   simulations together but later take them apart and analyze the
   results separately. The `-f`, `-l` and `-r` give the number of
   fragments, their length and the number of reads to simulate (the
   fragment length could be inferred from the input). In the output,
   the name of each simulated SMURF-seq read in `simulated.fa`
   contains a list of the original mapping locations of each fragment
   in the reference genome, so later we can use these to evaluate mapping
   performance.
5. Map the simulated reads: Select a mapping tools, set the parameters, and
   map the simulated reads in `simulated.fa` from the above step to the
   same reference genome used already. The output should be in SAM format,
   and we will assume it is named `simulated_out.sam`.

6. Evaluate the mapping performance: here we use a script from the `scripts`
   directory:
   ```
   ./evalSmurfSim.py -f simulated.fa -s simulated_out.sam -l 100
   ```
   This script computes precision and recall at the level of recovering
   fragments and at the level of correctly determining mapping location
   for each nucleotide in the simulated reads. In addition, the portion of
   each simulated read that is mapped in at least one identified fragment
   is reported, along with summary stats on the sizes of fragments identified
   by the mapper.

## Python

All Python scripts here are in Python3. The following non-standard
libraries are used: pysam and numpy.

## R

2. DNAcopy (R)

## Software tools

1. BWA
2. bedtools
2. samtools (1.9)

## SMURF-seq Scripts
1. Map the reads:

  ./map/smurfseq-map.sh [reference] [reads]

## Contacts and Bug Reports
Andrew D. Smith andrewds@usc.edu

Rishvanth K. Prabakar kaliappa@usc.edu

## Copyright and License Information
Copyright (C) 2018 University of Southern California

Authors: Rishvanth K. Prabakar

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
