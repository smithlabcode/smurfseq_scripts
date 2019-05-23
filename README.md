# SMURF-seq scripts

SMURF-seq is a protocol for sequencing short reads on a long-read
sequencer by randomly concatenating short fragments. This repo
contains the scripts we used to conduct initial analysis of SMURF-seq
data in the context of copy-number profiling (with sequencing done on
the Oxford MinION instrument), and to benchmark mapping methods for
performance in mapping SMURF-seq reads.

## Copy-number analysis

The copy-number analysis we performed using SMURF-seq reads closely
follows procedures already used in other publications. We first map
the SMURF-seq reads using BWA:
```
bwa mem -x ont2d -k 12 -W 12 \
    -A 4 -B 10 -O 6 -E 3 -T 120 bwa-mem/index/hg19.fa \
    smurf_reads.fa > mapped_smurf_reads.sam
```
The parameters for the Smith-Waterman scoring ('A', 'B', 'O' and 'E')
were determined using the simulation approach outined below (see also
manuscript and supp info). The 'T' flag gives the minimum alignment
score to output. The 'k' gives the size of k-mers to use for
seeds. The 'W' indicates to discard a chain if seed bases are shorter
than this value. The 'k' and 'W' are set to be liberal to catch and
evaluate as many candidate mappings as possible.

The mapped fragments are given to a script that filters ambiguously
mapped fragments:
```
./filterAlnScoreAndQual.py -i mapped_smurf_reads.sam \
    -o unambig_smurf_frags.sam -s 120 -q 1
```

The input file `mapped_smurf_reads.sam` is just the mapped reads
(e.g. with BWA). The output file `unambig_smurf_frags.sam`
contains mapped fragments with mapping quality greater than or equal
to 1.

Then the remaining fragments are given to a script that obtains
the counts of reads in bins:
```
./getBinCounts.py -i unambig_smurf_frags.sam -c hg19.chrom.sizes \
    -b bins_5k_hg19.bed -o bin_counts.bed -s bin_stats.txt
```
The input file `unambig_smurf_frags.sam` is the same as described above. 
The file `hg19.chrom.size` is the size of all chroms
in the reference genome. This file for the hg19 reference is supplied
in the `data` directory in this repo, and was obtained from the UCSC
Genome Browser's database. The pre-defined bins file `bins_5k_hg19.bed`
is also in the `data` directory of this repo, and defines the 5000
bins in the genome used for the CNV analysis. The first output file
`bin_counts.bed` is like bedgraph format: it has the chrom, start and
end given in `bins_5k_hg1.bed` but it also has two extra columns:
one is the count of reads in that bin, and the final column is those
counts divided by the average reads per bin. This information was
determined based on what is required in the next script.

In the next step we use an adaptation of a script originally due to
Timour et al. (Nat. Protocols, 2014). The script is run
as follows:
```
./cnvAnalysis.R bin_counts.bed SampleName bins_5k_hg19_gc.txt bins_5k_hg19_exclude.txt
```
The input file `bin_counts.bed` is the same as described above. The
input file `bins_5k_hg19_gc.txt` is the GC content of each bin. The
input `bins_5k_hg19_exclude.txt` is used to exclude certain parts of
the genome that attract an unusual amount of reads.  The format is
simply the line numbers, in the corresponding bed file, of the bins to
exclude from the CNV analysis. The first output is a PDF file
`SampleName.pdf` for the CNV profile. In addition, two tables are
saved: one table `SampleName.data.txt` with the information
(chromosome, genome position, GC content, bin count, segmented value)
for each bin, and the other table `SampleName.short.txt` summerizing
the breakpoints in the CNV profile.

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
   `scripts` directory as follows:
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

## Dependencies

* [Python:] All Python scripts here are in Python3 (we used 3.6.8). The
    following non-standard libraries are used: pysam (we used 0.15.0) and
    numpy (we used 1.15.0).

* [R:] The R script `cnvAnalysis.R` uses the [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html) library.
    > Seshan VE, Olshen A (2018). DNAcopy: DNA copy number data analysis. R package version 1.56.0.

* [Software tools:] For the simulations/valuations we require
    `bedtools` (we used v2.26.0). We also require the `deadzones`
    program from `http://github.com/smithlabcode/utils` but this could
    be substituted for any means of excluding unmappable regions.
    In our CNV analysis, we used `bwa` (0.7.17).

## Contacts and Bug Reports

- Andrew D. Smith andrewds@usc.edu
- Rishvanth K. Prabakar kaliappa@usc.edu

## Copyright and License Information
Copyright (C) 2019 The Authors

Authors: Rishvanth K. Prabakar and Andrew D. Smith

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
