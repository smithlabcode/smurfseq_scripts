# SMURF-seq
SMURF-seq is a protocol for sequencing short reads on a long-read
sequencer. This repo contains the scripts to map the reads from
a SMURF-seq run.

## Library Dependencies
1. pysam (python)
2. DNAcopy (R)

## Software Dependencies
1. BWA
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
