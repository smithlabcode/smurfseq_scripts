## To test mapping a simulated SMURF-seq read to a
## small simulated genome
../map/smurfseq-map.sh map/ref_sim.fa map/read_sim.fastq

## Generate CNV profile
../cnv/cnv-5k.sh map/read_sim.fastq.2x.sam test_cnv

## To test mapping a simulated-read to the human
## genome
../map/smurfseq-map.sh [path to human genome] map/read_hg.fastq
