#!/bin/bash

module load cutadapt

gene=450X950

# demultiplex

mkdir -p demux_$gene
zcat raw/2015_03_13_*$gene*reads.fastq.gz | cutadapt -j 8 -O 28 -e 0.1 -g file:MIDs.fa -o demux_$gene/{name}.fastq.gz -

# remove primers

for MID in {1..75}
do
	mkdir -p trimmed_$gene
	zcat demux_$gene/MID-$MID[FR].fastq.gz | cutadapt -j 8 -g 'AATTACCCAATCCTGACACAGG' --rc -O 21 -m 300 --trimmed-only - 2> trimmed_$gene/MID-$MID.report.txt \
		| cutadapt -j 8 -a 'CCAGGGTAATGATTAATAGGGAC' -O 22 -m 300 --trimmed-only -o trimmed_$gene/MID-$MID.fastq.gz - >> trimmed_$gene/MID-$MID.report.txt
done
