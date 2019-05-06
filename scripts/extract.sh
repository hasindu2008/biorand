#!/bin/bash

awk '{print $4}' inserts.bed > readid.txt

for file in fast5/LXBAB132606/1-E11-H11/reads/*/*; do echo $file; done > fast5fetcher.index
gzip fast5fetcher.index

python ~/installs/fast5_fetcher/fast5_fetcher.py -f inserts/readid.txt -s  LXBAB132606_summary.txt -i fast5fetcher.index.gz -o ./inserts/

 zcat ../all.fastq.gz | grep -A3 -F -f readid.txt  --no-group-separator > inserts.fastq

awk '{if(NR%4==1){printf $1":"$3"-"}; if(NR%4==2){ printf $2"\n"}}' pairs.tsv > reads.list

samtools faidx inserts.fastq -r reads.list > reads.extracted.faidx

	# awk '{print $1":"$2"-"$3}' inserts.bed > regions.list
	# samtools faidx /mnt/d/genome/hg38noAlt.fa -r regions.list > extracted_ref.fa