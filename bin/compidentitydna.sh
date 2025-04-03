#!/bin/bash

set -e
set -o pipefail

if [ "$#" -ne 3 ]; then
    echo "usage : $0 genome.fa/genome.idx reads1.fastq/reads1.bam reads1.fastq/reads1.bam"
        exit
fi

die() {
    local msg="${1}"
    echo "Error: ${msg}" >&2
    exit 1
}

test -e $1 || die "$1 does not exist"

FASTQ1=$2
test -e $FASTQ1 || die "$FASTQ1 does not exist"
extension1="${FASTQ1##*.}"
base1=$(basename $FASTQ1)
name1=${base1%%.*}
timestamp1=$(date +%Y%m%d_%H%M%S)
output1="${name1}_identity_scores_${timestamp1}.txt"
test -e $output1 && die "$output1 already exists, please remove it first"

FASTQ2=$3
test -e $FASTQ2 || die "$FASTQ2 does not exist"
extension1="${FASTQ2##*.}"
base2=$(basename $FASTQ2)
name2=${base2%%.*}
timestamp2=$(date +%Y%m%d_%H%M%S)
output2="${name2}_identity_scores_${timestamp2}.txt"
test -e $output2 && die "$output2 already exists, please remove it first"

OUTPUT="${name1}_${timestamp1}_vs_${name2}_${timestamp2}.txt"
test -e $OUTPUT && die "$OUTPUT already exists, please remove it first"

minimap2 --version > /dev/null || die "minimap2 not found in PATH"
datamash --version > /dev/null || die "datamash not found in PATH"
samtools --version > /dev/null || die "datamash not found in PATH"

if [ "$extension1" == "bam" ] ||  [ "$extension1" == "sam" ] ; then
    samtools fastq ${FASTQ1} | minimap2 -ax map-ont $1 -t8 --secondary=no - | samtools view -F 2308 -h - | paftools.js sam2paf - | awk '{print $1"\t"$10/$11}' | sort -k1,1 > ${output1}
else
    minimap2 -ax map-ont $1 -t8 --secondary=no ${FASTQ1} | samtools view -F 2308 -h - | paftools.js sam2paf - | awk '{print $1"\t"$10/$11}' | sort -k1,1 > ${output1}
fi
test -s "$output1" || die "Mapping failed"

if [ "$extension2" == "bam" ] ||  [ "$extension2" == "sam" ] ; then
    samtools fastq ${FASTQ2} | minimap2 -ax map-ont $1 -t8 --secondary=no - | samtools view -F 2308 -h - | paftools.js sam2paf - | awk '{print $1"\t"$10/$11}' | sort -k1,1 > ${output2}
else
    minimap2 -ax map-ont $1 -t8 --secondary=no ${FASTQ2} | samtools view -F 2308 -h - | paftools.js sam2paf - | awk '{print $1"\t"$10/$11}' | sort -k1,1 > ${output2}
fi
test -s "$output2" || die "Mapping failed"

join ${output1} ${output2} | awk '{print $1"\t"$2"\t"$3}' > ${OUTPUT} || die "join failed"
test -s "$OUTPUT" || die "join failed"

CORR=$(cat ${OUTPUT} | awk '{print $2"\t"$3}' | datamash -W ppearson 1:2)
test -z "$CORR" && die "datamash for correlation failed"

echo "Correlation: ${CORR}"
echo "Joined file: ${OUTPUT}"
