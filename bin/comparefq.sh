#!/bin/bash

set -e
set -o pipefail

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ $# -ne 2 ]
then
	echo "Usage: ${0} a.fastq b.fastq"
	exit
fi

A=$1
B=$2

test -e $A || die "$A not present."
test -e $B || die "$B not present."

test -e a.fastq && die "temporary file a.fastq is present. Remove that first!"
test -e b.fastq && die "temporary file b.fastq is present. Remove that first!"

cat $A  | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > a.fastq
cat $B  | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > b.fastq
diff -q a.fastq b.fastq || die "fastq differ"

rm -f a.fastq b.fastq
