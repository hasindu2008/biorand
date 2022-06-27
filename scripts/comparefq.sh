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
	echo "Usage: ${0} guppy_basecall_dir1 guppy_basecall_dir2"
	exit
fi

A=$1
B=$2

test -e $A || die "$A not present."
test -e $B || die "$B not present."

cat $A  | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > a.fastq
cat $B  | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > b.fastq
diff -q a.fastq b.fastq || die "fastq differ"

