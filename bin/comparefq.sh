#!/bin/bash

set -e
set -o pipefail

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

skip=0

# Get the options
while getopts ":s" option; do
   case $option in
      s) # 
         skip=1
   esac
done
shift $(expr $OPTIND - 1 )

if [ $# -ne 2 ]
then
        echo "Usage: ${0} [OPTIONS] a.fastq b.fastq"
        echo "-s        do not compare fastq header comments"
        exit
fi


A=$1
B=$2

test -e $A || die "$A not present."
test -e $B || die "$B not present."

test -e a.fastq && die "temporary file a.fastq is present. Remove that first!"
test -e b.fastq && die "temporary file b.fastq is present. Remove that first!"

if [ ${skip} == "1" ];
then
	echo "FASTQ header comments will not be compared."
	cat $A  | awk '{if(NR%4==1){print $1}else{print$0}}' | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > a.fastq
	cat $B  | awk '{if(NR%4==1){print $1}else{print$0}}' | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > b.fastq
else
	cat $A  | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > a.fastq
	cat $B  | paste - - - -  | sort -k1,1 -T . | tr '\t' '\n' > b.fastq
fi 
diff -q a.fastq b.fastq || die "fastq differ"

rm -f a.fastq b.fastq
