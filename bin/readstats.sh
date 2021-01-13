#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "usage : $0 in.fastq/in.fastq.gz/in.fasta/in/fasta.gz"
        exit
fi


die() {
    local msg="${1}"
    echo "Error: ${msg}" >&2
    exit 1
}

filename=$1
extension="${filename##*.}"
CAT=cat

if [ $extension = "gz"  ]
then
    extension=`basename $filename .gz`
    extension="${extension##*.}"
    CAT=zcat

fi

if [ $extension = "fastq"  -o $extension = "fq" ]
then
    #Get number of bases from a FASTQ
    echo "Total bases:"
    $CAT $filename | awk 'BEGIN{sum=0}{if(NR%4==2) {sum=sum+length($0); }}END{print sum}'

    #Get Number of reads from a fastq
    echo "Total reads:"
    echo "$($CAT $filename | wc -l)/4" | bc

    #Get mean and max rlends from a fastq
    echo "mean and max read lengths"
    $CAT $filename | awk '{if(NR%4==2) {print length($0)}}'  | datamash mean 1 max 1

elif [ $extension = "fasta"  -o $extension = "fa" ]
then
    #Get number of bases from a fasta
    echo "Total bases:"
    $CAT $filename | grep -v ">"  | wc | awk '{print $3-$1}'

    #Get number of reads from a fasta
    echo "Total reads:"
    $CAT $filename | grep ">"  | wc -l

else
    die "Unknown file extension $extension"
fi
