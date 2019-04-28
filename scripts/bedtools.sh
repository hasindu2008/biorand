#!/bin/bash

FILENAME=candidates

function get_gc_stat {
	bedtools sort -i $FILENAME.bed > $FILENAME.sorted.bed
	awk '{print $1":"$2"-"$3}' $FILENAME.sorted.bed > regions.list
	samtools faidx /mnt/d/genome/hg38noAlt.fa -r regions.list > extracted_ref.fa
	seqtk comp extracted_ref.fa > gcstat.tsv
	rm extracted_ref.fa regions.list
	#chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts

	cat  /mnt/d/genome/hg38noAlt.fa.fai | awk '{print $1"\t0\t"$2}' > hg38noAlt.bed
	bedtools sort -i hg38noAlt.bed > hg38noAlt.sorted.bed
	awk '{print $1"\t"$3}' hg38noAlt.sorted.bed > hg38noAlt.bed
	rm hg38noAlt.sorted.bed
	bedtools complement -i $FILENAME.sorted.bed -g hg38noAlt.bed > complement.bed
	awk '{print $1":"$2"-"$3}' complement.bed > complement.list

	samtools faidx /mnt/d/genome/hg38noAlt.fa -r complement.list > extracted_complement.fa
	seqtk comp extracted_complement.fa > gcstatcomplement.tsv
	rm extracted_complement.fa complement.list hg38noAlt.bed
}


function get_rmsk_overlap {

	REP=/mnt/d/genome/rmsk.bed

		
	echo "Intersection with repeat mask"
	bedtools intersect -a $FILENAME.bed -b $REP -wb  -loj > t_with_rmsk.tsv
	echo "uniq intersections"
	bedtools intersect -a $FILENAME.bed -b $REP -u | wc -l
	echo "no intersections"
	awk '{if($10==".") print $0}'  t_with_rmsk.tsv | wc -l
	awk '{print $10}' t_with_rmsk.tsv | sort | uniq -c | sort -nrk1,1 > rep_hist.txt
	echo ""

}


function get_rmsk_class_overlap {

	echo "Intersection with repeat mask class"
	REP=/mnt/d/genome/rmsk_class.bed
	bedtools intersect -a $FILENAME.bed -b $REP -wb  -loj > t_with_rmsk_class.tsv
	echo "uniq intersections"
	bedtools intersect -a $FILENAME.bed -b $REP -u | wc -l
	echo "no intersections"
	awk '{if($10==".") print $0}'  t_with_rmsk_class.tsv | wc -l
	awk '{print $10}' t_with_rmsk_class.tsv | sort | uniq -c | sort -nrk1,1 > rep_hist_class.txt

	
	echo "Intersection of complement with repeat mask class"
	REP=/mnt/d/genome/rmsk_class.bed
	bedtools intersect -a complement.bed -b $REP -wb  -loj > complements_with_rmsk_class.tsv
	echo "uniq intersections"
	bedtools intersect -a complement.bed -b $REP -u | wc -l
	echo "no intersections"
	awk '{if($7==".") print $0}'  complements_with_rmsk_class.tsv | wc -l
	awk '{print $7}' complements_with_rmsk_class.tsv | sort | uniq -c | sort -nrk1,1 > complement_rep_hist_class.txt
	echo ""
}

#get number of bases overlapped
function get_rmsk_class_overlap_bases {
	bedtools intersect -a $FILENAME.bed -b $REP  -wao > t_with_rmsk_class_bases.tsv
	bedtools intersect -a complement.bed -b $REP  -wao > complement_with_rmsk_class_bases.tsv
}

function get_odds_ratio {

	REPTYPE=$1
	martian_reps=$(awk -v var="$REPTYPE" '{if($10==var) sum=sum+$NF}END{print sum}' t_with_rmsk_class_bases.tsv)
	martian_sum=$(awk '{sum=sum+($3-$2)}END{print sum}' $FILENAME.bed)
	martian_non_reps=$(echo "$martian_sum-$martian_reps" | bc)

	non_martian_reps=$(awk -v var="$REPTYPE" '{if($7==var) sum=sum+$NF}END{print sum}' complement_with_rmsk_class_bases.tsv)
	non_martian_sum=$(awk '{sum=sum+($3-$2)}END{print sum}' complement.bed)
	non_martian_non_reps=$(echo "$non_martian_sum-$non_martian_reps" | bc)

	echo "repeat type : $1"
	echo "Martian regions reps : $martian_reps"
	echo "Martian regions total sum : $martian_sum"
	echo "Non Martian regions reps : $non_martian_reps"
	echo "Non Martian regions total sum : $non_martian_sum"
	echo -n "odds ratio : "
	echo "($martian_reps/$non_martian_reps)/($martian_non_reps/$non_martian_non_reps)" | bc -l
	
	echo ""
}

get_gc_stat
get_rmsk_overlap
get_rmsk_class_overlap
get_rmsk_class_overlap_bases
get_odds_ratio "Simple_repeat"
get_odds_ratio "SINE"
get_odds_ratio "LINE"
get_odds_ratio "LTR"
get_odds_ratio "Satellite"
get_odds_ratio "Low_complexity"
get_odds_ratio "Retroposon"



