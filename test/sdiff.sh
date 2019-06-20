for w in test1 test2;  
do    
	samtools view -F 2308 $w.sam | awk '{print $1,$3,$4}' > sam.$w.tmp
done

sdiff sam.test1.tmp sam.test2.tmp -w 2048 | grep "|" > mismatch.sdiff ; 


matched=$(sdiff sam.test1.tmp sam.test2.tmp -w 2048 | grep -v "|\|<\|>" | wc -l);
mismatches=$(sdiff sam.test1.tmp sam.test2.tmp -w 2048 | grep "|" | wc -l);     
otheronly=$(sdiff sam.test1.tmp sam.test2.tmp -w 2048 | grep ">" | wc -l);     
goldonly=$(sdiff sam.test1.tmp sam.test2.tmp -w 2048 | grep "<" | wc -l); 


#awk 'begin{}{if($1==$5 && $2==$6 && $3==$7) print "error"}' mismatch.$w.sdiff

#perform error correction
error=$(awk '{if($1!=$5) print $1,$5}' mismatch.sdiff | wc -l)
otheronly=$(echo "$otheronly+$error" | bc)
goldonly=$(echo "$goldonly+$error" | bc)
mismatches=$(echo "$mismatches-$error" | bc)


#error check
entries_in_gold=$(cat sam.test1.tmp | wc -l)
entries_in_sam=$(cat sam.test2.tmp | wc -l)
computed_entries_in_gold=$(echo "$matched+$mismatches+$goldonly" | bc )
computed_entries_in_sam=$(echo "$matched+$mismatches+$otheronly" | bc )

#each should be equal
[ "$computed_entries_in_gold" -ne "$entries_in_gold" ] && echo "Error" && exit
[ "$computed_entries_in_sam" -ne  "$entries_in_sam" ] && echo "Error" && exit

echo "mismatches test1 multipartmerge"
echo "$mismatches $goldonly $otheronly";

