#makes the quert matrix and combines it with the reference matrix

grep ">" query_renamed_100_Mlist.pep | sed 's/>//g' | sort | uniq > id_del

while read p;do
	grep -w $p query_renamed_100_Mlist.pep.targetP   | awk '{print $6","$7}' >> 1_del
	grep -w $p query_renamed_100_Mlist.pep.mitofates.result | awk -F"\t" '{print $2","$3}' >> 2_del
	grep -w $p query_renamed_100_Mlist.pep.mitofates.result | awk -F"\t" '{print $4}' | awk -F"(" '{print $1}' >> 3_del
done < id_del
paste -d, id_del 1_del 2_del 3_del > query_mts.matrix_del

cp ../orthology/query_renamed_100.pep .
grep ">" query_renamed_100.pep | sed 's/>//g' | sort | uniq > ids

while read p;do
	if grep -w --quiet $p query_mts.matrix_del
		then
			grep -w $p query_mts.matrix_del | sed 's/ /_/g' >> 8_del
		else
			echo "$p,N,N,N,N,N" >> 8_del
		fi
done < ids

cat reference_mts.matrix 8_del > final_mts.matrix

rm *_del
rm query_renamed_100.pep
rm ids


