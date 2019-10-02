grep ">" query_renamed_100.pep | sed 's/>//g' | sort | uniq > query_ids_del

while read p;do
	grep -w $p query.pfam | awk '{print $1}' >> query_ids_with_domains_with_score_del
done < domains_with_score

sort query_ids_with_domains_with_score_del | uniq > query_ids_with_domains_with_score_uniq_del

while read p;do
	grep -w $p query.pfam >> query.pfam.filtered_del	
done < query_ids_with_domains_with_score_uniq_del

while read p;do
	
	echo $p >> domscore1_del
	
	if grep -w --quiet "$p" query.pfam.filtered_del
		then
    			grep -w $p query.pfam | awk '{print $6}' > temp1
			
			while read q;do
				grep -w $q score.matrix | awk -F, '{print $2}' >> temp2 
			done < temp1

			sort -nr temp2 | head -1 >> domscore2_del

			rm temp1
			rm temp2

		else
    			echo "0"  >> domscore2_del
	fi
	
done < query_ids_del

paste -d, domscore1_del domscore2_del > domscore.matrix
rm *del


