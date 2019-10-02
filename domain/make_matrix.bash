#make query domain matrix

grep ">" query_renamed_100.pep | sed 's/>//g' | sort | uniq > query.all_del

cat *.pfam | sort | uniq | sort > all_pfam_del

while read p;do	
	echo $p >> id_dom_del

	if grep -w --quiet "$p" all_pfam_del
		then
    			grep -w $p all_pfam_del | awk '{print $7}' | sort | uniq | paste -sd";" | sed '/^$/d' >> dom_del
		else
    			echo "No_domain_detected"  >> dom_del
	fi

done < query.all_del

paste -d, id_dom_del dom_del > query_domain_del

#concatanate the query matrix and reference domain matrix
echo "pid,domain" >> header_del
cat reference_domain.matrix query_domain_del > domain_matrix_del
cat header_del domain_matrix_del > domain_final.matrix


#delete files
rm *_del
