grep ">" query_renamed_100.pep | sed 's/>//g' | sort | uniq > query.all

cat *.all > id_list_del

cat *mito_prot | sort | uniq > mito_prot_del

cat *.pfam | sort | uniq | sort > all_pfam_del

while read p;do
	
	echo $p >> id_dom_del

	if grep -w --quiet "$p" mito_prot_del
		then
    			echo "mito_dom"     >> mdom_del
		else
    			echo "nonmito_dom"  >> mdom_del
	fi
	
	
	if grep -w --quiet "$p" all_pfam_del
		then
    			grep -w $p all_pfam_del | awk '{print $7}' | sort | uniq | paste -sd";" | sed '/^$/d' >> dom_del
		else
    			echo "No_domain_detected"  >> dom_del
	fi
	
done < id_list_del

echo "ID,Mito_domain,Domain_combination" >> file_del

paste -d, id_dom_del mdom_del dom_del > domain.matrix_del

cat file_del domain.matrix_del > domain.matrix

rm *_del


