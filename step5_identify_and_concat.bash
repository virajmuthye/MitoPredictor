#run this as ./step5_identify_and_concat.bash <dom cutoff used> <cutoff for mitochondrial proteins>
mkdir concatanate
cd concatanate

#copy matrices
cp ../orthology/orthology.matrix .
cp ../orthology/query_renamed_100.pep .
cp ../domain/domain_final.matrix .
cp ../mts/final_mts.matrix .

#copy list of all mitochondrial proteins
cp ../update_final_matrix/all_mito_list .

#identify mitochondrial proteins
grep ">" query_renamed_100.pep | sed 's/>//g' | sort | uniq | sort > query_ids_del

##################################################################
#ortho score
##################################################################
#identify all mitochondrial OGs
while read p;do
	grep -w $p orthology.matrix | awk -F, '{print $1}' >> mito_ogs_del
done < all_mito_list
wait
sort mito_ogs_del | uniq | sort > mito_ogs_uniq_del

#find all query proteins from those OGs
#orthology approach results are in qery_ortho_sorted_del
while read p;do
	grep -w $p orthology.matrix | grep "qery" | awk -F, '{print $2}' >> qery_ortho_del
done < mito_ogs_uniq_del
sort qery_ortho_del | uniq > qery_ortho_sorted_del

#now, assign ortho score
while read p;do
	if grep -w --quiet $p qery_ortho_sorted_del
		then
			echo "1" >> ortho_score_del
		else
			echo "0" >> ortho_score_del
	fi
done < query_ids_del

##################################################################
#mts score
##################################################################
grep "qery" final_mts.matrix > query_mts_matrix_del
grep ",M," query_mts_matrix_del | awk -F, '{print $1}' > m_matrix_list_del
grep -v ",M," query_mts_matrix_del | awk -F, '{print $1}' > nm_matrix_list_del

while read p;do

	if grep -w --quiet $p m_matrix_list_del
		then
			grep -w $p query_mts_matrix_del  | awk -F, '{print 6-$3+$4}' >> mts_score_del
		else
			grep -w $p  query_mts_matrix_del | awk -F, '{print 0+$4}' >> mts_score_del
	fi

done < query_ids_del

##################################################################
#domain score
##################################################################
# need the cut-off used for domain analysis
grep ",mito_dom" domain_final.matrix | awk -F, '{print $1}' | sort | uniq > dom_list_del

while read p;do
	if grep -w --quiet $p dom_list_del
		then
			echo $1 >> domain_score_del
		else
			echo "0" >> domain_score_del
	fi
done < query_ids_del

##################################################################
#find overall mitoscore
##################################################################
paste -d, query_ids_del ortho_score_del mts_score_del domain_score_del > mito_score_matrix_del
while read p;do
	grep -w $p mito_score_matrix_del | awk -F, '{print $2+$3+$4}' >> total_score_del
done < query_ids_del

paste -d, mito_score_matrix_del total_score_del > mito_score_matrix_del_2
echo "pid,orthoscore,mtsscore,domainscore,totalscore" > header_mito_del
cat header_mito_del mito_score_matrix_del_2 > mito_score_matrix
rm *_del

awk -F, -v a=$2 '$5>=a {print $1}' mito_score_matrix | grep -v "pid" > query.mito.list

cat all_mito_list query.mito.list > all_mito_list_updated

##################################################################
#concatanate the matrices
##################################################################
for file in *.matrix;do
	head -1 $file > "$file".header_del
	tail -n +2 $file | sort -t, -k1 -n > "$file"_del 
	cut -d, -f2- "$file"_del > "$file"_truncated_del 
done

awk -F, '{print $1}' final_mts.matrix_del > all_del

while read p;do
	if grep -w --quiet $p all_mito_list_updated
		then
			echo "mito" >> mito_del
		else 
			echo "nonmito" >> mito_del
	fi

	if grep -w --quiet $p orthology.matrix
		then
			grep -w $p orthology.matrix | awk -F, '{print $1","$3}' >> og_del
		else 
			echo "no_og,no_og" >> og_del
	fi

	if grep -w --quiet $p domain_final.matrix
		then
			grep -w $p domain_final.matrix | cut -d, -f2- >> domain_del
		else 
			echo "no_domain" >> domain_del
	fi

	if grep -w --quiet $p final_mts.matrix
		then
			grep -w $p final_mts.matrix | cut -d, -f2- >> mts_del
		else 
			echo "no_mts" >> mts_del
	fi


done < all_del

paste -d, all_del mito_del og_del mts_del domain_del > FINAL.MATRIX_almost_del
awk -F, '{print $1}' FINAL.MATRIX_almost_del | sed 's/[0-9]//g' > species_del
paste -d, species_del FINAL.MATRIX_almost_del > FINAL.MATRIX_del
echo "species,pid,mito,og,geneid,tp,rc,mfprob,mf,mtslen,domcomb,mitodom" >> header_del
cat header_del FINAL.MATRIX_del > FINAL.MATRIX
sed -i 's/ /_/g' FINAL.MATRIX
cp FINAL.MATRIX ../.
cp mito_score_matrix ../.
rm *_del



  






