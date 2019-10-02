mkdir concatanate
cd concatanate

#copy matrices
cp ../prep/id_renamedid.map .
cp ../orthology/orthology.matrix .
cp ../orthology/query_renamed_100.pep .
cp ../domain/domain_final.matrix .
cp ../domain/domscore.matrix .
cp ../mts/final_mts.matrix .
cp ../mts/*mitofates.result .
cp ../mts/*targetP .
cp ../rf.R .
cp ../all.training.csv .
cp ../mts/query_renamed_100_Mlist.pep .
cp ../gene.csv .

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
grep ">" query_renamed_100_Mlist.pep | sed 's/>//g' | sort | uniq > query_renamed_100_Mlist.pep_del


while read p;do

	if grep -w --quiet $p query_renamed_100_Mlist.pep_del
		then
			grep -w $p query_renamed_100_Mlist.pep.targetP | awk '{print $3","$4","$5","$6","$7}' >> targetp_del
			grep -w $p query_renamed_100_Mlist.pep.mitofates.result | awk -F"\t" '{print $2","$5}' >> mf_del
		else
			echo "0.000,0.000,0.000,0,0" >> targetp_del
			echo "0.000,0.000" >> mf_del	
	fi

done < query_ids_del

sed -i 's/,S/,0/g' targetp_del
sed -i 's/,C/,0/g' targetp_del
sed -i 's/,_/,0/g' targetp_del
sed -i 's/,M/,1/g' targetp_del

##################################################################
#domain score
##################################################################
while read p;do
	if grep -w --quiet $p domscore.matrix
		then
			grep -w $p domscore.matrix | awk -F, '{print $2}' | bc -l | xargs printf "%.3f\n" >> domain_score_del
		else
			echo "0.000" >> domain_score_del
	fi
done < query_ids_del

##################################################################
#make matrix
##################################################################
paste -d, query_ids_del ortho_score_del targetp_del mf_del domain_score_del > mito_score_matrix_del
echo "pid,orthoscore,mtp,stp,other,pred,rc,mfprob,charge,domainscore" > header_mito_del
cat header_mito_del mito_score_matrix_del > mito_score_matrix_temp
rm *_del

awk -F, '{print $1","$9","$2","$8","$10","$3","$4","$5","$6","$7}' mito_score_matrix_temp > mito_score_matrix

rm mito_score_matrix_temp

Rscript rf.R

sed -i 's/"//g' rf.csv 

awk -F, '$12==1 {print $2}' rf.csv | sort | uniq > query.mito.list 

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
			echo "og" >> og2_del
			grep -w $p orthology.matrix | awk -F, '{print $1","$3}' >> og_del
		else 
			echo "no_og" >> og2_del
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

	if grep -w --quiet $p id_renamedid.map
		then
			grep -w $p id_renamedid.map | cut -d, -f1 >> mapcsv_2_del
		else 
			echo "nd" >> mapcsv_2_del
	fi

done < all_del

paste -d, all_del mapcsv_2_del mito_del og2_del og_del mts_del domain_del > FINAL.MATRIX_almost_del

awk -F, '{print $1}' FINAL.MATRIX_almost_del | sed 's/[0-9]//g' > species_del

paste -d, species_del FINAL.MATRIX_almost_del > FINAL.MATRIX_del

echo "species,pid,id,mito,og_present,og,geneid,tp,rc,mfprob,mf,mtslen,domcomb" >> header_del

cat header_del FINAL.MATRIX_del > FINAL.MATRIX

sed -i 's/ /_/g' FINAL.MATRIX

rm *_del


##################################################################

awk -F, '{print $1}' gene.csv > gene_temp1_del

while read p;do
	grep -w $p FINAL.MATRIX | awk -F, '{print $5}' >> gene_temp2_del
done < gene_temp1_del

paste -d, gene.csv gene_temp2_del > gene.map.csv_del

echo "renamed_proteinid,proteinid,og_number" > gene_header_del

cat gene_header_del gene.map.csv_del > gene.map.csv

cp FINAL.MATRIX ../.
cp mito_score_matrix ../.
cp rf.csv ../.
cp gene.map.csv ../.
rm *_del
