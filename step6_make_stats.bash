mkdir stats
cd concatanate

#Total number of mitochondrial proteins identified
echo "Total number of mitochondrial proteins" >> Stats.txt
wc -l query.mito.list >> Stats.txt
echo " " >> Stats.txt
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' query.mito.list query_renamed_100.pep > query_mitochondrial_proteome.fasta 

#Where to find the mitochondrial proteins
echo "You can find the mitochondrial proteins from the query file in query_mitochondrial_proteome.fasta" >> Stats.txt
echo " " >> Stats.txt

#Number of proteins with orthologs to known mitochondrial proteins
echo "Number of proteins with orthologs to known mitochondrial proteins" >> Stats.txt

cp ../orthology/mito.proteinortho.OG .

awk -F, '{print $1}' mito.proteinortho.OG | sort | uniq | sort > all_ogs

while read p;do 
	grep -w $p mito.proteinortho.OG | awk -F, '{print $1}' >> mito_ogs
done < all_mito_list

while read p;do 
	grep -w $p mito.proteinortho.OG | awk -F, '{print $1}' >> query_mito_ogs
done < query.mito.list

sort mito_ogs | uniq > mito_ogs_uniq
sort query_mito_ogs | uniq > query_mito_ogs_uniq

while read p;do
	grep -w $p mito.proteinortho.OG | grep "hsap" | awk -F, '{print $2}' >> hsap_query_temp
	grep -w $p mito.proteinortho.OG | grep "mmus" | awk -F, '{print $2}' >> mmus_query_temp
	grep -w $p mito.proteinortho.OG | grep "cele" | awk -F, '{print $2}' >> cele_query_temp
	grep -w $p mito.proteinortho.OG | grep "dmel" | awk -F, '{print $2}' >> dmel_query_temp
	grep -w $p mito.proteinortho.OG | grep "scer" | awk -F, '{print $2}' >> scer_query_temp
done < query_mito_ogs_uniq

for file in *query_temp;do
	sort $file | uniq > "$file".uniq
done

comm -12 hsap_query_temp.uniq all_mito_list > hsap_query
comm -12 mmus_query_temp.uniq all_mito_list > mmus_query
comm -12 cele_query_temp.uniq all_mito_list > cele_query
comm -12 dmel_query_temp.uniq all_mito_list > dmel_query
comm -12 scer_query_temp.uniq all_mito_list > scer_query

cp ../update_final_matrix/*newmap .

for file in *_query;
do
	sort $file | uniq > "$file".sorted
done


while read p;do
	grep -w $p human.newmap | awk -F, '{print $2}' >> hsap_query_mitoorthologs.txt
done < hsap_query

while read p;do
	grep -w $p mmus.newmap | awk -F, '{print $2}' >> mmus_query_mitoorthologs.txt
done < mmus_query

while read p;do
	grep -w $p cele.newmap | awk -F, '{print $2}' >> cele_query_mitoorthologs.txt
done < cele_query

while read p;do
	grep -w $p dmel.newmap | awk -F, '{print $2}' >> dmel_query_mitoorthologs.txt
done < dmel_query

while read p;do
	grep -w $p yeast.newmap | awk -F, '{print $2}' >> scer_query_mitoorthologs.txt
done < scer_query


grep "hsap" all_mito_list | sort  > hsap_mito
grep "mmus" all_mito_list | sort  > mmus_mito
grep "cele" all_mito_list | sort  > cele_mito
grep "dmel" all_mito_list | sort  > dmel_mito
grep "scer" all_mito_list | sort  > scer_mito


comm -23 hsap_mito hsap_query.sorted > hsap_query_no_mitoorthologs
comm -23 mmus_mito mmus_query.sorted > mmus_query_no_mitoorthologs
comm -23 cele_mito cele_query.sorted > cele_query_no_mitoorthologs
comm -23 dmel_mito dmel_query.sorted > dmel_query_no_mitoorthologs
comm -23 scer_mito scer_query.sorted > scer_query_no_mitoorthologs


while read p;do
	grep -w $p human.newmap | awk -F, '{print $2}' >> hsap_query_no_mitoorthologs.txt
done < hsap_query_no_mitoorthologs

while read p;do
	grep -w $p mmus.newmap | awk -F, '{print $2}' >> mmus_query_no_mitoorthologs.txt
done < mmus_query_no_mitoorthologs

while read p;do
	grep -w $p cele.newmap | awk -F, '{print $2}' >> cele_query_no_mitoorthologs.txt
done < cele_query_no_mitoorthologs

while read p;do
	grep -w $p dmel.newmap | awk -F, '{print $2}' >> dmel_query_no_mitoorthologs.txt
done < dmel_query_no_mitoorthologs

while read p;do
	grep -w $p yeast.newmap | awk -F, '{print $2}' >> scer_query_no_mitoorthologs.txt
done < scer_query_no_mitoorthologs

comm -12 mito_ogs_uniq query_mito_ogs_uniq > query_mito_ogs_with_atleast_1_species

while read p;do
	grep -w $p mito.proteinortho.OG | awk -F, '{print $2}' | grep "qery" >> query_mito_ortho
done < query_mito_ogs_with_atleast_1_species
sort query_mito_ortho | uniq > query_mito_ortho_uniq
comm -12 query.mito.list query_mito_ortho_uniq | sort > query_proteins_with_orthologs_to_known_mitochondrial_proteins.txt
comm -12 query.mito.list query_mito_ortho_uniq | wc -l >> Stats.txt
echo " " >> Stats.txt


#Number of proteins with MTS
echo "Number of proteins with MTS" >> Stats.txt
awk -F, '$4=="mito" && $8=="M" && $11=="Possessing_mitochondrial_presequence" {print $2}' FINAL.MATRIX | grep "qery" | sort > query_with_MTS.txt
wc -l query_with_MTS.txt >> Stats.txt
echo " " >> Stats.txt

echo "Number of proteins with MTS and no ortholog to reference mitochondrial proteins" >> Stats.txt
comm -23 query_with_MTS.txt query_proteins_with_orthologs_to_known_mitochondrial_proteins.txt | sort > query_with_MTS_and_no_orthologs_to_known_mitochondrial_proteins.txt
comm -23 query_with_MTS.txt query_proteins_with_orthologs_to_known_mitochondrial_proteins.txt | wc -l >> Stats.txt
echo " " >> Stats.txt
echo "Number of proteins with ortholog to reference mitochondrial proteins and no MTS" >> Stats.txt
comm -13 query_with_MTS.txt query_proteins_with_orthologs_to_known_mitochondrial_proteins.txt | sort > query_with_no_MTS_and_orthologs_to_known_mitochondrial_proteins.txt
comm -13 query_with_MTS.txt query_proteins_with_orthologs_to_known_mitochondrial_proteins.txt | wc -l >> Stats.txt
echo " " >> Stats.txt

#Species-specific mitochondrial proteins
echo "Number of species-specific mitochondrial proteins" >> Stats.txt
awk -F, '$4=="mito" && $5=="no_og" {print $2}' FINAL.MATRIX | grep "qery" | sort | uniq > query_species_specific.txt
wc -l query_species_specific.txt >> Stats.txt
echo " " >> Stats.txt

#Where to find them
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' query_species_specific.txt query_renamed_100.pep > query_species_specific.fasta 

#Species-specific mitochondrial protein matrix
while read p;do
	grep -w $p FINAL.MATRIX >> species_specific_matrix.txt_temp
done < query_species_specific.txt

cut -d, -f2,8- species_specific_matrix.txt_temp > species_specific_matrix.txt_temp2

awk -F, '{print $1}' species_specific_matrix.txt_temp2 > species_specific_list

while read p;do 
	grep -w $p mito_score_matrix | cut -d, -f2- >> species_specific_list_score
done < species_specific_list

paste -d, species_specific_matrix.txt_temp2 species_specific_list_score > species_specific_matrix.txt_temp3

echo "protein,tp,rc,mf_prob,mf_pred,mpp_site,domain_comb,orthoscore,mtsscore,domainscore,totalscore" > header

cat header species_specific_matrix.txt_temp3 > species_specific_matrix.txt


#mitochondrial domains
while read p;do
	grep -w $p FINAL.MATRIX | awk -F, '{print $13}' | sed 's/;/\n/g' >> query_mito_domain_list
done < query.mito.list

echo "Number of protein domains in mitochondrial proteins" >> Stats.txt
sort query_mito_domain_list | uniq | wc -l >> Stats.txt
echo " " >> Stats.txt

cp ../domain/query.pfam .


while read p;do
	grep -w $p query.pfam | awk '{print $6}' >> query_mito_domain_id
	grep -w $p query.pfam | awk '{print $7}' >> query_mito_domain_name
done < query.mito.list

sort query_mito_domain_id | uniq > query_mito_domain_id.uniq
sort query_mito_domain_name | uniq > query_mito_domain_name.uniq

while read p;do
	grep -w $p query.pfam >> query.mito.pfam
done < query.mito.list

echo "Number of mitochondrial proteins with at least 1 domain" >> Stats.txt
awk '{print $1}' query.mito.pfam | sort | uniq | wc -l >> Stats.txt
echo " " >> Stats.txt

while read p;do
	echo $p >> dom_count
	grep -w $p query.mito.pfam | awk '{print $1}' | sort | uniq | wc -l >> dom_count2
done < query_mito_domain_name.uniq
paste -d, dom_count dom_count2 > domain_count.matrix
echo " " >> Stats.txt

echo "Top 10 most abundant protein domains in mitochondrial proteins" >> Stats.txt
echo " " >> Stats.txt
echo "protein domain,number of proteins with the domain" >> Stats.txt 
sort -t, -n -r -k2 domain_count.matrix | head -10 >> Stats.txt


echo "protein,orthoscore,mtsscore,domainscore" >> query_mito_score.matrix
while read p;do 
	grep -w $p mito_score_matrix >> query_mito_score.matrix;
done < query.mito.list

mv Stats.txt ../stats/.
mv query_mitochondrial_proteome.fasta ../stats/.
mv query_proteins_with_orthologs_to_known_mitochondrial_proteins.txt ../stats/.
mv query_with_MTS.txt ../stats/.
mv query_with_MTS_and_no_orthologs_to_known_mitochondrial_proteins.txt ../stats/.
mv query_with_no_MTS_and_orthologs_to_known_mitochondrial_proteins.txt ../stats/.
mv query_species_specific.txt ../stats/.
mv query_species_specific.fasta ../stats/.
mv species_specific_matrix.txt ../stats/.
mv domain_count.matrix ../stats/.
mv query_mito_score.matrix ../stats/.
mv *_query_mitoorthologs.txt ../stats/.
mv *_query_no_mitoorthologs.txt ../stats/.

echo "Check the stat folder for results of this script"
