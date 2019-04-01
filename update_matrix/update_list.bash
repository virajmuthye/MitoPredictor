# use as follows
# if you want to add a new human protein to the mitochondrial list, 
# copy the final matrix file in this folder (if not here already)
# save the list of new proteins (Uniprot IDS) in the file named "newmito.list"
# run this script ./update_list.bash hsap
# the four letter abbreviations are as follows 1] human : hsap 2] mouse : mmus 3] c.elegans : cele 4] drosophila : dmel 5] yeast : scer

#extracts the renamed protein ids of the added mitochondrial proteins
while read p;do
	grep -w $p "$1".newmap | awk -F, '{print $1}' >> newmito_to_update
done < newmito.list

#adds those to existing list of mito proteins 
cat newmito_to_update all_mito_list | sort | uniq | sort > new_all_mito_list

# changes the final matrix using the new mito list
# change matrix name for final matrix if needed
# if column2 is renamed protein names
# mito should be last column

awk -F, '{print $2}' finalmatrix > matrix_del

while read p;do
	if grep "$q" new_all_mito_list
		then
    			echo "mito"    >> col_del
		else
    			echo "nonmito" >> col_del
	fi
done < matrix_del

#delete laste column (mito/nonmito)
awk -F, 'NF{NF--};1' < finalmatrix > finalmatrix2

#delete earlier matrix
rm finalmatrix

#add new mito/nonmito column to new matrix
paste -d, finalmatrix2 col_del > finalmatrix

#remove unnecessary files
rm col_del
rm matrix_del
rm all_mito_list

