grep "hsap" mito.proteinortho.OG | sort | uniq | sort > human_og_del
awk -F, '{print $1}' human_og_del > 1_del
awk -F, '{print $2}' human_og_del > 2_del

while read p;do
	grep -w $p human.map | awk -F, '{print $2}' >> human_geneid_del
done < 2_del
wait
paste -d, 1_del human_geneid_del > og_geneid_temp
awk -F, '{print $1}' og_geneid_temp | sort | uniq | sort > 5_del
while read p;do
	grep -w $p og_geneid_temp | head -1 >> og_geneid
done < 5_del
awk -F, '{print $1}' mito.proteinortho.OG > 3_del
while read p;do
if grep -w $p og_geneid
	then 
		grep -w $p og_geneid | awk -F, '{print $2}' >> 4_del
	else
		echo "N" >> 4_del
	fi
done < 3_del
paste -d, mito.proteinortho.OG 4_del > orthology.matrix_del
echo "og,pid,geneid" > header_del
cat header_del orthology.matrix_del > orthology.matrix
rm *_del
rm og_geneid
rm og_geneid_temp





