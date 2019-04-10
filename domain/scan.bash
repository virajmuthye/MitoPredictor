for file in *.pfam;do

	while read p;do
		grep -w $p $file | awk '{print $1}' >> "$file".mito_prot_temp
	done < mito_domains

sort "$file".mito_prot_temp | uniq > "$file".mito_prot

rm "$file".mito_prot_temp 

done


