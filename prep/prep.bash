#renames all sequences to a four letter abbreviation "qery"
#this can be changes simply by changing the "qery" to any other abbreviation
for file in *fasta;do
	awk '/^>/{print ">qery" ++i; next}{print}' < $file > renamed.pep
done

# removes redundancy by running CD-hit
# currently removes all proteins at 98% similarity, this can be changed by changing the -c option
cd-hit -i renamed.pep -o renamed.pep.98 -c 0.98 -n 5

#takes proteins above 100 amino-acids
cut -f1 -d" " renamed.pep.98 > "$file".pep2
rm renamed.pep.98
for file in *.pep2;
do
	mv $file "$file".temp
	awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' "$file".temp | awk '$2<100 {print $1}' | sort > "$file".below100.seqlengths
	grep ">" "$file".temp | sed 's/>//g' | sort > "$file".all
	comm -13 "$file".below100.seqlengths "$file".all > "$file".above100
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "$file".above100 "$file".temp > query_renamed_100.pep

	rm "$file".all
	rm "$file".temp
	rm "$file".below100.seqlengths
	rm "$file".above100

	cp query_renamed_100.pep ../orthology/.
	cp query_renamed_100.pep ../domain/.
done


#filters the sequences to get those beginning with M for mts-analysis
for file in *renamed_100.pep;do
	grep '^M' -B 1 --no-group-separator $file | grep ">" | sed 's/>//g' | sort | uniq > "$file".Mlist
	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "$file".Mlist $file > "$file".Mlist.fasta.delete
	awk '/^>/{f=!d[$1];d[$1]=1}f' "$file".Mlist.fasta.delete > query_renamed_100_Mlist.pep
	rm *.delete
	cp query_renamed_100_Mlist.pep ../mts/.
done

rm *.clstr
rm *Mlist

