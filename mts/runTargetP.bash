# takes in the proteome in a FASTA format and runs TargetP
# requires that you have renamed the Transdecoder pep files as a four letter abbreviation
# example, kudoa iwatai was renamed as "kiwa". All proteins the transdecoder pep files will haave names kiwa1, kiwa2 ....
# run as ./runTargetP.sh <four letter abbreviation> in the folder where the predicted proteome from Transdecoder are kept
# the files were split into smaller files and the software was run due to issues at ISU and is NOT a necessary step 

for f in *.fasta;
do
awk '/^>/ {OUT=substr($0,2) ".fsa"}; {print >> OUT; close(OUT)}' $f

wait

mkdir "$f".temp_targetp
mv *.fsa "$f".temp_targetp/
cd "$f".temp_targetp


	for file in *; 
	do

    		targetp -N $file > "$file".target_del 
    		tail -n +9 "$file".target_del | head -n -2 > "$file".target
    		rm "$file".target_del

		if [ ! -f "$file".target ]
		then
        		echo $file;
		else
        		rm $file; 
		fi

	done

cat *.target > "$f".targetP
awk '$6=="M" {print $1}' "$f".targetP > "$f".rc12345
awk '$6=="M" && $7<=3 {print $1}' "$f".targetP > "$f".rc123

rm *.target
rm *.temp

cp "$f".targetP ../.
cp "$f".rc12345 ../.
cp "$f".rc123 ../.
cd ..

rm -rf "$f".temp_targetp

done
echo "done"