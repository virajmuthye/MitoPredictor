#takes in the proteome and runs TargetP
#this requires that you have renamed the Transdecoder pep files as a four letter abbreviation
#example, kudoa iwatai was renamed as "kiwa". All proteins the transdecoder pep files will be names kiwa1, kiwa2 ....
#run as ./runTargetP.sh <four letter abbreviation> in the folder where the predicted proteome from Transdecoder are kept


for f in *Mlist.pep;
do
awk '/^>/ {OUT=substr($0,2) ".fsa"}; {print >> OUT; close(OUT)}' $f

wait

mkdir "$f".temp_targetp
mv *.fsa "$f".temp_targetp/
cd "$f".temp_targetp


for file in *; 
do
    /work/LAS/dlavrov-lab/mitoproteome/TARGETP/1.1/targetp -N $file > "$file".target

if [ ! -f "$file".target ]
then
        echo $file;
else
        rm $file; 
fi

done

cat *.target > "$f".trial.targetP.temp

grep $1 "$f".trial.targetP.temp > "$f".targetP

awk '$6=="M" {print $1}' "$f".targetP > "$f".rc12345

awk '$6=="M" && $7<=3 {print $1}' "$f".targetP > "$f".rc123

rm  *.target
rm *.temp

cp "$f".targetP ../.
cp "$f".rc12345 ../.
cp "$f".rc123 ../.
cd ..

rm -rf "$f".temp_targetp

done
echo "done"




