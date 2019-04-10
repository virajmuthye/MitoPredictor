#takes in the proteome and runs MitoFates 
#run as ./runMitoFates.sh in the folder where the predicted proteome from Transdecoder are kept

for f in *Mlist.pep;
do
awk '/^>/ {OUT=substr($0,2) ".fsa"}; {print >> OUT; close(OUT)}' $f
wait

mkdir "$f".temp_mitofates
mv *.fsa "$f".temp_mitofates/
cd "$f".temp_mitofates

for file in *; 
do
    perl /work/LAS/dlavrov-lab/mitoproteome/MitoFates/MitoFates.pl $file metazoa > "$file".mitofates

if [ ! -f "$file".mitofates ]
then

        echo $file;
else
        rm $file; 
fi

done

cat *.mitofates | grep -v "^Sequence" > "$f".mitofates.result

grep "Possessing" "$f".mitofates.result | awk '{print $1}' > "$f".mitofates.result.mitochondrial.list

rm  *.mitofates

cp "$f".mitofates.result ../.
cp "$f".mitofates.result.mitochondrial.list ../.

cd ..

rm -rf "$f".temp_mitofates

done



