# preps the file for further analysis
for file in *.fasta;do
cp $file prep/.
done
cd prep
./prep.bash
