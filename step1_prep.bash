# preps the file for further analysis
module use /opt/rit/modules

for file in *.fasta;do
cp $file prep/.
done
cd prep
./prep.bash
