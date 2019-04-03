# takes in the proteome in a FASTA format and runs MitoFates 
# run as ./runMitoFates.sh in the folder where the predicted proteome from Transdecoder are kept

for file in *.fasta;
do
	MitoFates.pl $file metazoa > "$file".mitofates
done