#This script takes the output from ProteinOrthov5.16b and maps every protein to a Orthologous Group Number
#Usage ./proteinortho_parser.bash $file
#Results will be in <Filename>.OG in format OGnumber,Sequence header e.g. OG1,avas1


#First, remove the header from the file
grep -v "#" $1 > "$1".temp2

#Clean the file for further analysis
cat -n "$1".temp2 | cut -f1,5- | sed 's/*//g' | sed 's/ //g' | sed 's/\t/,/g' > "$1".temp

#Save the OG numbers in one file, to be used later
awk -F "," '{print $1}' "$1".temp > filetrial

while read p;
do

	grep -w $p "$1".temp | awk -F "," '{$1=""; print $0}' | sed 's/ /\n/g' | sed '/^\s*$/d' | sed "s/^/OG$p,/g" >> "$1".OG 

done < filetrial

#Remove unnecessary files
rm "$1".temp
rm "$1".temp2
rm filetrial

echo "Done!"






