#rename transcriptome
#requires Transdecoder to be installed
#file should be in FASTA format, and have extension ".fasta"
#we keep both maps for future use:
#map1- original transcript name, tspt
#map2- map between transcript and predicted protein 
for file in *.fasta;do
	awk '/^>/{print ">tspt" ++i; next}{print}' < $file > "$file".renamed.pep
	grep ">" $file   | sed 's/>//g' > del1
	grep ">" "$file".renamed.pep | sed 's/>//g' > del2
	paste -d, del1 del2 > transcriptome.MAP
	rm del1
	rm del2
done

#translates proteins using transdecoder
for file in *.renamed.pep;do
	TransDecoder.LongOrfs -t $file
	TransDecoder.Predict -t $file
done

#rename predicted protein set
for file in *.transdecoder.pep;do
	awk '/^>/{print ">tspt" ++i; next}{print}' < $file > query.fasta
	grep ">" $file   	   | sed 's/>//g' > del1
	grep ">" query.fasta   | sed 's/>//g' > del2
	paste -d, del1 del2 > transcript_protein.MAP
	rm del1
	rm del2
done
#the predicted protein set will be in query.fasta
