# Identify protein domains in the query proteome
# If you have PFAM results of the proteome, you can comment out the next line
pfam_scan.pl -fasta query_renamed_100.pep -dir . -outfile query.pfam

wait

# Identify mitochondrial proteins with mito_only domains
# Mito_only domains were identified by filter.bash
while read p;do
	grep -w $p query.pfam | awk '{print $1}' >> query_mito_prot_temp
done < mito_domains
sort query_mito_prot_temp | uniq > query_mito_prot
rm query_mito_prot_temp







