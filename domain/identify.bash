# Identify protein domains in the query proteome
# If you have PFAM results of the proteome, you can comment out the next line
module use /opt/rit/modules
module load pfamscan

pfam_scan.pl -fasta query_renamed_100.pep -dir . -outfile query.pfam







