module use /opt/rit/modules
module load ncbi-blast/2.4.0+
module load perl
module load python

perl /work/LAS/dlavrov-lab/mitoproteome/kma_tmp/january_nonbilaterian/earlier_analysis_RBBH_against_mitoproteomereference/portho/proteinortho_v5.16b/proteinortho5.pl -project=mito *.pep

rm *.psq
rm *.pin
rm *.phr






