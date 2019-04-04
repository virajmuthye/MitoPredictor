Files in this folder:

1] .dom files
These are dom files which contain the following columns
domain id, #mitochondrial proteins with domain, #non-mitochondrial proteins with domain, dscore
dscore = [#mitochondrial proteins with domain\ #total proteins with domain]
dscore is an indicator for how mitochondria-specific the domain is; for eg, higher the value, more mitochondrial proteins have that domain than non-mitochondrial. 1: found only in mitochondrial proteins, 0: found only in non-mitochondrial protein. Value for dscore ranges between 0 and 1. 

2] .pfam files
PFAM results of all reference proteomes : hsap, mmus, cele, dmel and scer
PFAM libraries should be downloaded.

3] Scripts
identify.bash
 Uses pfam_scan.pl to identify domains in the query proteome
 output : query.pfam

filter.bash
 Identifies domains from 3/4 animal proteomes where the dscore is the cut-off specifief
 used as ./filter.bash 0.9
 the list is a file called mito_domains : mitochondria-specific domains at the cut-off specified

scan.bash
 Identifies proteins containing the mito-specific protein domains in the proteomes

make_matrix.bash
 Makes the final domain matrix and saves it as "domain_final.matrix"
 columns are "protein id, domain combination, contains mito-specific domain [Y/N]"


