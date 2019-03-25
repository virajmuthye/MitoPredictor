#comment out the following line if you already have a set of proteins from transcriptome
./translate.bash

######################
# Orthology analysis #
######################

#identify OGs from Proteinortho
./orthology.bash
#convert Proteinortho results in a parseable format
#results will be stored in mito.proteinortho.OG
./proteinortho_parser.bash mito.proteinortho


######################
# MTS analysis #
######################


