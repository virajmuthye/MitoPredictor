#Proteinortho and perl should be installed
#Make sure all complete proteomes are in a folder called test
#add "-single" option if you want to report proteins WITHOUT orthologs

perl -e=0.00001 -sim=0.95 -identity=25 -cov=50 -purity 1 -conn 0.1 proteinortho5.pl -project=mito $file *pep






