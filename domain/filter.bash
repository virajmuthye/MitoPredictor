# We identify domains, for which XXX % of proteins are mitochondrial. so, here, if a==1, it means it is
# found only in mitochondrial proteins, whereas a>0 means that it is present in at least 1 mitochondrial protein
# the higher the value the more stric the search
# Use as follows : ./filter.bash <cutoff>   (cut-off is between 0 and 1)
# Mitodomains are stored in mito_domains

for file in *.dom;do
	awk -F, -v a=$1 '$4>=a {print $1}' $file > "$file".filtered
done

wait

cat *.filtered | sort | uniq -c | awk '$1>2 {print $2}' > mito_domains

rm *filtered
