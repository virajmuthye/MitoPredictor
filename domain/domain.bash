# Calls the four scripts for domain analysis
# run as ./domain.bash dscore
# check README.txt for determining dscore, use 1 if not sure
./identify.bash
wait
./filter.bash $1
wait
./scan.bash
wait
./make_matrix.bash
echo "Done!"
