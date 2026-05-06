./test_3_res.sh

# Replace ./CpG_With_Strand.bed in CpG_intersect.sh with CpG positions file before running
./CpG_intersect.sh triple_0 triple_1 triple_2 triple_3 triple_4
./rename_files_3_res.sh
./generate_true_3_res.sh triple_0
./generate_true_3_res.sh triple_1
./generate_true_3_res.sh triple_2
./generate_true_3_res.sh triple_3
./generate_true_3_res.sh triple_4
./test_results_3_res.sh

files=`ls triple_*/intersect_*[TP,FP,FN].bed`
for file in $files; do 
echo $file
./test_results_CpG.sh $file
done
