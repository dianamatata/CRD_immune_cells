DIR=/home/users/a/avalosma/scratch/10_CRD_QTLs/mixedCRDs/LD

analyse_file=$DIR/file_sizes.txt
echo 'nbr lines	file'> $analyse_file
for file in $DIR/*mean.txt.ld; do
	echo $(cat $file | wc -l) $(echo $file | rev| cut -d'/' -f1 | rev) >> $analyse_file
done
	
