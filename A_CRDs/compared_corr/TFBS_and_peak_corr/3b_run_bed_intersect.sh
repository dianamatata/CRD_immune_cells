FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"
for chr in {1..22}
do
  for threshold in 10 50 100 150
  do
 	cmd3="source $FOLDER/TFBS_and_peak_corr/3bis_bed_intersect.sh $chr $threshold"
	echo $cmd3
	wsbatch -J intersect2 \
	--partition=shared-cpu \
	--time=12:00:00 \
	--mem=10000 \
	-o $FOLDER/log/intersect3.out \
	-e $FOLDER/log/intersect3.err \
	--wrap="$cmd3"
	 
  done
	  
done

