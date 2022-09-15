FOLDER="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/compared_corr"

cmd="source $FOLDER/merge_sub_files.sh"

wsbatch -J merge.job \
--partition=shared-bigmem  \
--time=12:00:00 \
--mem=100000  \
-o $FOLDER/log/merge.out  \
-e $FOLDER/log/merge.err  \
--wrap="$cmd"

