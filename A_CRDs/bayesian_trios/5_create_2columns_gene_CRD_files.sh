#!/bin/bash
DIR=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1000
DIR=/home/users/a/avalosma/scratch/5_CRDgene/merged_nominal_1
OUTDIR=$DIR/2_columns
mkdir -p $OUTDIR

for file in $DIR/*"mean"*.txt.gz ; do
        f=$(echo $file | rev | cut -d "/" -f1 | cut -c4- | rev)
#	zcat $file | cut -d' ' -f1,8 | sed -e 's/ /\t/g' > $OUTDIR/$f
        zcat $file | cut -d' ' -f1,8 > $OUTDIR/$f
done;  
  

  
