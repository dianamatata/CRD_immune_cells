#!/bin/bash

OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans/significants

analysis_file=analysis_transCRD.txt

echo 'number of significant hits for transCRDs' > $analysis_file

for file in $OUTDIR/*0.01.txt ; do
	echo $(cat $file | wc -l) $(echo $file | rev | cut -d '/' -f1 | rev) >> $analysis_file
done

for file in $OUTDIR/*0.05.txt ; do
        echo $(cat $file | wc -l) $(echo $file | rev | cut -d '/' -f1 | rev) >> $analysis_file
done
