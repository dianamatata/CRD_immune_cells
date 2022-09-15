#!/bin/bash

FOLDER=/home/users/a/avalosma/scratch/2_CRD/
mkdir -p $FOLDER/CRD_names

for file in $FOLDER/quantify_ALL/*gz ; do
	f="$(basename -- $file)"
	zcat $file | cut -f4 > $FOLDER/CRD_names/${f: : -3}
done
