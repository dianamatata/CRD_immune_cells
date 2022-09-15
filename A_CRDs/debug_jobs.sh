#!/bin/bash

# loop across folder, jobs with .err
# check if canceled due to time limit or space

directory=/home/users/a/avalosma/scratch/2_CRD/OUT
error_message='canceled'
error_message='home'
#  Failed to open file
for file in $directory/*.err; do
	egrep -w 'warning|error|critical|Failed|canceled' $file
done

	if grep -Rq $error_message $file
	then
		echo $file
	fi
done

# grep list of error words: gfep -wFm1
egrep -w 'warning|error|critical|Failed|canceled' $file
        echo $file
        cat $file | grep $error_message
