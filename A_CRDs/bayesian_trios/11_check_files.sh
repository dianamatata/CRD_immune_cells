#!/bin/bash
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/not_signif/mapping_1000
DIR2=/home/users/a/avalosma/scratch/12_TRIPLETS/mapping_1000

analysis_file=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/A_CRDs/bayesian_trios/4.2_analysis_chunks.txt
echo ' ' > $analysis_file


echo '$DIR' > $analysis_file

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		cd $DIR
		name=${data_type}_${cell_type}
		nbr_chunks=$(ls | grep $name | wc -l)
		echo $name $nbr_chunks >> $analysis_file
	done
done

echo '$DIR2' >> $analysis_file

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                cd $DIR2
                name=${data_type}_${cell_type}
                nbr_chunks=$(ls | grep $name | wc -l)
                echo $name $nbr_chunks >> $analysis_file
        done
done

DIR3=/home/users/a/avalosma/scratch/12_TRIPLETS/significants
echo '$DIR3' >> $analysis_file

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                cd $DIR3
                name=${data_type}_${cell_type}
                nbr_chunks=$(ls | grep $name | wc -l)
                echo $name $nbr_chunks >> $analysis_file
        done
done
