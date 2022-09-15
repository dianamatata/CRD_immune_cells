#!/bin/bash
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/significants
DIR2=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_signif
mkdir -p $DIR2

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
		file=$DIR/FDR_0.05_${data_type}_${cell_type}_CRD_gene_var_ALL.significant.txt
		variant=$(cat $file | cut -d' ' -f8)
		gene=$(cat $file | cut -d' ' -f1 | cut -d';' -f1)
		CRD=$(cat $file | cut -d' ' -f1 | cut -d';' -f2-)
                triplet=$(paste <(echo "$variant") <(echo "$gene") <(echo "$CRD") --delimiters ';' )
		echo $triplet > $DIR2/${data_type}_${cell_type}_triplet.txt
	done
done

#                 triplet=$(paste <(echo "$variant") <(echo "$gene") <(echo "$CRD") --delimiters '\t' )


# for nominal 1
DIR=/home/users/a/avalosma/scratch/12_TRIPLETS/not_signif/merged_nominal1
DIR2=/home/users/a/avalosma/scratch/12_TRIPLETS/triplets_all

mkdir -p $DIR2

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                file=$DIR/${data_type}_${cell_type}_CRD_gene_var_ALL.txt.gz
                zcat $file | cut -d' ' -f8 > $DIR/${data_type}_${cell_type}_variants.txt
		zcat $file | cut -d' ' -f1 | cut -d';' -f1 > $DIR/${data_type}_${cell_type}_genes.txt
		zcat $file | cut -d' ' -f1 | cut -d';' -f2- > $DIR/${data_type}_${cell_type}_CRDs.txt
                paste $DIR/${data_type}_${cell_type}_variants.txt \
		$DIR/${data_type}_${cell_type}_genes.txt \
		$DIR/${data_type}_${cell_type}_CRDs.txt > $DIR2/${data_type}_${cell_type}_triplet.txt
        done
done


