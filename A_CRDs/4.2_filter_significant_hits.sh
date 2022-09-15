#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/4_CRD_residualized
OUT_aCRDs=$DATADIR/mapping_aCRDs_many_PCs
OUT_sCRDs=$DATADIR/mapping_sCRDs_many_PCs

PVALUE_THRESHOLD=0.05
SUMMARY_CIS_QTL_FILE=$DATADIR/summary_significant_cisQTLs.txt
echo >$SUMMARY_CIS_QTL_FILE 
echo "cell_type folder PC QTLs significant_QTLs" >$SUMMARY_CIS_QTL_FILE

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' 'methyl_neut' 'methyl_mono' 'methyl_tcell' ; do
	for folder in $OUT_aCRDs $OUT_sCRDs
		for PC in {0..50..2}; do
			FILE=$folder/permut_${cell_type}_$PC.txt
			SIGNIFICANT_CISQTL_FILE=$folder/permut_signif_${cell_type}_$PC.txt
			cmd="cat $FILE | tr ' ' '\t' | awk '{ if (\$19 > $PVALUE_THRESHOLD ) print \$6 \$19}' > $SIGNIFICANT_CISQTL_FILE " 
			eval $cmd
			COUNT_CIS_QTLS=$(cat $SIGNIFICANT_CISQTL_FILE | wc -l)
                        echo $cell_type $(echo $folder  | cut -d"_" -f4) $PC $(cat $FILE | wc -l) $COUNT_CIS_QTLS
			echo "$cell_type $(echo $folder  | cut -d"_" -f4) $PC $(cat $FILE | wc -l) $COUNT_CIS_QTLS" >> $SUMMARY_CIS_QTL_FILE
		done
	done
done

# $(cat $FILE | wc -l) is always the same for each cell_type, because there is a fixed number of genes tested per cell type

<<'COMMENTS'
                        COUNT_CIS_QTLS=0
			while read -r line; do
				PVALUE=$(echo "$( echo "$line "| tr ' ' '\t' | cut -f19)")
				if [[ "$PVALUE"<"$PVALUE_THRESHOLD" ]]; then
				    COUNT_CIS_QTLS=$(expr $COUNT_CIS_QTLS + 1)
				    echo $line>>$SIGNIFICANT_CISQTL_FILE
				fi
		    	done < "$FILE"
			echo $COUNT_CIS_QTLS
		    	echo -e "$PC\t$COUNT_CIS_QTLS">>$SUMMARY_CIS_QTL_FILE
COMMENTS
