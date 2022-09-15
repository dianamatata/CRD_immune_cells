# we want to add the start stop position of the gene in the files xxx# we want to add the start stop position of the gene in the files FDR used in 10...sh to create subfiles


# Directories ---------------------------------

DIR_gene=/home/users/a/avalosma/scratch/5_CRDgene/merged_1000
DIR_signif=/home/users/a/avalosma/scratch/12_TRIPLETS/significants
DIR_signif_out=/home/users/a/avalosma/scratch/12_TRIPLETS/significants/dist_corrected


# Debugs ---------------------------------

data_type='hist'
cell_type='neut'

# Main ---------------------------------

for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                name=${data_type}_${cell_type}
		file_FDR=$DIR_signif/FDR_0.05_${name}_1000000_CRD_gene_var_ALL.significant.txt
		file_out=$DIR_signif_out/FDR_0.05_${name}_1000000_CRD_gene_var_ALL.significant.txt
		file_gene=$DIR_gene/${name}_mean_mapping_CRD_gene_ALL.txt.gz
		touch $file_out
		
		while read -r line; do
			gene=$(echo $line | cut -d';' -f1)
			gene=$(echo $line | cut -d';' -f1)
			gene_start=$(zcat $file_gene | grep $gene | cut -d' ' -f3)
                        gene_stop=$(zcat $file_gene | grep $gene | cut -d' ' -f4)
			echo $line | awk "{\$3=\"$gene_start\" ; print}" |  awk "{\$4=\"$gene_stop\" ; print}" >> $file_out
		done < $file_FDR
	done
done

# #               line=$(cat $file_FDR | head -1)
