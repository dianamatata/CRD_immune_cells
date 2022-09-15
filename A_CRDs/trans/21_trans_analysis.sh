#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
CLOMICs=/home/users/a/avalosma/bin/clomics_nikos/bin/clomics
DATACRD=/home/users/a/avalosma/scratch/2_CRD
OUTDIR=/home/users/a/avalosma/scratch/7_CRD_Trans
mkdir -p $OUTDIR $OUTDIR/OUTsimple


for data_type in  'methyl' 'hist' ; do
        for cell_type in 'neut' 'mono' 'tcell' ; do
                for module in 'mean' 'loom' ; do                        
			MOD=$DATACRD/quantify_ALL/${data_type}_${cell_type}.ALLchr.${module}.txt.gz
			for c0 in $(seq 1 22); do
				for c1 in $(seq $(($c0 + 1)) 22); do
					OUT=$OUTDIR/OUTsimple/${data_type}_${cell_type}_${module}_chr_${c0}_${c1}_trans
					cmd="$CLOMICs trans --bed $MOD --regions $c0 $c1 --out ${OUT}.txt.gz --full"
					echo $cmd
					eval $cmd
				done
			done
		done
	done
done


# expected 2772 files

# verification 
cd $OUTDIR
count=$(ls | cut -d '_' -f1-3 | uniq -c)
count2=0
for c0 in $(seq 1 22); do
	for c1 in $(seq $(($c0 + 1)) 22); do
		count2=$((count2+1))	
	done
done

echo "$count" 
echo "expected: $count2" # somme de 1 a n= n*(n+1)/2 with n=21, result=231
