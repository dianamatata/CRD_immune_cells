#!/bin/bash

DATADIR3=/home/users/a/avalosma/scratch/3_CRD

bedtools merge -i $DATADIR3/70_72.ALLchr.module.intersect.bed > $DATADIR3/70_72.ALLchr.module.merged.intersect.bed
bedtools merge -i $DATADIR3/70_73.ALLchr.module.intersect.bed > $DATADIR3/70_73.ALLchr.module.merged.intersect.bed
bedtools merge -i $DATADIR3/72_73.ALLchr.module.intersect.bed > $DATADIR3/72_73.ALLchr.module.merged.intersect.bed
bedtools merge -i $DATADIR3/70_72_73.ALLchr.module.intersect.bed > $DATADIR3/70_72_73.ALLchr.module.merged.intersect.bed
