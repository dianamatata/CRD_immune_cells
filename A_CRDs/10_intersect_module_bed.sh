#!/bin/bash

BEDTOOLS=bedtools
DATADIR=/home/users/a/avalosma/scratch/3_CRD

$BEDTOOLS intersect -a $DATADIR/EGAD00001002670.ALLchr.module.merged.bed -b $DATADIR/EGAD00001002672.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR/70_72.ALLchr.module.intersect.bed

$BEDTOOLS intersect -a $DATADIR/EGAD00001002670.ALLchr.module.merged.bed -b $DATADIR/EGAD00001002673.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR/70_73.ALLchr.module.intersect.bed

$BEDTOOLS intersect -a $DATADIR/EGAD00001002673.ALLchr.module.merged.bed -b $DATADIR/EGAD00001002672.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR/72_73.ALLchr.module.intersect.bed

$BEDTOOLS intersect -a $DATADIR/EGAD00001002670.ALLchr.module.merged.bed -b $DATADIR/73_72.ALLchr.module.intersect.bed  | sort -V -k1,1 -k2,2n > $DATADIR/70_72_73.ALLchr.module.intersect.bed

# $BEDTOOLS intersect -a $DATADIR/EGAD00001002670.ALLchr.module.bed -b $DATADIR/EGAD00001002672.ALLchr.module.bed > $DATADIR/70_72.ALLchr.module.intersect.bed
