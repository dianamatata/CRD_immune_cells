# Genetic variation in correlated regulatory region of Immunity

https://www.biorxiv.org/content/10.1101/2022.07.21.500922v1


code to replicate the analysis performed in this preprint.



## Data

Data downloaded from the European Genome-Phenome Archive

website: https://ega-archive.org/datasets

datasets used: 

EGAD00001002663 Illumina HiSeq 2000, 193 samples 

EGAD00010000850 DNA methylation profiles of monocytes, neutrophils and T cells from 525 healthy donors 

EGAD00001002675 RNA-Seq data for 205 mature neutrophil sample(s)

EGAD00001002670 ChIP-Seq data for 182 mature neutrophil sample(s). 

EGAD00001002671 RNA-Seq data for 212 CD4-positive, alpha-beta T cell sample(s). 

EGAD00001002673 ChIP-Seq data for 154 CD4-positive, alpha-beta T cell sample(s). 

EGAD00001002672 ChIP-Seq data for 172 CD14-positive, CD16-negative classical monocyte sample(s). 

EGAD00001002674 RNA-Seq data for 197 CD14-positive, CD16-negative classical monocyte sample(s).



and  Hi-C data from Javierre et al., Cell 2016

## Softwares
Publicly available code used in this study :
BCFtools v1.8 (based on HTSlib v1.8)
PLINK v1.90b5
R v3.5.1
HOMER v4.9 (webpage: http://homer.ucsd.edu/homer/ngs/peaks.html)
Clomics v1.0 (https://github.com/odelaneau/clomics)
QTLTools v1.3.1 (https://qtltools.github.io/qtltools)
GOrilla for identifying and visualizing enriched GO terms (http://cbl-gorilla.cs.technion.ac.il/)
