# RNAseq_analysis_scripts
Scripts for RNA seq analysis. They can be used directly using the outputs from [*IARCbioinfo's RNA seq workflow*].(https://github.com/IARCbioinfo/RNAseq-nf)

## Unsupervised analysis: *RNAseq_unsupervised.R*

This script performs unsupervised analyses (Principal Component Analysis and clustering) from htseq-count outputs. 

## Prerequisites
This R script requires the following packages:
- ConsensusClusterPlus
- ade4
- AnnotationDbi
- org.Hs.eg.db
- DESeq2
- biomaRt
- fpc


### Usage
```bash
Rscript cluster.R [options]
```

| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*-f* | . | folder with count files |
*-o* | out | output directory name |
*-p* |  count.txt | pattern for count file names |
*-n* | 500 | number of genes to use for clustering |
*-t* | auto | count transformation method; 'rld', 'vst', or 'auto' |
*-c* | hc | clustering algorithm to be passed to ConsensusClusterPlus|
*-l*  | complete | method for hierarchical clustering to be passed to ConsensusClusterPlus|
*-h*    |  | Show help message and exit|

For example, one can type
```bash
Rscript cluster.R -f input -p count -t rld -o output/ -n 500
```

### Details
The script involves 3 steps
- **Data transformation** using either variance-stabilization log2-tranform or the r-log tranform from package *DESeq2*
- **PCA** of tranformed counts with package *ade4*
- **Clustering** of transformed counts with package *ConsensusClusterPlus*

### Output
- a .RData file with 4 objects: transformed counts *di*, the *ade4* pca object *pca*,the ConsensusClusterPlus object *clusters*, and the list of cluster and item consensus *icl*

#### PCA
- PCA plots with the first 2 PC, with colors corresponding to the clusters from ConsensusClusterPlus, for K=1,2,...,6 clusters.
- loadings plot for the first 2 PC
- csv file with the genes with the greatest loadings in the first PC
#### Clustering
- a folder with plots from ConsensusClusterPlus
