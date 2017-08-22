# RNAseq_analysis_scripts
Scripts for RNA seq analysis. They can be used directly using the outputs from [*IARCbioinfo's RNA seq workflow*].(https://github.com/IARCbioinfo/RNAseq-nf)

## Unsupervised analysis: *RNAseq_unsupervised.R*

This script performs unsupervised analyses (Principal Component Analysis and clustering) from htseq-count outputs. 

## Prerequisites
This R script requires the following packages:
- ConsensusClusterPlus
- ade4
- DESeq2
- fpc


### Usage
```bash
Rscript RNAseq_unsupervised.R [options]
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
Rscript RNAseq_unsupervised.R -f input -p count -t rld -o output/ -n 500
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

## Compare an unsupervised analysis with a list of variables: *RNAseq_unsupervised_compare.R*

This script compares the result of an unsupervised analyses (Principal Component Analysis and clustering) obtained for example using script *RNAseq_unsupervised_compare.R* with an arbitrary number of variables (categorical or continuous).

## Prerequisites
This R script requires the following packages:
- ade4
- permute


### Usage
```bash
Rscript RNAseq_unsupervised_compare.R [options]
```

| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*-R* | . | .RData file with results from clustering in variable *clusters* and results from PCA in variable *pca* |
*-i* | . | name of input file with variables in column and variable names as first line |
*-o* | out | output file preffix |
*-h*    |  | Show help message and exit|

For example, one can type
```bash
Rscript RNAseq_unsupervised_compare.R -R RNAseq_unsupervise.RData -i variables.txt -o output/
```

### Details
For each clustering present in variable *cluster*, the script involves 3 steps
- **Plotting** the variables on the first two PC axes
- **Computing the best matching** between clusters and the levels of the variable *(note: if the variable is continuous, categories are formed by subdividing the range of the variable into intervals of the same size)*
- **Testing the significance of the best matching** by computing a null distribution of matchings across 1000 permutations

### Output
For each column (i.e., variable) of the input table, a .pdf file with *K* rows, where *K* is the number of clusterings in variable *cluster*, and 3 columns: 
- Column 1 represents the first 2 PCs of the PCA; colored ellipses correspond to the clusters of the *clusters* variable from the .RData file, and point colors correspond to the levels of the variable
- Column 2 represents the best matching between the clustering and the variable
- Column 3 represents the null distribution of the best matching (gray), the observed best matching and the p-value of the best matching (red)

## Supervised analysis: *RNAseq_supervised.R*

This script performs supervised analyses (Differential Expression Analysis) from htseq-count outputs. 

## Prerequisites
This R script requires the following package:
- DESeq2

Depending on the options used, the following packages are also required:
- BiocParallel (with option -c)
- IHW (with option -m)

### Usage
```bash
Rscript RNAseq_supervised.R [options]
```

| **PARAMETER** | **DEFAULT** | **DESCRIPTION** |
|-----------|--------------:|-------------| 
*-f* | . | folder with count files |
*-g* | . | file with sample groups |
*-o* | out | output directory name |
*-p* |  count.txt | pattern for count file names |
*-c*   |  1 | number of cores for statistical computation |
*-q* | 0.1 | False Discovery Rate |
*-m* | FALSE | Use Independent Hypothesis Weighting for multiple-testing procedure |
*-h*    |  | Show help message and exit|



For example, one can type
```bash
Rscript RNAseq_supervised.R -f input -g groups.txt -o output/
```

### Details
The script performs DE analysis of gene count data under a Poisson glm with package *DESeq2*

### Output
- a plot of fold changes as a function of normalized counts
- a plot of normalized counts of the most significant DE gene, as a function of the group
- a .csv file with gene names, fold changes, and p-values
- a .RData file with 1 object: the deseq2 results table
