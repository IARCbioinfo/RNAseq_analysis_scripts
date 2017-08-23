####################################################
### Script to perform Clustering of RNA seq data ###
####################################################

# get options
library("optparse")

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", help="folder with count files [default= %default]", metavar="character"),
  make_option(c("-g", "--group_file"), type="character", default=".", help="file with sample groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="count.txt", help="pattern for count file names [default= %default]", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, help="number of cores for statistical computation [default= %default]", metavar="numeric"),
  make_option(c("-q", "--FDR"), type="numeric", default=0.1, help="False Discovery Rate [default= %default]", metavar="numeric"),
  make_option(c("-m", "--IHW"), type="logical", default=FALSE, help="Use Independent Hypothesis Weighting for multiple-testing procedure [default= %default]", metavar="logical")
); 

## add something to specify class of group file columns

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### custom htseqcount read for compatibility with htseq 0.8 
DESeqDataSetFromHTSeqCount2 <- function( sampleTable, directory=".", design, ignoreRank=FALSE, ...){
  if (missing(design)) stop("design is missing")
  l <- lapply( as.character( sampleTable[,2] ), function(fn) read.table( file.path( directory, fn ), fill=T ) )
  if( ! all( sapply( l, function(a) all( a$V1 == l[[1]]$V1 ) ) ) )
    stop( "Gene IDs (first column) differ between files." )
  oldSpecialNames <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
  # either starts with two underscores
  # or is one of the old special names (htseq-count backward compatability)
  specialRows <- (substr(l[[1]]$V1,1,1) == "_") | l[[1]]$V1 %in% oldSpecialNames
  tbl <- sapply( l, function(a) a[,ncol(a)] ) # changed
  colnames(tbl) <- sampleTable[,1]
  rownames(tbl) <- l[[1]][,ncol(l[[1]])-1] # changed
  rownames(sampleTable) <- sampleTable[,1]
  tbl <- tbl[ !specialRows, , drop=FALSE ]
  object <- DESeqDataSetFromMatrix(countData=tbl,colData=sampleTable[,-(1:2),drop=FALSE],design=design,ignoreRank, ...)
  return(object)
}   

# load libraries
library(DESeq2)

# define some nice colors
prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5))

## load group file
groups = read.table(opt$group_file,h=T)

## build count table from count files
directory <- opt$folder
outdir    <- opt$out
dir.create(outdir, showWarnings = FALSE)
sampleFiles <- grep(opt$pattern,list.files(directory),value=TRUE) # find names
sampleTable <- data.frame(sampleName = paste( "sample", 1:length(sampleFiles),sep = "") ,fileName = sampleFiles, groups) # table for DESeq
head(sampleTable)
ddsHTSeq <- DESeqDataSetFromHTSeqCount2(sampleTable = sampleTable,directory = directory, design= as.formula( paste('~',colnames(groups)) ) )
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ] # filter out rows with no counts
setwd(outdir)
#dir.create("other", showWarnings = FALSE)

#dds$condition <- relevel(dds$condition, ref = "untreated") # to change ref level of factors

# perform differential expression analysis
para = FALSE
if(opt$cores>1){
  library("BiocParallel")
  register(MulticoreParam(opt$cores))
  para = TRUE
}
dds = DESeq(ddsHTSeq,parallel = para) 

require(gtools)
rna = resultsNames(dds)[-1]
per = combinations(nrow(per),2,rna)
res = vector("list",nrow(per))
if(opt$IHW){
  library("IHW")
  for(i in 1:length(res)) res[[i]] <- results(dds,parallel = para, alpha=opt$FDR,filterFun = ihw, contrast = list(per[i,1],per[i,2]) )
}else{
  for(i in 1:length(res)) res[[i]] <- results(dds,parallel = para, alpha=opt$FDR, contrast = list(per[i,1],per[i,2]) )
}

resOrdered <- lapply(res, function(x) x[order(x$padj),] )
#summary(res)

# get gene names
#resOrdered$genes <- 

# plot
pdf("LogFoldChanges_DE.pdf",h=3.5,w=3.5*length(res))
par(mfrow=c(1,length(res)),family="Times")
for(i in 1:length(res)) plotMA(res[[i]], ylim=c(-2,2),main=paste(per[i,1],"vs",per[i,2]))
dev.off()

pdf("Counts_smallestpval.pdf",h=3.5,w=3.5*length(res))
par(mfrow=c(1,length(res)),family="Times")
for(i in 1:length(res)) plotCounts(dds, gene=which.min(res[[i]]$padj), intgroup=colnames(groups)[1],main=paste(per[i,1],"vs",per[i,2]) )
dev.off()

# save top gene loadings
for(i in 1:length(resOrdered)) write.csv(resOrdered[[i]],file=paste("Genes_DE_",per[i,1],"_vs_",per[i,2],".csv",sep="") )

# save results
save(resOrdered, file = "RNAseq_supervised.RData")

