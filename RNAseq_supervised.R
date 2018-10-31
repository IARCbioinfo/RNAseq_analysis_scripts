#!/usr/bin/env Rscript

####################################################
### Script to perform Clustering of RNA seq data ###
####################################################

library("optparse")
# get options

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", help="folder with count files [default= %default]", metavar="character"),
  make_option(c("-g", "--group_file"), type="character", default=".", help="file with sample groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="count.txt", help="pattern for count file names [default= %default]", metavar="character"),
  make_option(c("-c", "--cores"), type="numeric", default=1, help="number of cores for statistical computation [default= %default]", metavar="numeric"),
  make_option(c("-q", "--FDR"), type="numeric", default=0.1, help="False Discovery Rate [default= %default]", metavar="numeric"),
  make_option(c("-G", "--genespans_file"), type="character", default=NULL, help="Gene spans file to compute FPKMs [default= %default]", metavar="character"),
  make_option(c("-t", "--thres"), type="numeric", default=1, help="Threshold gene expression (FPM or FPKM) [default= %default]", metavar="numeric"),
  make_option(c("-T", "--theta"), type="numeric", default=0, help="Threshold LFC to test [default= %default]", metavar="numeric"),
  make_option(c("-r", "--row.names"), type="character", default=NULL, help="Row names for group file (passed to read.table) [default= %default]", metavar="character"),
  make_option(c("-i", "--maxit"), type="numeric", default=100, help="Number of iterations for DESeq algorithm (passed to DESeq) [default= %default]", metavar="numeric"),
  make_option(c("-m", "--IHW"), type="logical", default=FALSE, help="Use Independent Hypothesis Weighting for multiple-testing procedure [default= %default]", metavar="logical")
); 



## add something to specify class of group file columns

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(DESeq2)
require(gtools)

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
  tbl2 <- l[[1]][,1:(ncol(l[[1]])-1)] # changed
  colnames(tbl) <- sampleTable[,1]
  rownames(tbl) <- l[[1]][,1] # changed
  rownames(sampleTable) <- sampleTable[,1]
  tbl <- tbl[ !specialRows, , drop=FALSE ]
  tbl2<- tbl2[ !specialRows, , drop=FALSE ]
  print(dim(tbl))
  print(dim(tbl2))
  print(dim(sampleTable[,-(1:2),drop=FALSE]))
  object <- DESeqDataSetFromMatrix(countData=tbl,colData=sampleTable[,-(1:2),drop=FALSE],design=design,ignoreRank=ignoreRank, ...)
  rowData(object) <- as.data.frame(tbl2)
  return(object)
}   


# define some nice colors
prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5))

## load group file
if(opt$row.names=="NULL"){
  opt$row.names=NULL
}else{
  opt$row.names=as.numeric(opt$row.names)
}
print( opt$row.names )
print(opt$group_file)
groups = read.table(opt$group_file,h=T,row.names = opt$row.names,sep=" ")
print(groups)

## build count table from count files
directory <- opt$folder
outdir    <- opt$out
dir.create(outdir, showWarnings = FALSE)
sampleFiles <- grep(opt$pattern,list.files(directory),value=TRUE) # find names
sampleTable <- data.frame(sampleName = gsub(".txt","", sampleFiles) ,fileName = sampleFiles, groups) # table for DESeq
head(sampleTable)

if(!is.null(opt$row.names)) print( mean(rownames(groups)==rownames(sampleTable)) )

print( as.formula( paste('~',paste(colnames(groups),collapse = " + ")  ) ) )

if(sum( apply(!is.na(groups),1,prod)==0 )>0){
  print(paste(sum( apply(!is.na(groups),1,prod)==0 ),"samples with missing data removed:") )
  print(sampleFiles[apply(!is.na(groups),1,prod)==0])
}
IDs = apply(!is.na(groups),1,prod)>0
ddsHTSeq  <- DESeqDataSetFromHTSeqCount2(sampleTable = sampleTable[IDs,],directory = directory, design= as.formula( paste('~',paste(colnames(groups),collapse = " + ")  ) ) )
groups = subset.data.frame(groups,IDs)

if(!is.null(opt$genespans_file) ){
  genespans = read.table(opt$genespans_file)
  gs2 = merge(rowData(ddsHTSeq)$V1,genespans,by.x=1, by.y=1,all.x=T,all.y=F,sort=F)
  print(head(gs2))
  gr1 <- GRanges(gs2[,2],IRanges(gs2[,3],gs2[,4]),strand = gs2[,5] )
  rnames = rownames(ddsHTSeq)
  rowRanges(ddsHTSeq) <- gr1
  rownames(ddsHTSeq) <- rnames
  fpkmdds = fpkm(ddsHTSeq)
  ddsHTSeq  <- ddsHTSeq[ rowSums(fpkmdds>opt$thres)>=2, ] # alternative filterinf based on fpm : keep if at least 2 samples with fpm > 1
}else{
  ddsHTSeq  <- ddsHTSeq[ rowSums(fpm(ddsHTSeq)>opt$thres)>=2, ] # alternative filterinf based on fpm : keep if at least 2 samples with fpm > thres
  print(ddsHTSeq)
}
#ddsHTSeq2 <- ddsHTSeq[,apply(!is.na(groups),1,prod)>0]
print(ddsHTSeq)
#ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ] # filter out rows with no counts
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
dds = DESeq(ddsHTSeq,parallel = para) #, maxit=as.numeric(opt$maxit) ) 

Vnames = colnames(groups)
res = vector("list",length(Vnames))
per = vector("list",length(Vnames))
for(i in 1:length(Vnames)){
  grtmp = groups[,i]
  if(class(grtmp)=="factor"){
    rna = levels(groups[,i])
    per[[i]] = combinations(length(rna),2,rna)
    res[[i]] = vector("list",nrow(per[[i]]))
  }else{
    rna = matrix(c(1,0),ncol=2)
    per[[i]] = rna
    res[[i]] = vector("list",1)
  }
}
if(opt$IHW){
  library("IHW")
  for(i in 1:length(res)){
    if(class(groups[,i])=="factor"){
      for(j in 1:length(res[[i]])) res[[i]][[j]] <- results(dds,parallel = para, alpha=opt$FDR,filterFun = ihw, contrast = c(Vnames[i],per[[i]][j,1],per[[i]][j,2]), lfcThreshold = opt$theta, altHypothesis = "greaterAbs" )
    }else{
      res[[i]][[1]] <- results(dds,parallel = para, alpha=opt$FDR, contrast = list(Vnames[i]) , lfcThreshold = opt$theta, altHypothesis = "greaterAbs" )
    }
  }
}else{
  for(i in 1:length(res)){
    if(class(groups[,i])=="factor"){
      for(j in 1:length(res[[i]])) res[[i]][[j]] <- results(dds,parallel = para, alpha=opt$FDR, contrast = c(Vnames[i],per[[i]][j,1],per[[i]][j,2]) , lfcThreshold = opt$theta, 
                                                            altHypothesis = "greaterAbs" )
    }else{
      res[[i]][[1]] <- results(dds,parallel = para, alpha=opt$FDR, contrast = list(Vnames[i]) , lfcThreshold = opt$theta, altHypothesis = "greaterAbs" )
    }
  }
}

resOrdered <- lapply(res, function(x) lapply(x, function(y) y[order(y$padj),] ) )
#summary(res)

# get gene names
#resOrdered$genes <- 

# plot
for(i in 1:length(res)){
  pdf(paste("LogFoldChanges_DE_",Vnames[i],".pdf",sep=""),h=3.5,w=3.5*length(res[[i]]))
  par(mfrow=c(1,length(res[[i]])),family="Times")
  for(j in 1:length(res[[i]])) plotMA(res[[i]][[j]], alpha = opt$FDR, ylim=c(-2,2),main=paste(per[[i]][j,1],"vs",per[[i]][j,2]))
  dev.off()
}

for(i in 1:length(res)){
  if(class(groups[,i])=="factor"){
      pdf(paste("Counts_smallestpval_",Vnames[i],".pdf",sep=""),h=3.5,w=3.5*length(res[[i]]))
    par(mfrow=c(1,length(res[[i]])),family="Times")
    for(j in 1:length(res[[i]])) plotCounts(dds, gene=which.min(res[[i]][[j]]$padj), intgroup=Vnames[i],main=paste(rownames(res[[i]][[j]][which.min(res[[i]][[j]]$padj),]),",",per[[i]][j,1],"vs",per[[i]][j,2]) )
    dev.off()
  }
}
print("done plotting")

# find genes where algo did not converge
DEGres = mcols(dds)
noconvl = as.character(DEGres[which( !DEGres$betaConv ),2])

# save DEG
for(i in 1:length(resOrdered)){
  for(j in 1:length(resOrdered[[i]])){
    resOrdered[[i]][[j]][which( rownames(resOrdered[[i]][[j]])%in%noconvl ),6] = "No_convergence"
    write.csv(resOrdered[[i]][[j]],file=paste("Genes_DE_",Vnames[i],"_",per[[i]][j,1],"_vs_",per[[i]][j,2],".csv",sep="") )
  }
}
print("done saving DEG")

# save results
write.table(rbind(names(opt),unlist(opt)),"options.txt",col.names=F,row.names=F,quote=F)
save(dds,res,resOrdered, file = "RNAseq_supervised.RData")

print("done saving results")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


