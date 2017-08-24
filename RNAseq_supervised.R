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
  make_option(c("-F", "--fpm"), type="numeric", default=1, help="Minimum gene FPM [default= %default]", metavar="numeric"),
  make_option(c("-m", "--IHW"), type="logical", default=FALSE, help="Use Independent Hypothesis Weighting for multiple-testing procedure [default= %default]", metavar="logical")
); 

## add something to specify class of group file columns

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### custom dataset load function to enable missing data in numeric variable
DESeqDataSet2 <- function (se, design, ignoreRank = TRUE) {
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    }
    else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }
  if (is.null(assayNames(se)) || assayNames(se)[1] != "counts") {
    message("renaming the first element in assays to 'counts'")
    assayNames(se)[1] <- "counts"
  }
  if (any(is.na(assay(se)))) 
    stop("NA values are not allowed in the count matrix")
  if (any(assay(se) < 0)) {
    stop("some values in assay are negative")
  }
  if (!is.integer(assay(se))) {
    if (any(round(assay(se)) != assay(se))) {
      stop("some values in assay are not integers")
    }
    message("converting counts to integer mode")
    mode(assay(se)) <- "integer"
  }
  if (all(assay(se) == 0)) {
    stop("all samples have 0 counts for all genes. check the counting script.")
  }
  if (all(rowSums(assay(se) == assay(se)[, 1]) == ncol(se))) {
    warning("all genes have equal values for all samples. will not be able to perform differential analysis")
  }
  if (any(duplicated(rownames(se)))) {
    warning(sum(duplicated(rownames(se))), " duplicate rownames were renamed by adding numbers")
    rnms <- rownames(se)
    dups <- unique(rnms[duplicated(rnms)])
    for (rn in dups) {
      idx <- which(rnms == rn)
      rnms[idx[-1]] <- paste(rnms[idx[-1]], c(seq_len(length(idx) - 
                                                        1)), sep = ".")
    }
    rownames(se) <- rnms
  }
  designVars <- all.vars(design)
  if (!all(designVars %in% names(colData(se)))) {
    stop("all variables in design formula must be columns in colData")
  }
  designVarsClass <- sapply(designVars, function(v) class(colData(se)[[v]]))
  if (any(designVarsClass == "character")) {
    warning("some variables in design formula are characters, converting to factors")
    for (v in designVars[designVarsClass == "character"]) {
      colData(se)[[v]] <- factor(colData(se)[[v]])
    }
  }
  if (length(designVars) == 1) {
    var <- colData(se)[[designVars]]
    if (all(var == var[1])) {
      stop("design has a single variable, with all samples having the same value.\n  use instead a design of '~ 1'. estimateSizeFactors, rlog and the VST can then be used")
    }
  }
  designVarsNumeric <- sapply(designVars, function(v) is.numeric(colData(se)[[v]]))
  if (any(designVarsNumeric)) {
    warnIntVars <- FALSE
    for (v in designVars[designVarsNumeric]) {
      if (all(colData(se)[[v]] == round(colData(se)[[v]]),na.rm = T)) {
        warnIntVars <- TRUE
      }
    }
    if (warnIntVars) {
      message("the design formula contains a numeric variable with integer values,\n  specifying a model with increasing fold change for higher values.\n  did you mean for this to be a factor? if so, first convert\n  this variable to a factor using the factor() function")
    }
  }
  designFactors <- designVars[designVarsClass == "factor"]
  missingLevels <- sapply(designFactors, function(v) any(table(colData(se)[[v]]) == 
                                                           0))
  if (any(missingLevels)) {
    message("factor levels were dropped which had no samples")
    for (v in designFactors[missingLevels]) {
      colData(se)[[v]] <- droplevels(colData(se)[[v]])
    }
  }
  singleLevel <- sapply(designFactors, function(v) all(colData(se)[[v]] == 
                                                         colData(se)[[v]][1]))
  if (any(singleLevel)) {
    stop("design contains one or more variables with all samples having the same value,\n  remove these variables from the design")
  }
  modelMatrix <- stats::model.matrix.default(design, data = as.data.frame(colData(se)))
  #if (!ignoreRank) {
  #  checkFullRank(modelMatrix)
  #}
  lastDV <- length(designVars)
  if (length(designVars) > 0 && designVarsClass[lastDV] == 
      "factor") {
    lastDVLvls <- levels(colData(se)[[designVars[lastDV]]])
    controlSynonyms <- c("control", "Control", "CONTROL")
    for (cSyn in controlSynonyms) {
      if (cSyn %in% lastDVLvls) {
        if (cSyn != lastDVLvls[1]) {
          message(paste0("it appears that the last variable in the design formula, '", 
                         designVars[lastDV], "',\n  has a factor level, '", 
                         cSyn, "', which is not the reference level. we recommend\n  to use factor(...,levels=...) or relevel() to set this as the reference level\n  before proceeding. for more information, please see the 'Note on factor levels'\n  in vignette('DESeq2')."))
        }
      }
    }
  }
  mcolsCols <- DataFrame(type = rep("input", ncol(colData(se))), 
                         description = rep("", ncol(colData(se))))
  mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
    mcolsCols
  }
  else if (all(names(mcols(colData(se))) == c("type", "description"))) {
    mcolsCols
  }
  else {
    cbind(mcols(colData(se)), mcolsCols)
  }
  object <- new("DESeqDataSet", se, design = design)
  mcolsRows <- DataFrame(type = rep("input", ncol(mcols(object))), 
                         description = rep("", ncol(mcols(object))))
  mcols(mcols(object)) <- if (is.null(mcols(mcols(object)))) {
    mcolsRows
  }
  else if (all(names(mcols(mcols(object))) == c("type", "description"))) {
    mcolsRows
  }
  else {
    cbind(mcols(mcols(object)), mcolsRows)
  }
  metadata(object)[["version"]] <- packageVersion("DESeq2")
  return(object)
}

### custom DESeq from matrix that uses DESeqDataSet2
DESeqDataSetFromMatrix2 <- function (countData, colData, design, tidy = FALSE, ignoreRank = FALSE, ...){
  print(dim(countData))
  print(dim(colData))
  if (tidy) {
    stopifnot(ncol(countData) > 1)
    rownms <- as.character(countData[, 1])
    countData <- countData[, -1, drop = FALSE]
    rownames(countData) <- rownms
  }
  stopifnot(ncol(countData) == nrow(colData))
  countData <- as.matrix(countData)
  if (is(colData, "data.frame")) 
    colData <- as(colData, "DataFrame")
  if (!is.null(rownames(colData)) & !is.null(colnames(countData))) {
    if (all(sort(rownames(colData)) == sort(colnames(countData)))) {
      if (!all(rownames(colData) == colnames(countData))) {
        stop(paste("rownames of the colData:\n  ", paste(rownames(colData), 
                                                         collapse = ","), "\n  are not in the same order as the colnames of the countData:\n  ", 
                   paste(colnames(countData), collapse = ",")))
      }
    }
  }
  if (is.null(rownames(colData)) & !is.null(colnames(countData))) {
    rownames(colData) <- colnames(countData)
  }
  se <- SummarizedExperiment(assays = SimpleList(counts = countData), 
                             colData = colData, ...)
  object <- DESeqDataSet2(se, design = design, ignoreRank) # change here
  return(object)
}

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
  print(dim(tbl))
  print(dim(sampleTable[,-(1:2),drop=FALSE]))
  object <- DESeqDataSetFromMatrix2(countData=tbl,colData=sampleTable[,-(1:2),drop=FALSE],design=design,ignoreRank=ignoreRank, ...)
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
ddsHTSeq  <- DESeqDataSetFromHTSeqCount2(sampleTable = sampleTable,directory = directory, design= as.formula( paste('~',paste(colnames(groups),collapse = " + ")  ) ) )
ddsHTSeq  <- ddsHTSeq[ rowSums(fpm(ddsHTSeq)>opt$fpm)>=2, ] # alternative filterinf based on fpm : keep if at least 2 samples with fpm > 1
ddsHTSeq2 <- ddsHTSeq[,apply(!is.na(groups),1,prod)>0]
if(sum( apply(!is.na(groups),1,prod)==0 )>0){
  print(paste(sum( apply(!is.na(groups),1,prod)==0 ),"samples with missing data removed:") )
  print(sampleFiles[apply(!is.na(groups),1,prod)==0])
}
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
dds = DESeq(ddsHTSeq2,parallel = para) 

require(gtools)
Vnames = colnames(groups)
res = vector("list",length(Vnames))
per = vector("list",length(Vnames))
for(i in 1:length(Vnames)){
  rna = levels(groups[,i])
  per[[i]] = combinations(length(rna),2,rna)
  res[[i]] = vector("list",nrow(per[[i]]))
}
if(opt$IHW){
  library("IHW")
  for(i in 1:length(res)){
    for(j in 1:length(res[[i]])) res[[i]][[j]] <- results(dds,parallel = para, alpha=opt$FDR,filterFun = ihw, contrast = c(Vnames[i],per[[i]][j,1],per[[i]][j,2]) )
  }
}else{
  for(i in 1:length(res)){
    for(j in 1:length(res[[i]])) res[[i]][[j]] <- results(dds,parallel = para, alpha=opt$FDR, contrast = c(Vnames[i],per[[i]][j,1],per[[i]][j,2]) )
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
  for(j in 1:length(res[[i]])) plotMA(res[[i]][[j]], ylim=c(-2,2),main=paste(per[[i]][j,1],"vs",per[[i]][j,2]))
  dev.off()
}

for(i in 1:length(res)){
  pdf(paste("Counts_smallestpval_",Vnames[i],".pdf",sep=""),h=3.5,w=3.5*length(res[[i]]))
  par(mfrow=c(1,length(res[[i]])),family="Times")
  for(j in 1:length(res[[i]])) plotCounts(dds, gene=which.min(res[[i]][[j]]$padj), intgroup=Vnames[i],main=paste(per[[i]][j,1],"vs",per[[i]][j,2]) )
  dev.off()
}

# save top gene loadings
for(i in 1:length(resOrdered)){
  for(j in 1:length(resOrdered[[i]])) write.csv(resOrdered[[i]][[j]],file=paste("Genes_DE_",Vnames[i],"_",per[[i]][j,1],"_vs_",per[[i]][j,2],".csv",sep="") )
}

# save results
save(resOrdered, file = "RNAseq_supervised.RData")

