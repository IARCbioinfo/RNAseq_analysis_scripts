#!/usr/bin/env Rscript

#############################################################
### Script to perform transcript-level DE of RNA seq data ###
#############################################################

# get options
library("optparse")

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", help="folder with count files [default= %default]", metavar="character"),
  make_option(c("-g", "--group_file"), type="character", default=".", help="file with sample groups [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="count.txt", help="pattern for count file names [default= %default]", metavar="character"),
  make_option(c("-t", "--thres"), type="numeric", default=1, help="Threshold variance in gene expression (FPKM) [default= %default]", metavar="numeric"),
  make_option(c("-r", "--row.names"), type="character", default=NULL, help="Row names for group file (passed to read.table) [default= %default]", metavar="character"),
  make_option(c("-c", "--covar"), type="numeric", default=2, help="Column index of covariable to use for regression; other columns are treated as adjustment variables [default= %default]", metavar="character")
); 

## add something to specify class of group file columns

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load libraries
library(ballgown)
library(genefilter)
library(dplyr)

## load group file
print( opt$row.names )
if(is.null(opt$row.names)){
  row = NULL
}else row = as.numeric(opt$row.names)

groups = read.table(opt$group_file,h=T,row.names = row)
print(groups)

## load ballgown files
bg_data = ballgown(dataDir = opt$folder, samplePattern = opt$pattern, pData=groups)
print(bg_data)

if(!is.null(opt$row.names)) print( mean(rownames(groups)==rownames(sampleTable)) )

# filter low-abundance genes
bg_data_filt = subset(bg_data,"rowVars(texpr(bg_data))>1",genomesubset=TRUE)
print(bg_data_filt)

outdir    <- opt$out
dir.create(outdir, showWarnings = FALSE)
setwd(outdir)

# perform differential expression analysis
covar = as.numeric(opt$covar)
results_transcripts = stattest(bg_data_filt,feature="transcript",covariate=colnames(groups)[covar],adjustvars = colnames(groups)[-c(1,covar)], getFC=FALSE, meas="FPKM")
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_data_filt), geneIDs=ballgown::geneIDs(bg_data_filt), results_transcripts)
results_transcripts = arrange(results_transcripts,pval)

write.csv(results_transcripts, "transcript_results.csv", row.names=FALSE)
print("done saving DEG")

signif = subset(results_transcripts,results_transcripts$qval<0.05)

# plot
fpkm = texpr(bg_data,meas="FPKM")
fpkm = log2(fpkm+1)
#boxplot(fpkm,col=as.numeric(groups[,covar]),las=2,ylab='log2(FPKM+1)')
#plotTranscripts(ballgown::geneIDs(bg_data)[1729], bg_data, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

for(i in 1:nrow(signif)){
  pdf(paste("Transcript_DE_",signif$geneNames[i],"_",signif$geneIDs[i],".pdf",sep=""),h=3.5*floor(sqrt(length(unique(groups[,covar])))),w=3.5*ceiling(sqrt(length(unique(groups[,covar])))) )
  par(family="Times",las=1)
  plotMeans(as.character(signif$geneIDs[i]), bg_data, groupvar=colnames(groups)[covar], meas='FPKM', colorby='transcript',labelTranscripts = TRUE)
  dev.off()
}

print("done plotting")

# save results
write.table(rbind(names(opt),unlist(opt)),"options.txt",col.names=F,row.names=F,quote=F)
save(bg_data, file = "RNAseq_supervised_transcript.RData")
print("done saving results")

