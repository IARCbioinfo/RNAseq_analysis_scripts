#!/usr/bin/env Rscript

####################################################
### Script to perform Clustering of RNA seq data ###
####################################################

# get options
library("optparse")

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=".", help="folder with count files [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="count.txt", help="pattern for count file names [default= %default]", metavar="character"),
  make_option(c("-n", "--nbgenes"), type="numeric", default=500, help="number of genes to use for clustering [default= %default]", metavar="number"),
  make_option(c("-s", "--nsub"), type="numeric", default=1000, help="the number of genes to subsample for vst [default= %default]", metavar="number"),
  make_option(c("-t", "--transform"), type="character", default="auto", help="count transformation method; 'rld', 'vst', or 'auto' [default= %default]", metavar="character"),
  make_option(c("-c", "--clusteralg"), type="character", default="hc", help="clustering algorithm to be passed to ConsensusClusterPlus; 'km' (k-means on data matrix), 'kmdist' (k-means on distances), 'hc' (hierarchical clustering), 'pam' (paritioning around medoids) [default= %default]", metavar="character"),
  make_option(c("-l", "--linkage"), type="character", default="complete", help="method for hierarchical clustering to be passed to ConsensusClusterPlus; 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid' [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load libraries
library(ConsensusClusterPlus)
library(ade4)
library(DESeq2)
library(fpc)
library(cluster)

# define some nice colors
prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5),rgb(0,124/255,146/255),rgb(178/255,111/255,22/255),rgb(234/255,171/255,0),rgb(83/255,40/255,79/255))

### custom htseqcount read for compatibility with htseq 0.8 
DESeqDataSetFromHTSeqCount2 <- function( sampleTable, directory=".", design, ignoreRank=FALSE, ...) 
{
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

## build count table from count files
directory <- opt$folder
n = opt$nbgenes
outdir    <- opt$out
dir.create(outdir, showWarnings = FALSE)
sampleFiles  <- grep(opt$pattern,list.files(directory),value=TRUE) # find names
sampleTable <- data.frame(sampleName = paste( "sample", 1:length(sampleFiles),sep = "") ,fileName = sampleFiles) # table for DESeq
ddsHTSeq <- DESeqDataSetFromHTSeqCount2(sampleTable = sampleTable,directory = directory, design= ~ 1)
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ] # filter out rows with no counts
setwd(outdir)
dir.create("PCA", showWarnings = FALSE)
dir.create("clustering", showWarnings = FALSE)

# compute tranformation
if(opt$transform=="rld"){
  di <- rlog(ddsHTSeq, blind = TRUE) # rlog transformation
  tran = "r-log"
}else{
  if(opt$transform=="vst"){
    di  <- vst(ddsHTSeq, blind = TRUE, nsub = opt$nsub)  #vst transformation
    tran = "vst"
  }else{
    rld <- rlog(ddsHTSeq, blind = TRUE) # rlog transformation
    vsd <- vst(ddsHTSeq, blind = TRUE, nsub = opt$nsub)  #vst transformation
    
  ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
  df <- bind_rows( as_data_frame(log2(counts(ddsHTSeq, normalized=TRUE)[, 1:2]+1)) %>% mutate(transformation = "log2(x + 1)"),
                 as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
                 as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")) #plot 2 first samples
  colnames(df)[1:2] <- c("x", "y")  
  pdf("Compare_transformations.pdf",h=4,w=4*3)
  par(mfrow=c(1,1),family="Times")
  print(ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation) )# find best dist
  dev.off()
  if( length(sampleFiles)<30 ){
    di = rld
    tran = "r-log"
  }else{
    di = vsd
    tran = "vst"
  }
  }
}

# plot tranformation effect
l2t = log(assay(ddsHTSeq)+1,2)
qdi = apply(assay(di),2,quantile,c(0.025,0.5,0.975))
mal = max(l2t)
mav = max(assay(di))
ma = max(c(mal,mav))
ql2t = apply(l2t,2,quantile,c(0.025,0.5,0.975))
png("Transformation_effect.png",h=400,w=400*5)
par(mfrow=c(1,2),family="Times",las=2)
boxplot( log(counts(ddsHTSeq,normalized=F)+1,2),ylim=c(0,ma) ,xlab="Samples",ylab = "log2-transformed counts")
lines(ql2t[1,],col=2)
lines(ql2t[2,],col=2)
lines(ql2t[3,],col=2)
boxplot(assay(di),ylim=c(0,ma)  ,xlab="Samples", ylab= paste(tran,"counts") )
lines(qdi[1,],col=2)
lines(qdi[2,],col=2)
lines(qdi[3,],col=2)
dev.off()

# PCA 
pca <- dudi.pca(t(assay(di)),scannf = F,nf = 10,center = T, scale = F)

# genes that contribute the most to variance in 1st PC
c1 = (pca$c1[,1])**2
idgc1   = sort.int(c1[c1>1/length(c1)],decreasing = T,index.return = T)
genes = data.frame(gene_name=rownames(pca$c1)[c1>1/length(c1)][idgc1$ix],stringsAsFactors = F)
genes$gene_loading_PC1 = pca$c1[,1][c1>1/length(c1)][idgc1$ix]
if(ncol(pca$c1)>1) genes$gene_loading_PC2 = pca$c1[,2][c1>1/length(c1)][idgc1$ix]

# plot PCA loadings
pdf("PCA/PCA_loadings.pdf",h=4.5,w=4.5)
par(mfrow=c(1,1),family="Times",las=1)
cc1 = sort.int(c1,decreasing = T,index.return = T)
plot(cumsum( cc1$x) ,type="l",xlab="Gene",ylab="Cumulative distribution of loadings",lwd=2)
lines(cumsum((pca$c1[cc1$ix,2])**2) ,col=2,lwd=2)
polygon( c(0,length( idgc1$ix ),length( idgc1$ix ),0), c(0,0,rep(sum(c1[c1>1/length(c1)][idgc1$ix]),2))  ,border = NA,col=rgb(1,0,0,0.5))
legend("bottomright",legend = c("PC1","PC2"),col=1:2,lwd=2)
text(length( idgc1$ix )+8000,0.05,"loading > 1/n",col=rgb(1,0,0,0.5),offset=1)
dev.off()

# get gene names
#shortnames = sapply(genes$ensembl_gene_id, function(x) strsplit(x,split = "\\.")[[1]][1] )
#genes$gene_name <-  #mapIds(org.Hs.eg.db, keys=shortnames,column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# save top gene loadings
write.csv(genes,file="PCA/Genes_PC1.csv")

# clustering
# select genes with largest variation across samples
d = assay(di)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:n],] 
d = sweep(d,1, apply(d,1,median,na.rm=T))

# compute clusters and consensus
maxK = min(10,ncol(d)-2)
clusters = ConsensusClusterPlus(d,maxK = maxK,title=paste("clustering/ConsensusClusterPlus",opt$clusteralg,opt$linkage,sep="_"),innerLinkage=opt$linkage,clusterAlg=opt$clusteralg,distance="euclidean" ,seed=1,plot="png")
icl     = calcICL(clusters,title=paste("clustering/ConsensusClusterPlus",opt$clusteralg,opt$linkage,sep="_"),plot="png")

# compute clustering stats
stcl   = lapply(2:maxK, function(i) cluster.stats(dist(t(d)),clusters[[i]]$consensusClass) )
ldunn = sapply(1:(maxK-1), function(i) stcl[[i]]$dunn )
lwbr  = sapply(1:(maxK-1), function(i) stcl[[i]]$wb.ratio ) #c(stcl.B$wb.ratio,stcl.hc$wb.ratio,stcl.Whc$wb.ratio,stcl.W2hc$wb.ratio,stcl.km$wb.ratio)
lch   = sapply(1:(maxK-1), function(i) stcl[[i]]$ch ) #c(stcl.B$ch,stcl.hc$ch,stcl.Whc$ch,stcl.W2hc$ch,stcl.km$ch)
lsil = vector("list",(maxK-1))
pdf(paste("clustering/Silhouette_",opt$clusteralg,"_",opt$linkage,".pdf",sep=""),h=4,w=4*(maxK-1))
par(mfrow=c(1,(maxK-1)))
for(i in 2:maxK){
  sil = silhouette(clusters[[i]]$consensusClass,dist(t(d),method = "euclidean"))
  sizes = table(clusters[[i]]$consensusClass)
  plot( sil ,col=rep( clusters[[i]]$clrs[[3]],rep=sizes) ,main=paste("K=",i))
  lsil[[i-1]]=sil
}
dev.off()

msil = sapply(1:(maxK-1), function(i) mean( lsil[[i]][,3] ) )
cdl = lapply(2:maxK, function(i) as.dist(1-clusters[[i]]$consensusMatrix ) )
md = dist( t(d),method = "euclidean")
corl =sapply(cdl, cor,md)

# plot clustering stats
pdf(paste("clustering/Cluster_separation_stats_",opt$clusteralg,"_",opt$linkage,".pdf",sep=""),h=3.5,w=3.5*(maxK-1))
par(mfrow=c(1,(maxK-1)),family="Times")
co = rep(1,(maxK-1))
co[which.max(ldunn)]=2
barplot(ldunn ,names.arg = paste("K =",2:maxK) ,las=2,ylab="Dunn index",col= co)
co = rep(1,(maxK-1))
co[which.min(lwbr)]=2
barplot(lwbr ,names.arg = paste("K =",2:maxK) ,las=2,ylab="Within-between SS ratio",col= co)
co = rep(1,(maxK-1))
co[which.max(lch)]=2
barplot(lch ,names.arg = paste("K =",2:maxK) ,las=2,ylab="Calinski-Harabasz index",col= co)
co = rep(1,(maxK-1))
co[which.max(msil)]=2
barplot(msil,names.arg = paste("K =",2:maxK) ,las=2,ylab="Mean Silhouette",col= co)
co = rep(1,(maxK-1))
co[which.max(corl)]=2
barplot(corl,names.arg = paste("K =",2:maxK) ,las=2,ylab="Mean cophenetic distance",col= co)
dev.off()

# plot PCA with clusters
pdf(paste("PCA/PCA_",opt$clusteralg,"_",opt$linkage,".pdf",sep=""),h=3,w=3*(maxK-1))
par(mfrow=c(1,(maxK-1)),family="Times")
for(i in 2:(maxK)) s.class(pca$li,as.factor(clusters[[i]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
dev.off()

# save results
save(di,pca,clusters,icl, file = paste("RNAseq_unsupervised_",opt$clusteralg,"_",opt$linkage,".RData",sep="") )
