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
  make_option(c("-t", "--transform"), type="character", default="auto", help="count transformation method; 'rld', 'vst', or 'auto' [default= %default]", metavar="character"),
  make_option(c("-c", "--clusteralg"), type="character", default="hc", help="clustering algorithm to be passed to ConsensusClusterPlus; 'km' (k-means on data matrix), 'kmdist' (k-means on distances), 'hc' (hierarchical clustering), 'pam' (paritioning around medoids) [default= %default]", metavar="character"),
  make_option(c("-l", "--linkage"), type="character", default="complete", help="method for hierarchical clustering to be passed to ConsensusClusterPlus; 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', or 'centroid' [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load libraries
require(ConsensusClusterPlus)
require(xlsx)
library(ade4)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(biomaRt)
library(fpc)

# define some nice colors
prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5))

## build count table from count files
directory <- opt$folder
n = opt$nbgenes
print(n)
outdir    <- opt$out
dir.create(outdir, showWarnings = FALSE)
sampleFiles  <- grep(opt$pattern,list.files(directory),value=TRUE) # find names
sampleTable <- data.frame(sampleName = paste( "sample", 1:length(sampleFiles),sep = "") ,fileName = sampleFiles) # table for DESeq
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory, design= ~ 1)
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
    di  <- vst(ddsHTSeq, blind = TRUE)  #vst transformation
    tran = "vst"
  }else{
    rld <- rlog(ddsHTSeq, blind = TRUE) # rlog transformation
    vsd <- vst(ddsHTSeq, blind = TRUE)  #vst transformation
    
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
genes = data.frame(ensembl_gene_id=rownames(pca$c1[c1>1/length(c1),])[idgc1$ix],stringsAsFactors = F)
genes$gene_loading_PC1 = pca$c1[,1][c1>1/length(c1)][idgc1$ix]
genes$gene_loading_PC2 = pca$c1[,2][c1>1/length(c1)][idgc1$ix]

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
shortnames = sapply(genes$ensembl_gene_id, function(x) strsplit(x,split = "\\.")[[1]][1] )
genes$ensembl_symbol <- mapIds(org.Hs.eg.db, keys=shortnames,column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# save top gene loadings
write.csv(genes,file="PCA/Genes_PC1.csv")

# clustering
# select genes with largest variation across samples
d = assay(di)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:n],] 
d = sweep(d,1, apply(d,1,median,na.rm=T))

# compute clusters and consensus
clusters = ConsensusClusterPlus(d,maxK = 6,title=paste("clustering/ConsensusClusterPlus",opt$clusteralg,opt$linkage,sep="_"),innerLinkage=opt$linkage,clusterAlg=opt$clusteralg,distance="euclidean" ,seed=1,plot="png")
icl     = calcICL(clusters,title=paste("clustering/ConsensusClusterPlus",opt$clusteralg,opt$linkage,sep="_"),plot="png")

# compute clustering stats
stcl   = lapply(2:6, function(i) cluster.stats(dist(t(d)),clusters[[i]]$consensusClass) )
ldunn = sapply(1:5, function(i) stcl[[i]]$dunn )
lwbr  = sapply(1:5, function(i) stcl[[i]]$wb.ratio ) #c(stcl.B$wb.ratio,stcl.hc$wb.ratio,stcl.Whc$wb.ratio,stcl.W2hc$wb.ratio,stcl.km$wb.ratio)
lch   = sapply(1:5, function(i) stcl[[i]]$ch ) #c(stcl.B$ch,stcl.hc$ch,stcl.Whc$ch,stcl.W2hc$ch,stcl.km$ch)

# plot clustering stats
pdf("clustering/Cluster_separation_stats.pdf",h=3.5,w=3.5*3)
par(mfrow=c(1,3),family="Times")
co = rep(1,5)
co[which.max(ldunn)]=2
par(mfrow=c(1,3),family="Times")
barplot(ldunn ,names.arg = paste("K =",2:6) ,las=2,ylab="Dunn index",col= co)
co = rep(1,5)
co[which.min(lwbr)]=2
barplot(lwbr ,names.arg = paste("K =",2:6) ,las=2,ylab="Within-between SS ratio",col= co)
co = rep(1,5)
co[which.max(lch)]=2
barplot(lch ,names.arg = paste("K =",2:6) ,las=2,ylab="Calinski-Harabasz index",col= co)
dev.off()

# plot PCA with clusters
pdf("PCA/PCA.pdf",h=3,w=3*5)
par(mfrow=c(1,5),family="Times")
s.class(pca$li,as.factor(clusters[[2]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(clusters[[3]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(clusters[[4]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(clusters[[5]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(clusters[[6]]$consensusClass),col=prettycolors,xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
dev.off()

# save results
save(di,pca,clusters,icl, file = "RNAseq_unsupervised.RData")
