####################################################
### Script to perform Clustering of RNA seq data ###
####################################################

require("ConsensusClusterPlus")
require(xlsx)
library("dplyr")
library("ggplot2")
library("ade4")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("DESeq2")
library("biomaRt")

args = commandArgs(trailingOnly=TRUE)
#1: folder with count files
#2: pattern for count files
#3: count transformation method; rld, vst, or automatic
#4: out directory
#5: number of genes to use for clustering

# define some nice colors
prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5))

## build count table from count files
directory <- args[1]
n = as.numeric(args[5])
print(n)
outdir    <- args[4]
dir.create(outdir, showWarnings = FALSE)
sampleFiles  <- grep(args[2],list.files(directory),value=TRUE) # find names
sampleTable <- data.frame(sampleName = paste( "sample", 1:length(sampleFiles),sep = "") ,fileName = sampleFiles) # table for DESeq
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory, design= ~ 1)
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ] # filter out rows with no counts

setwd(outdir)

## compute tranformation
if(args[3]=="rld"){
  di <- rlog(ddsHTSeq, blind = TRUE) # rlog transformation
}else{
  if(args[3]=="vst"){
    di  <- vst(ddsHTSeq, blind = TRUE)  #vst transformation
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
  if( length(sampleFiles)<30 ) di = rld
  else di = vsd
  }
}

### PCA 
pca <- dudi.pca(t(assay(di)),scannf = F,nf = 10,center = T, scale = F)

# genes that contribute the most to variance in 1st PC
c1 = (pca$c1[,1])**2
idgc1   = sort.int(c1[c1>1/length(c1)],decreasing = T,index.return = T)
genes = data.frame(ensembl_gene_id=rownames(pca$c1[c1>1/length(c1),])[idgc1$ix],stringsAsFactors = F)
genes$gene_loading_PC1 = pca$c1[,1][c1>1/length(c1)][idgc1$ix]
genes$gene_loading_PC2 = pca$c1[,2][c1>1/length(c1)][idgc1$ix]

pdf("PCA_loadings.pdf",h=4.5,w=4.5)
par(mfrow=c(1,1),family="Times",las=1)
cc1 = sort.int(c1,decreasing = T,index.return = T)
plot(cumsum( cc1$x) ,type="l",xlab="Gene",ylab="Cumulative distribution of loadings",lwd=2)
lines(cumsum((pca$c1[cc1$ix,2])**2) ,col=2,lwd=2)
polygon( c(0,length( idgc1$ix ),length( idgc1$ix ),0), c(0,0,rep(sum(c1[c1>1/length(c1)][idgc1$ix]),2))  ,border = NA,col=rgb(1,0,0,0.5))
legend("bottomright",legend = c("PC1","PC2"),col=1:2,lwd=2)
text(length( idgc1$ix )+8000,0.05,"loading > 1/n",col=rgb(1,0,0,0.5),offset=1)
dev.off()

shortnames = sapply(genes$ensembl_gene_id, function(x) strsplit(x,split = "\\.")[[1]][1] )
genes$ensembl_symbol <- mapIds(org.Hs.eg.db, keys=shortnames,column="SYMBOL", keytype="ENSEMBL", multiVals="first")

write.csv(genes,file="Genes_PC1.csv")

## clustering
if(is.na(n)) n = sum(c1>1/length(c1))
d = assay(di)[(c1>1/length(c1))[1:n],]
#mads=apply(d,1,mad)
#d=d[rev(order(mads))[1:400],] #use 400 most variable genes 
d = sweep(d,1, apply(d,1,median,na.rm=T))

resultshc = ConsensusClusterPlus(d,maxK = 6,title="ConsensusClusterPlus_hc",innerLinkage="complete",clusterAlg="hc",distance="euclidean" ,seed=1,plot="png")
resultskm = ConsensusClusterPlus(d,maxK = 6,title="ConsensusClusterPlus_km",innerLinkage="complete",clusterAlg="kmdist",distance="euclidean" ,seed=1,plot="png")
iclhc     = calcICL(resultshc,title="ConsensusClusterPlus_hc",plot="png")
iclkm     = calcICL(resultskm,title="ConsensusClusterPlus_km",plot="png")

save(resultshc,resultskm,iclhc,iclkm, file = "clusters.RData")

## plot pca with clusters
pdf("PCA_hc.pdf",h=3,w=3*5)
par(mfrow=c(1,5),family="Times")
s.class(pca$li,as.factor(resultshc[[2]]$consensusClass),col=prettycolors[as.numeric(resultshc[[2]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultshc[[3]]$consensusClass),col=prettycolors[as.numeric(resultshc[[3]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultshc[[4]]$consensusClass),col=prettycolors[as.numeric(resultshc[[4]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultshc[[5]]$consensusClass),col=prettycolors[as.numeric(resultshc[[5]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultshc[[6]]$consensusClass),col=prettycolors[as.numeric(resultshc[[6]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
dev.off()

pdf("PCA_km.pdf",h=3,w=3*5)
par(mfrow=c(1,5),family="Times")
s.class(pca$li,as.factor(resultskm[[2]]$consensusClass),col=prettycolors[as.numeric(resultskm[[2]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultskm[[3]]$consensusClass),col=prettycolors[as.numeric(resultskm[[3]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultskm[[4]]$consensusClass),col=prettycolors[as.numeric(resultskm[[4]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultskm[[5]]$consensusClass),col=prettycolors[as.numeric(resultskm[[5]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
s.class(pca$li,as.factor(resultskm[[6]]$consensusClass),col=prettycolors[as.numeric(resultskm[[6]]$consensusClass)],xax = 1,yax=2,addaxes = T,sub= paste("ConsensusClusterPlus", paste(paste("PC",1:2,": ",format(pca$eig[1:2]/sum(pca$eig)*100,digits=2),"%",sep=""),collapse = ", ")) )
dev.off()
