#!/usr/bin/env Rscript

####################################################
### Script to perform Clustering of RNA seq data ###
####################################################

# get options
library("optparse")

option_list = list(
  make_option(c("-r", "--rdata"), type="character", default=".", help=".RData file with clustering and PCA results [default= %default]", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=".", help="input file with groups [default= %default]", metavar="character"),
  make_option(c("-m", "--Kmin"), type="numeric", default=2, help="minimum number of clusters [default= %default]", metavar="numeric"),
  make_option(c("-M", "--Kmax"), type="numeric", default=5, help="maximum number of clusters [default= %default]", metavar="numeric"),
  make_option(c("-n", "--n"), type="numeric", default=5, help="number of permutations [default= %default]", metavar="numeric"),
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# load libraries
library(ade4)
#require(permute)

# define some nice colors
prettycolors = c(1,2,rgb(0,152/255,219/255),rgb(233/255,131/255,0),rgb(0,155/255,118/255),rgb(0.5,0.5,0.5),rgb(0,124/255,146/255),rgb(178/255,111/255,22/255),rgb(234/255,171/255,0),rgb(83/255,40/255,79/255))
contcol      = function(x){
  xx = x
  xx[is.na(x)] = 0
  xnorm = (xx-min(xx,na.rm=T))/(max(xx,na.rm=T)-min(xx,na.rm=T))
  return( rgb(1,0.9-xnorm/1.2,0.9-xnorm/1.2) )
}

# define useful functions
match2 = function(v1i,v2i){
  v1 = as.numeric(v1i)
  v2 = as.numeric(v2i)
  if(max(v1,na.rm = T)>=max(v2,na.rm = T)){ #smallest must be v2
    vtmp=v2
    v2=v1
    v1=vtmp
  }
  ncl = max(v2,na.rm=T)
  
  ltmp = lapply(1:ncl, function(i) 1:max(v1,na.rm=T))
  test = expand.grid(ltmp)
  perms4 = test[apply(test,1,function(x) length(unique(x))==max(v1,na.rm=T) ),]
  match = sum( v1==v2,na.rm=T)
  curr = v2
  for(i in 1:nrow(perms4)){
    tmpclass = v2
    for(j in 1:ncl) tmpclass[v2==j] = perms4[i,j]
    matchtmp = sum( v1==tmpclass,na.rm=T) #find best cluster correspondance
    #print(matchtmp)
    if(matchtmp>match){
      match = matchtmp
      curr = tmpclass
    }
  }
  if(max(as.numeric(v1i),na.rm=T)>=max(as.numeric(v2i),na.rm=T)){
    #print("rev")
    return(list(match,curr,v1))
  }else{
    #print("norev")
    return(list(match,v1,curr))
  }
}

plotmatch <-function(match,rmatch,main,pm,off=0,n){
  nn = length(rmatch)
  h = hist(rmatch/n*100,xlim=c(0,100),ylim=c(0,0.25),probability = T,breaks = seq(0,100,1),main=main,border = NA,col="grey",xlab="Match (%)")
  h$counts[h$mids< (match[[1]]/n*100) ] = 0 
  h$counts = h$counts/nn
  plot(h,probability=T,add=T,col=2,border=NA)
  abline(v=match[[1]]/n*100,col=2,lty=2,lwd=2)
  mtext(format(match[[1]]/n*100,digits = 3), 1,1.9,at = match[[1]]/n*100,col=2)
  text(match[[1]]/n*100+off,0.10,labels = paste("p =",pm),col=2)
}

#ordering
nicesort <- function(v1,v2){
  #1: by Histo
  order1 = sort.int(v1,index.return = T)$ix
  #2: within H, by bueno
  curr = v2
  order2 = order1
  curr1 = v1[order1]
  curr2 = v2[order1]
  for(i in 1:4){
    order2[curr1==i] = order1[curr1==i][sort.int(curr2[curr1==i],index.return = T)$ix]
  }
  return(order2)
}

# load data
#load("merged_TCGA_EGA/RNAseq_unsupervised_km_complete.RData")
load(opt$rdata)
minK = opt$Kmin
maxK = opt$Kmax

#variables = read.table("metadata_merged.txt",h=T)
variables = read.table(opt$input,h=T)
varnames  = colnames(variables)

lims = c(min(pca$li[,1:2]),max(pca$li[,1:2],na.rm=T))

for(k in 1:ncol(variables)){
pdf(paste(opt$out,"Compare_clust_",varnames[k],".pdf",sep=""),w=4*3,h=4*(maxK-1))
par(mfrow=c((maxK-1),3),family="Times")
for(i in minK:maxK){
  plot(-1000,-1000,xlim=lims,ylim=lims,xlab=paste("PC",1,": ",format(pca$eig[1]/sum(pca$eig)*100,digits=2),"%",sep=""), ylab=paste("PC",2,": ",format(pca$eig[2]/sum(pca$eig)*100,digits=2),"%",sep="") ,main="PCA")
  if(!(class(variables[,k]) %in% c("numeric","integer") ) ){
    grsvals = variables[,k]
    grs     = grsvals[!is.na(grsvals)] 
    cctmp   = clusters[[i]]$consensusClass[!is.na(grsvals)]
    mtmp    = match2(cctmp,grs)
    newclass1 = sapply(1:max(clusters[[i]]$consensusClass,na.rm=T), function(j) mtmp[[2]][cctmp==j ][1] )
    newclass2 = sapply(1:max(as.numeric(grs),na.rm=T), function(j) mtmp[[3]][as.numeric(grs)==j ][1] )
    
    #PCA
    s.class(pca$li,as.factor(clusters[[i]]$consensusClass),col=clusters[[i]]$clrs[[3]][newclass1],xax = 1,yax=2,addaxes = T,sub= "" ,add= T,cpoint = 0)
    points(pca$li[,1],pca$li[,2],col=prettycolors[newclass2][as.numeric(variables[,k])],pch=16)
    legend("topright",legend = levels(variables[,k]),col=prettycolors[newclass2][1:max(as.numeric(variables[,k]),na.rm=T)],pch=16)
    
    #clustering
    ordtmp = nicesort(cctmp,grs)
    
    #layout(m)
    plot(-10,-10,xlim=c(-6,16),ylim=c(-length(mtmp[[2]])/10,length(mtmp[[2]])),axes=F,xlab="",ylab="",main=paste("Matching Clusters/",varnames[k],sep="") )
    for(ii in 1:length(mtmp[[2]])) polygon(c(0,0,4.5,4.5),c(ii,ii+1,ii+1,ii),border = NA,col=clusters[[i]]$clrs[[3]][newclass1][cctmp][ordtmp][ii])
    for(ii in 1:length(mtmp[[2]])){if((clusters[[i]]$clrs[[3]][mtmp[[2]]][ordtmp][ii])!=(clusters[[i]]$clrs[[3]][mtmp[[3]]][ordtmp][ii]))  segments(4.6,ii+0.5,5.4 ,ii+0.5 ,col=rgb(1,0,0,0.5)) }
    for(ii in 1:length(mtmp[[2]])) polygon(c(5.5,5.5,10,10),c(ii,ii+1,ii+1,ii),border = NA,col=prettycolors[newclass2][grs][ordtmp][ii])
    legend("left", legend=paste("Cluster ",unique(cctmp),sep="") , fill =clusters[[i]]$clrs[[3]][newclass1],bty = 'n',border=NA)
    legend("right", legend=levels(grs) , fill =prettycolors[newclass2],bty = 'n',border=NA)
    legend("bottom", legend="Mismatch" , lty=1,lwd=2,col=rgb(1,0,0,0.5),bty = 'n',border=NA)
    
    #test
    rmatch = sapply(1:opt$n, function(j){match2(sample(cctmp),grs )[[1]] }  )
  }else{
    #clustering
    grsvals = variables[,k]
    grsvals2= grsvals[!is.na(grsvals)] 
    cuts = seq( min(grsvals2),max(grsvals2,na.rm=T),length.out = i+1)
    grs = cut( grsvals2 ,cuts ,include.lowest = T)
    meds = cuts[1:i] + (cuts[2]-cuts[1])/2
    idmeds = sapply(meds, function(x) which.min(abs(x-grsvals2) )) 
    cctmp   = clusters[[i]]$consensusClass[!is.na(grsvals)]
    
    #PCA
    s.class(pca$li,as.factor(clusters[[i]]$consensusClass),col=clusters[[i]]$clrs[[3]],xax = 1,yax=2,addaxes = T,sub= "" ,add= T,cpoint = 0)
    points(pca$li[,1],pca$li[,2],col=contcol(variables[,k]),cex=(variables[,k]/max(variables[,k],na.rm = T)*2) ,pch=16)
    legend("topright",legend = format(grsvals2[idmeds],digits=2) ,col=contcol(grsvals2)[floor(idmeds)],pch=16,pt.cex = (grsvals2/max(grsvals2)*2)[floor(idmeds)]  )
    
    #clustering
    mtmp   = match2(cctmp,grs)
    ordtmp = nicesort(mtmp[[2]],mtmp[[3]])
    
    #layout(m)
    plot(-10,-10,xlim=c(-6,16),ylim=c(-20,length(mtmp[[2]])),axes=F,xlab="",ylab="",main=paste("Matching Clusters/",varnames[k],sep="") )
    for(ii in 1:length(mtmp[[2]])) polygon(c(0,0,4.5,4.5),c(ii,ii+1,ii+1,ii),border = NA,col=clusters[[i]]$clrs[[3]][cctmp][ordtmp][ii])
    for(ii in 1:length(mtmp[[2]])){if((clusters[[i]]$clrs[[3]][mtmp[[2]]][ordtmp][ii])!=(clusters[[i]]$clrs[[3]][mtmp[[3]]][ordtmp][ii]))  segments(4.6,ii+0.5,5.4 ,ii+0.5 ,col=rgb(1,0,0,0.5)) }
    for(ii in 1:length(mtmp[[2]])) polygon(c(5.5,5.5,10,10),c(ii,ii+1,ii+1,ii),border = NA,col=contcol(meds)[grs][ordtmp][ii])
    legend("bottomleft", legend=paste("Cluster ",unique(cctmp),sep="") , fill =clusters[[i]]$clrs[[3]])
    legend("bottomright", legend=levels(grs) , fill =contcol(meds))
    legend("bottom", legend="Mismatch" , lty=1,lwd=2,col=rgb(1,0,0,0.5))
    
    rmatch = sapply(1:opt$n, function(j){match2(sample(cctmp),grs )[[1]] }  )
  }
  pm     = mean( rmatch>= mtmp[[1]] )
  plotmatch(mtmp,rmatch,paste("Pr of best matching Cluster/",varnames[k],sep=""),pm,n=length(cctmp),off=15)
}

dev.off()
}


