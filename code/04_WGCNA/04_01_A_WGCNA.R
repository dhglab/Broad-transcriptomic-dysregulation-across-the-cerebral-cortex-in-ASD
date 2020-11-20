##### 04_01_A_WGCNA.R
##### Perform WGCNA (Gene-Level)
##### November 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(WGCNA); library(gridExtra); library(ggplot2); library(lmtest)
library(gplots); library(cluster); library(flashClust); library(base)
library(Hmisc); library(limma); library(reshape2); library(grid)
library(statmod)

load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
rm(datExpr)
tdatExpr=t(datExpr.reg)

##### Pipeline #####

###(1) Choose parameters for rWGCNA
###(2) Perform rWGCNA
###(3) Identify module associations with covariates
###(4) Permute ME region-specific ASD effects
###(5) Gene biotype permutation

##### (1) Choose parameters for rWGCNA #####

### Calculate soft power threshold

powers = c(seq(1,9,by=1),seq(10,30,by=2))
powerTable = pickSoftThreshold(data= tdatExpr, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = 5000)

pdf("plots/04_WGCNA/04_01_A_01_CalculateSftPowerThresh.pdf",width=12)

par(mfrow=c(1,2))
sft = powerTable
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], main=names(powerTable), xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
plot.new()
par(mfrow=c(1,1))
grid.table(signif(sft$fitIndices,2),rows=NULL)

dev.off()

### choose soft threshold of 6

### Get TOM matrix 

  pow=6
  
  enableWGCNAThreads()
  net = blockwiseModules(datExpr=tdatExpr, maxBlockSize=25000,networkType="signed",corType="bicor",  
                         power = pow, mergeCutHeight= 0.1, minModuleSize= 50, pamStage=FALSE, reassignThreshold=1e-6, 
                         saveTOMFileBase="data_user/04_WGCNA/04_01_A_01_WGCNA_prelim", saveTOMs=TRUE, 
                         verbose = Inf, deepSplit=4)
  save(net,file="data_user/04_WGCNA/04_01_A_01_WGCNA_prelim_network.RData")
  
### Cut Tree and Plot

load("data_user/04_WGCNA/04_01_A_01_WGCNA_prelim_network.RData")

MEs = net$MEs;
moduleLabels = net$colors;
moduleColors = (moduleLabels)
Tree = net$dendrograms[[1]];

pdf(file="plots/04_WGCNA/04_01_A_01_PreliminaryTreeCut.pdf")

plotDendroAndColors(Tree, moduleColors,
                    "Prelim Modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Prelim ASDPan Dendrogram; CPM w/ cqn offset tech regressed")
dev.off()


### Recut Tree to optimize tree cutting parameters

load("data_user/04_WGCNA/04_01_A_01_WGCNA_prelim-block.1.RData")

geneTree = hclust(1-TOM, method = "average")
colors = vector(mode="list")
labels = vector(mode="list")

for (pam in c(FALSE)) {
  for (minModSize in c(50,100, 200)) {
    for (dthresh in c(0.1, 0.2)) {
      for(ds in c(1,2,3,4)) {
        print(paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,",PAM=",pam,sep=""))
        tree = cutreeHybrid(dendro = geneTree, minClusterSize= minModSize, pamStage=pam, cutHeight = 0.999, deepSplit=ds, distM=as.matrix(1-TOM))
        merged = mergeCloseModules(exprData= tdatExpr, colors = tree$labels, cutHeight=dthresh)
        colors = cbind(colors, labels2colors(merged$colors))
        
        labels = c(labels, paste("DS=", ds, ",MMS=", minModSize, ",DCOR=",dthresh,sep=""))
      }
    }
  }
}

pdf("plots/04_WGCNA/04_01_A_01_ASDPan_DiffParams_PrelimTreeCut.pdf",width=12)
plotDendroAndColors(geneTree, colors, groupLabels=labels, addGuide= TRUE,
                    dendroLabels=FALSE, main="Diff Params ASDPan", cex.colorLabels=0.5)
dev.off()
  
### Choose final tree cutting parameters

ds=4; minModSize=50; dthresh=0.1; pam= FALSE

### Load TOMs, make tree, cut tree, calculate module eigengenes and kME table
networks=vector(mode="list", length=1)
names(networks) = c("datExpr")
networks[[1]]$datExpr = tdatExpr; 
networks[[1]]$TOMfile = "data_user/04_WGCNA/04_01_A_01_WGCNA_prelim-block.1.RData"  

### Plot Tree and Visualize Chosen Tree Cutting Parameters

load("data_user/04_WGCNA/04_01_A_01_WGCNA_prelim_network.RData")

MEs = net$MEs;
moduleLabels = net$colors;
moduleColors = (moduleLabels)
Tree = net$dendrograms[[1]];

load("data_user/04_WGCNA/04_01_A_01_WGCNA_prelim-block.1.RData")

pdf("plots/04_WGCNA/04_01_A_01_ASDPan_TreeCut_ds4-mm50-dthresh01-pamF.pdf",width=10)

for(i in c(1)) {
  networks[[i]]$tree = hclust(1-TOM, method="average")
  networks[[i]]$cut = cutreeHybrid(dendro = networks[[i]]$tree, pamStage=pam, minClusterSize= minModSize, cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-TOM))
  networks[[i]]$merged = mergeCloseModules(exprData= networks[[i]]$datExpr, colors = networks[[i]]$cut$labels, cutHeight=dthresh)
  networks[[i]]$MEs = moduleEigengenes(networks[[i]]$datExpr, colors=networks[[i]]$merged$colors, softPower=6)
  networks[[i]]$kMEtable = signedKME(networks[[i]]$datExpr, datME = networks[[i]]$MEs$eigengenes,corFnc = "bicor")
  plotDendroAndColors(networks[[i]]$tree, colors=labels2colors(networks[[i]]$merged$colors), dendroLabels = F)
}

dev.off()

##### (2) Perform rWGCNA #####

### It is recommended to run this next section on a high-performance computing cluster.

for(ind in c(1:100)){
  
  load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
  rm(datExpr)
  load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
  
  set.seed(ind)
  
  indA=which(datMeta$Diagnosis=="ASD")
  indC=which(datMeta$Diagnosis=="CTL")
  indD=which(datMeta$Diagnosis=="Dup15q")
  
  autInd=sample(indA,length(indA),replace=T)
  conInd=sample(indC,length(indC),replace=T)
  dupInd=sample(indD,length(indD),replace=T)
  
  datMeta_model$sample_id = datMeta$sample_id
  
  grp1=datMeta_model[c(autInd,conInd,dupInd),]
  data1=datExpr.reg[,c(autInd,conInd,dupInd)]
  
  datExpr=as.data.frame(t(data1))
  
  ### final parameters for module discovery
  
  softPower=6
  ds=4
  mms=50
  dcor=0.1
  
  ### Network Analysis Starting
  
  adjacency = adjacency(datExpr, power = softPower, type = "signed",corFnc="bicor");
  
  ### Adjacency Calculation Done
  
  TOM = TOMsimilarity(adjacency);
  
  ### TOM Calculation Done
  
  dissTOM = 1-TOM
  geneTree = flashClust(as.dist(dissTOM), method = "average");
  
  filename1=paste("data_user/04_WGCNA/04_01_A_02_resampleTOMs/resampling_TOMs_",ind,".RData",sep="")
  save(TOM,adjacency,file=filename1)
  
  ### TOM calculation complete
  
  tree = cutreeHybrid(dendro = geneTree,
                      pamStage=F,
                      minClusterSize =mms,
                      cutHeight = 0.9999,
                      deepSplit = ds,
                      distM = as.matrix(dissTOM))
  
  merge <- mergeCloseModules(exprData = datExpr,colors = tree$labels, cutHeight = dcor)
  
  filename=paste("data_user/04_WGCNA/04_01_A_02_data_resampled/resampling_TOMs",ind,".RData",sep="")
  save(geneTree,merge,datExpr,grp1,ind,file=filename)				
  
  ### Run complete	
  
}

set.seed(8675309)
topDir=paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD",sep="")

for(ind in c(1:100)){
  
  ## Get the Reference TOM
  load("data_user/04_WGCNA/04_01_A_01_WGCNA_prelim-block.1.RData")

  TOM.mat <- as.matrix(TOM)
  nGenes <- nrow(TOM.mat)
  
  ## This is a tricky part - we can't load all TOMs simultaneously, so we need to load them piecewise and take the median
  bsize <- 5000
  
  ## Set the general parameters
  prob <- 0.95 ## Define the percentile (called prob) for defining the quantile
  nSubset <- 1000000 ## number of randomly sampled entries of the TOM matrix
  
  ## Choose the sampled TOM entries
  subsetEntries = sample(nGenes*(nGenes-1)/2, size = nSubset)
  
  ## Load the original TOM
  TOMsubset.main <- vectorizeMatrix(TOM.mat)[subsetEntries] ## Select the sampled TOM entries 
  quantile.TOM.main <- quantile(TOMsubset.main,probs=prob,type = 8) ## Calculate the quantile corresponding to prob  
  
  print("Scaling for Reference TOM Complete !")
  
  set=ind
  
  if (!file.exists(paste(topDir,"/data_user/04_WGCNA/04_01_A_02_piecewiseTOMs/TOM",set,"piece",ceil(nGenes/bsize),".RData",sep=""))) {
    setName=dir(paste(topDir,"/data_user/04_WGCNA/04_01_A_02_resampleTOMs/",sep=""))[set]
    load(paste(topDir,"/data_user/04_WGCNA/04_01_A_02_resampleTOMs/",setName,sep=""))
    tmpTOM <- as.matrix(TOM)
    TOMsubset <- vectorizeMatrix(tmpTOM)
    quantile.TOM <- quantile(TOMsubset,probs=prob,type = 8) ## Calculate the quantile corresponding to prob
    ## vector of quantiles of the individual TOM matrices
    beta.prob = log(quantile.TOM.main)/log(quantile.TOM) ## calculate the power of the adjacency function
    save(beta.prob,file=paste(topDir,"/data_user/04_WGCNA/04_01_A_02_ConsensusTOMscalingQuantiles/iteration",ind,".RData",sep=""))
    ## Scaling powers to equalize reference TOM values
    tmpTOM <- tmpTOM^(beta.prob) ## use the power adjacency function for calibrating the TOM
    
    for (i in 1:ceil(nGenes/bsize)) {
      if (i < ceil(nGenes/bsize)) {
        subTOM <- tmpTOM[,c((i-1)*bsize+1):(i*bsize)]
        save(subTOM,file=paste(topDir,"/data_user/04_WGCNA/04_01_A_02_piecewiseTOMs/TOM",set,"piece",i,".RData",sep=""))
      } else {
        subTOM <- tmpTOM[,c((i-1)*bsize+1):nGenes]
        save(subTOM,file=paste(topDir,"/data_user/04_WGCNA/04_01_A_02_piecewiseTOMs/TOM",set,"piece",i,".RData",sep=""))
      }
    }
  } else {
    print(dir(paste(topDir,"/data_user/04_WGCNA/04_01_A_02_resampleTOMs/",sep=""))[set])
    print("Already processed!")
  }
}

for(p in c(1:ceil(nGenes/bsize))){
  
  consSubTOM <- vector("list",100)
  
  for (set in 1:100) {
    gc()
    load(paste(topDir,"/data_user/04_WGCNA/04_01_A_02_piecewiseTOMs/TOM",set,"piece",p,".RData",sep=""))
    consSubTOM[[set]]$subTOM <- subTOM
  }
  consensusTOM <- matrix(NA,nrow=nrow(subTOM),ncol=ncol(subTOM))
  for (i in 1:nrow(subTOM) ) {
    gc()
    if (i%%100 == 0) { print(paste("On gene",i,sep=" ")) }          
    for (j in 1:ncol(subTOM)) {
      thisvec <- c(consSubTOM[[1]]$subTOM[i,j], consSubTOM[[2]]$subTOM[i,j], consSubTOM[[3]]$subTOM[i,j], consSubTOM[[4]]$subTOM[i,j], consSubTOM[[5]]$subTOM[i,j], consSubTOM[[6]]$subTOM[i,j], consSubTOM[[7]]$subTOM[i,j], consSubTOM[[8]]$subTOM[i,j], consSubTOM[[9]]$subTOM[i,j], consSubTOM[[10]]$subTOM[i,j], consSubTOM[[11]]$subTOM[i,j], consSubTOM[[12]]$subTOM[i,j], consSubTOM[[13]]$subTOM[i,j], consSubTOM[[14]]$subTOM[i,j], consSubTOM[[15]]$subTOM[i,j], consSubTOM[[16]]$subTOM[i,j], consSubTOM[[17]]$subTOM[i,j], consSubTOM[[18]]$subTOM[i,j], consSubTOM[[19]]$subTOM[i,j], consSubTOM[[20]]$subTOM[i,j], consSubTOM[[21]]$subTOM[i,j], consSubTOM[[22]]$subTOM[i,j], consSubTOM[[23]]$subTOM[i,j], consSubTOM[[24]]$subTOM[i,j], consSubTOM[[25]]$subTOM[i,j], consSubTOM[[26]]$subTOM[i,j], consSubTOM[[27]]$subTOM[i,j], consSubTOM[[28]]$subTOM[i,j], consSubTOM[[29]]$subTOM[i,j], consSubTOM[[30]]$subTOM[i,j], consSubTOM[[31]]$subTOM[i,j], consSubTOM[[32]]$subTOM[i,j], consSubTOM[[33]]$subTOM[i,j], consSubTOM[[34]]$subTOM[i,j], consSubTOM[[35]]$subTOM[i,j], consSubTOM[[36]]$subTOM[i,j], consSubTOM[[37]]$subTOM[i,j], consSubTOM[[38]]$subTOM[i,j], consSubTOM[[39]]$subTOM[i,j], consSubTOM[[40]]$subTOM[i,j], consSubTOM[[41]]$subTOM[i,j], consSubTOM[[42]]$subTOM[i,j], consSubTOM[[43]]$subTOM[i,j], consSubTOM[[44]]$subTOM[i,j], consSubTOM[[45]]$subTOM[i,j], consSubTOM[[46]]$subTOM[i,j], consSubTOM[[47]]$subTOM[i,j], consSubTOM[[48]]$subTOM[i,j], consSubTOM[[49]]$subTOM[i,j], consSubTOM[[50]]$subTOM[i,j], consSubTOM[[51]]$subTOM[i,j], consSubTOM[[52]]$subTOM[i,j], consSubTOM[[53]]$subTOM[i,j], consSubTOM[[54]]$subTOM[i,j], consSubTOM[[55]]$subTOM[i,j], consSubTOM[[56]]$subTOM[i,j], consSubTOM[[57]]$subTOM[i,j], consSubTOM[[58]]$subTOM[i,j], consSubTOM[[59]]$subTOM[i,j], consSubTOM[[60]]$subTOM[i,j], consSubTOM[[61]]$subTOM[i,j], consSubTOM[[62]]$subTOM[i,j], consSubTOM[[63]]$subTOM[i,j], consSubTOM[[64]]$subTOM[i,j], consSubTOM[[65]]$subTOM[i,j], consSubTOM[[66]]$subTOM[i,j], consSubTOM[[67]]$subTOM[i,j], consSubTOM[[68]]$subTOM[i,j], consSubTOM[[69]]$subTOM[i,j], consSubTOM[[70]]$subTOM[i,j], consSubTOM[[71]]$subTOM[i,j], consSubTOM[[72]]$subTOM[i,j], consSubTOM[[73]]$subTOM[i,j], consSubTOM[[74]]$subTOM[i,j], consSubTOM[[75]]$subTOM[i,j], consSubTOM[[76]]$subTOM[i,j], consSubTOM[[77]]$subTOM[i,j], consSubTOM[[78]]$subTOM[i,j], consSubTOM[[79]]$subTOM[i,j], consSubTOM[[80]]$subTOM[i,j], consSubTOM[[81]]$subTOM[i,j], consSubTOM[[82]]$subTOM[i,j], consSubTOM[[83]]$subTOM[i,j], consSubTOM[[84]]$subTOM[i,j], consSubTOM[[85]]$subTOM[i,j], consSubTOM[[86]]$subTOM[i,j], consSubTOM[[87]]$subTOM[i,j], consSubTOM[[88]]$subTOM[i,j], consSubTOM[[89]]$subTOM[i,j], consSubTOM[[90]]$subTOM[i,j], consSubTOM[[91]]$subTOM[i,j], consSubTOM[[92]]$subTOM[i,j], consSubTOM[[93]]$subTOM[i,j], consSubTOM[[94]]$subTOM[i,j], consSubTOM[[95]]$subTOM[i,j], consSubTOM[[96]]$subTOM[i,j], consSubTOM[[97]]$subTOM[i,j], consSubTOM[[98]]$subTOM[i,j], consSubTOM[[99]]$subTOM[i,j], consSubTOM[[100]]$subTOM[i,j])
      consensusTOM[i,j] <- median(thisvec)
    }
  }
  save(consensusTOM,file=paste(topDir,"/data_user/04_WGCNA/04_01_A_02_piecewiseTOMs/consTOMs/ConsensusTOM_nSets100_piece",p,".RData",sep=""))
  
  cat('Done .... \n')
  
}

topDir=paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/data_user/04_01_A_02_piecewiseTOMs/consTOMs",sep="")
setwd(topDir)

consensusTOM_final=matrix(NA,nrow=24836,ncol=1)

for (i in 1:5) ##5 comes from ceil(nGenes/bsize)
{
  filename=paste("ConsensusTOM_nSets100_piece",i,".RData",sep="")
  load(filename);
  cat(dim(consensusTOM),'\n')
  consensusTOM_final=cbind(consensusTOM_final,consensusTOM)
  cat(dim(consensusTOM_final),'\n')
  
}

consensusTOM_final= consensusTOM_final[,-1]

dim(consensusTOM_final)  ###should be a square matrix

save(consensusTOM_final,file="04_01_A_02_ConsensusTOM_final.RData")

### Analyze consensus TOM

rm(list=ls())

### First, get resampled expression and metadata

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/data_user/04_01_A_02_data_resampled",sep=""))

nSets=100 ### from resampling

### Loading Expression Data

multiExpr=vector(mode="list",length=nSets)
for (i in 1:nSets){
  load(dir()[i])
  multiExpr[[i]] = list(data=as.data.frame(datExpr)) 
  names(multiExpr[[i]]$data)=colnames(datExpr)
  rownames(datExpr) = make.names(rownames(datExpr), unique=TRUE)
  rownames(multiExpr[[i]]$data)=rownames(datExpr)
  cat('Done for data ....',i,'\n')	
}

### Loading Meta Data

multiMeta=vector(mode="list",length=nSets)
for (i in 1: nSets){
  load(dir()[i])
  multiMeta[[i]] = list(data=as.data.frame(grp1)) 
  names(multiMeta[[i]]$data)=colnames(grp1)
  rownames(grp1) = make.names(rownames(grp1), unique=TRUE)
  rownames(multiMeta[[i]]$data)=rownames(grp1)
  cat('Done for data ....',i,'\n')	
}

checkSets(multiExpr) # check data size
checkSets(multiMeta) # check data size

setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD",sep=""))

setwd("../")
rm(datExpr,geneTree,i,ind,merge,grp1)

save(list=ls(),file="data_user/04_WGCNA/01_01_A_02_consTOM_results/MultiData_Resample.RData")

### Make dendrogram with all 100 different resampled TOMs as rows underneath

rm(list=ls())

wkdir="C:/Users/jillh/Dropbox/GitHub/"
load("data_user/04_WGCNA/04_01_A_02_piecewiseTOMs/consTOMs/04_01_A_02_ConsensusTOM_final.RData")
load("data_user/04_WGCNA/04_01_A_02_consTOM_results/MultiData_Resample.RData")

consTree = hclust(as.dist(1-consensusTOM_final), method = "average");

nsets=100
minModSize = 50
ds = 4
dthresh = 0.1
mColorh <- mLabelh <- colorLabels <- NULL

### This code gives us a reference module to use when matching modules later

print("getting tree")
treeCons = cutreeHybrid(dendro = consTree, pamStage=FALSE,
                        minClusterSize = minModSize, cutHeight = 0.9999,
                        deepSplit = ds, distM = as.matrix(1-consensusTOM_final))

## consMerged takes awhile
print("merging close modules")
consMerged <- mergeCloseModules(exprData = multiExpr,colors = treeCons$labels,
                                cutHeight = dthresh)

consMergedCol = labels2colors(consMerged$colors)

mColorh <- cbind(mColorh,consMergedCol)
mLabelh <- c(mLabelh,"ConsensusMods")

### Find all 100 different module distributions from 100 resamples TOMs

print("getting modules from 100 resampled TOMs")

for (i in 1:100){
  print(paste("Loading Data for Set ",i,sep=""))
  load(paste("data_user/04_WGCNA/04_01_A_02_data_resampled/resampling_TOMs",i,".RData",sep=""))
  print(paste("Starting TOM ",i,sep=""))
  consColors <- matchLabels(labels2colors(merge$colors),consMergedCol)
  mColorh <- cbind(mColorh,consColors)
  mLabelh <- c(mLabelh,paste("TOM",i,sep=" "))
  print("Done")
}

save(consTree,treeCons,consMerged,mColorh,mLabelh,file="data_user/04_WGCNA/04_01_A_02_consTOM_results/ConsensusWGCNAData.RData")

### Plot results

pdf("plots/04_WGCNA/04_01_A_02_TOMs_Dendrogram_MergeCols_wConsMods.pdf",height=50,width=20)

plotDendroAndColors(consTree,mColorh,groupLabels=mLabelh,addGuide=TRUE,
                    dendroLabels=FALSE,main="Dendrogram with 100 Resampled TOMs",
                    cex.colorLabels=0.5)

dev.off()

png("plots/04_WGCNA/04_01_A_02_TOMs_Dendrogram_MergeCols_wConsMods.png",width = 800, height = 1250, units = "px", pointsize = 12)

plotDendroAndColors(consTree,mColorh,groupLabels=mLabelh,addGuide=TRUE,
                    dendroLabels=FALSE,main="Dendrogram with 100 Resampled TOMs",
                    cex.colorLabels=0.5)

dev.off()

## also make dendro with just the cons TOM

pdf("plots/04_WGCNA/04_01_A_02_ConsTOM_Dendrogram_MergeCols.pdf",width=20,height=15)

plotDendroAndColors(consTree,mColorh[,1],groupLabels=mLabelh[1],addGuide=TRUE,
                    dendroLabels=FALSE,main="Dendrogram with Consensus TOM",
                    cex.colorLabels=1)

dev.off()

##### (3) Identify module associations with covariates #####

rm(list=ls())
wkdir="C:/Users/jillh/Dropbox/GitHub/"
load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

load("data_user/04_WGCNA/04_01_A_02_consTOM_results/ConsensusWGCNAData.RData")

colors = consMerged$colors
color_names = labels2colors(colors)

MEs = moduleEigengenes(t(datExpr.reg),colors=colors, softPower=6)
MEs = MEs$eigengenes

### make a key - will name modules accordingly. This is important since some modules were merged.

module_key = data.frame("Color"=color_names,
                        "Number_Original"=colors)

module_key = module_key[-which(duplicated(module_key$Color)==TRUE),]
module_key = module_key[order(module_key$Number_Original),]
module_key$Number_Finalized = c(0:35)
module_key$Module_Name = paste("M",module_key$Number_Finalized,"_",module_key$Color,sep="")
module_key$ME_Name = paste("ME",module_key$Number_Finalized,"_",module_key$Color,sep="")

labs <- as.numeric(gsub("ME","",colnames(MEs)))
cols <- labels2colors(labs)
colnames(MEs) <- paste("ME",module_key$Number_Finalized[match(labs,module_key$Number_Original)],"_",cols,sep = "")
rownames(MEs) <- colnames(datExpr.reg)

kME = signedKME(t(datExpr.reg), datME = MEs,corFnc = "bicor")

modules = list()

for(col in unique(color_names)){
  modules[[col]]=rownames(datExpr.reg)[which(color_names==col)]
}

save(colors,color_names,MEs,kME,modules,module_key,file="data_user/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")


### Cluster MEs

# Calculate dissimilarity of module eigengenes
MEs = MEs[,-1]
MEDiss = 1-bicor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

pdf(file="plots/04_WGCNA/04_01_A_03_ME-Tree.pdf",width=12,height=10)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

ME_order = METree$order

### get whole cortex and region specific asd effects, as well as effects of age and sex 

library(limma)
design = model.matrix(~ 0 + DxReg + Sex + Ancestry + Age + Age_sqd, data = datMeta_model)

corfit= duplicateCorrelation(t(MEs),design,block=datMeta_model$Subject)
lm = lmFit(t(MEs), design,block=datMeta_model$Subject, correlation = corfit$consensus)

this_contrast_asd = makeContrasts(contrast=(DxRegASD_BA17 + DxRegASD_BA20_37 + DxRegASD_BA24 + DxRegASD_BA3_1_2_5 +
                                              DxRegASD_BA38 + DxRegASD_BA39_40 + DxRegASD_BA4_6 + DxRegASD_BA41_42_22 + 
                                              DxRegASD_BA44_45 + DxRegASD_BA7 + DxRegASD_BA9 - DxRegCTL_BA17 - 
                                              DxRegCTL_BA20_37 - DxRegCTL_BA24 - DxRegCTL_BA3_1_2_5 - DxRegCTL_BA38 - 
                                              DxRegCTL_BA39_40 - DxRegCTL_BA4_6 - DxRegCTL_BA41_42_22 - DxRegCTL_BA44_45 - 
                                              DxRegCTL_BA7 - DxRegCTL_BA9)/11, levels=design)

this_contrast_dup15q = makeContrasts(contrast=(DxRegDup15q_BA17 + DxRegDup15q_BA20_37 + DxRegDup15q_BA24 + DxRegDup15q_BA3_1_2_5 + 
                                                 DxRegDup15q_BA38 + DxRegDup15q_BA39_40 + DxRegDup15q_BA4_6 + DxRegDup15q_BA41_42_22 + 
                                                 DxRegDup15q_BA44_45 + DxRegDup15q_BA7 + DxRegDup15q_BA9 - DxRegCTL_BA17 - 
                                                 DxRegCTL_BA20_37 - DxRegCTL_BA24 - DxRegCTL_BA3_1_2_5 - DxRegCTL_BA38 - 
                                                 DxRegCTL_BA39_40 - DxRegCTL_BA4_6 - DxRegCTL_BA41_42_22 - DxRegCTL_BA44_45 - 
                                                 DxRegCTL_BA7 - DxRegCTL_BA9)/11, levels=design)

### Make tables for association of MEs with all traits

### first: whole cortex

assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,ncol(MEs),5))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-  c("SexF","Age","Age_sqd","Diagnosis_ASD","Diagnosis_Dup15q")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  colnames(MEs)

fit = eBayes(lm,trend=T,robust=T)
coefs = colnames(fit$coefficients)

for(j in c(1:ncol(assoc_table_p))){
  if(length(grep("Diagnosis",colnames(assoc_table_p)[j])) > 0){
    if(j==4){
      fit_asd = contrasts.fit(lm, this_contrast_asd)
      fit2_asd= eBayes(fit_asd,trend = T, robust = T)
      tt_ASD_Region= topTable(fit2_asd, coef=1, number=Inf, sort.by = 'none')
      for(me in colnames(MEs)){
        assoc_table_p[me,j] = signif(tt_ASD_Region$adj.P.Val[which(rownames(tt_ASD_Region)==me)],3)
        assoc_table_beta[me,j] = signif(tt_ASD_Region$logFC[which(rownames(tt_ASD_Region)==me)],3)
      }
    }else if(j==5){
      fit_dup = contrasts.fit(lm, this_contrast_dup15q)
      fit2_dup= eBayes(fit_dup,trend = T, robust = T)
      tt_Dup15q_Region= topTable(fit2_dup, coef=1, number=Inf, sort.by = 'none')
      for(me in colnames(MEs)){
        assoc_table_p[me,j] = signif(tt_Dup15q_Region$adj.P.Val[which(rownames(tt_Dup15q_Region)==me)],3)
        assoc_table_beta[me,j] = signif(tt_Dup15q_Region$logFC[which(rownames(tt_Dup15q_Region)==me)],3)
      }    
    }
  }else{
    idx=grep(paste("\\b",colnames(assoc_table_p)[j],"\\b",sep=""),coefs)
    tt_tmp = topTable(fit,coef=idx,number=Inf, sort.by = 'none')
    for(me in colnames(MEs)){
      assoc_table_p[me,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_tmp)==me)],3)
      assoc_table_beta[me,j] = signif(tt_tmp$logFC[which(rownames(tt_tmp)==me)],3)
    }
  }
}

assoc_table_p_wc = assoc_table_p
assoc_table_beta_wc = assoc_table_beta

### next: just get region specific effects (one matrix for ASD, and one for dup15q)

## ASD

assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,ncol(MEs),11))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-  c("BA9","BA24","BA44_45","BA4_6","BA3_1_2_5","BA7",
                                                            "BA39_40","BA17","BA41_42_22","BA20_37","BA38")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  colnames(MEs)

for(j in c(1:ncol(assoc_table_p))){
  reg = colnames(assoc_table_p)[j]
  form = paste("DxRegASD_",reg," - DxRegCTL_",reg,sep="")
  this_contrast = makeContrasts(contrast=form,levels=design)
  fit = contrasts.fit(lm, this_contrast)
  fit2= eBayes(fit,trend = T, robust = T)
  tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
  for(me in colnames(MEs)){
    assoc_table_p[me,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_ASD_Region)==me)],3)
    assoc_table_beta[me,j] = signif(tt_tmp$logFC[which(rownames(tt_ASD_Region)==me)],3)
  }
}

assoc_table_p_asdReg = assoc_table_p
assoc_table_beta_asdReg = assoc_table_beta

## Dup15q

assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,ncol(MEs),11))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-  c("BA9","BA24","BA44_45","BA4_6","BA3_1_2_5","BA7",
                                                            "BA39_40","BA17","BA41_42_22","BA20_37","BA38")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  colnames(MEs)

for(j in c(1:ncol(assoc_table_p))){
  reg = colnames(assoc_table_p)[j]
  form = paste("DxRegDup15q_",reg," - DxRegCTL_",reg,sep="")
  this_contrast = makeContrasts(contrast=form,levels=design)
  fit = contrasts.fit(lm, this_contrast)
  fit2= eBayes(fit,trend = T, robust = T)
  tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
  for(me in colnames(MEs)){
    assoc_table_p[me,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_ASD_Region)==me)],3)
    assoc_table_beta[me,j] = signif(tt_tmp$logFC[which(rownames(tt_ASD_Region)==me)],3)
  }
}

assoc_table_p_dup15Reg = assoc_table_p
assoc_table_beta_dup15Reg = assoc_table_beta

assoc_table_list = list("WC_P"=assoc_table_p_wc,"WC_B"=assoc_table_beta_wc,
                        "ASDReg_P"=assoc_table_p_asdReg,"ASDReg_B"=assoc_table_beta_asdReg,
                        "Dup15Reg_P"=assoc_table_p_dup15Reg,"Dup15Reg_B"=assoc_table_beta_dup15Reg)

save(assoc_table_list,file="data_user/04_WGCNA/04_01_A_03_ME_Diagnosis_by_Region_Effects.RData")


##### (4) Permute ME region-specific ASD effects #####

rm(list=ls())
wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

### Run Permutation

for(iter in c(1:1000)){
  
  load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
  load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
  load("data_user/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")
  ### or: load("data_provided/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")
  
  MEs = MEs[,-1]
  
  datMeta_model$Diagnosis = factor(datMeta$Diagnosis,levels=c("CTL","ASD","Dup15q"))
  datMeta_model$Region = as.factor(datMeta$region)
  
  ### scramble regional identity
  
  if(iter==1){
    datMeta_model_perm = datMeta_model
  }else{
    set.seed(as.numeric(iter-1))
    datMeta_model_perm = datMeta_model
    for(sub in unique(datMeta_model$Subject)){
      idx = which(datMeta_model$Subject==sub)
      datMeta_model_perm$Region[idx]=sample(datMeta_model_perm$Region[idx])
    }
  }
  
  rm(datMeta_model)
  
  library(limma)
  
  datMeta_model_perm$Region=gsub("-",".",datMeta_model_perm$Region)
  datMeta_model_perm$DxRegion = factor(paste0(datMeta_model_perm$Diagnosis, ".", datMeta_model_perm$Region))
  design = model.matrix(~ 0 + DxRegion + Sex + Ancestry + Age + Age_sqd, data = datMeta_model_perm)
  
  corfit= duplicateCorrelation(t(MEs),design,block=datMeta_model_perm$Subject)
  lm = lmFit(t(MEs), design,block=datMeta_model_perm$Subject, correlation = corfit$consensus)
  
  this_contrast_asd = makeContrasts(contrast=(DxRegionASD.BA17 + DxRegionASD.BA20.37 + DxRegionASD.BA24 + DxRegionASD.BA3.1.2.5 +
                                                DxRegionASD.BA38 + DxRegionASD.BA39.40 + DxRegionASD.BA4.6 + DxRegionASD.BA41.42.22 +
                                                DxRegionASD.BA44.45 + DxRegionASD.BA7 + DxRegionASD.BA9 - DxRegionCTL.BA17 -
                                                DxRegionCTL.BA20.37 - DxRegionCTL.BA24 - DxRegionCTL.BA3.1.2.5 - DxRegionCTL.BA38 -
                                                DxRegionCTL.BA39.40 - DxRegionCTL.BA4.6 - DxRegionCTL.BA41.42.22 - DxRegionCTL.BA44.45 -
                                                DxRegionCTL.BA7 - DxRegionCTL.BA9)/11, levels=design)
  
  ### Make tables for association of MEs with all traits
  
  ### just get region specific effects of ASD
  
  assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,ncol(MEs),11))
  colnames(assoc_table_p) <- colnames(assoc_table_beta) <- c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                                                             "BA39.40","BA17","BA41.42.22","BA20.37","BA38")
  rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  colnames(MEs)
  
  for(j in c(1:ncol(assoc_table_p))){
    reg = colnames(assoc_table_p)[j]
    form = paste("DxRegionASD.",reg," - DxRegionCTL.",reg,sep="")
    this_contrast = makeContrasts(contrast=form,levels=design)
    fit = contrasts.fit(lm, this_contrast)
    fit2= eBayes(fit,trend = T, robust = T)
    tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
    for(me in colnames(MEs)){
      assoc_table_p[me,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_tmp)==me)],3)
      assoc_table_beta[me,j] = signif(tt_tmp$logFC[which(rownames(tt_tmp)==me)],3)
    }
  }
  
  save(assoc_table_p,assoc_table_beta,file=paste("data_user/04_WGCNA/04_01_A_04_permutations/ME_permutation_",iter,".RData",sep=""))
  
}

### (2) Compile Permutation Results

### iter=1 is the real result

rm(list=ls())
wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))
load("data_user/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")

MEs = MEs[,-1]

mes = colnames(MEs)

me_reg_effect_list = list()

for(i in c(1:10001)){
  if (i%%100 == 0) {print(paste(i,"/",length(list.files("permutations/")),sep=""))}
  load(paste("data_user/04_WGCNA/04_01_A_04_permutations/ME_permutation_",i,".RData",sep=""))
  for(me in mes){
    for(reg in c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                 "BA39.40","BA17","BA41.42.22","BA20.37","BA38")){
      if(i==1){
        me_reg_effect_list[[me]][[reg]]$p=assoc_table_p[me,reg]
        me_reg_effect_list[[me]][[reg]]$beta=assoc_table_beta[me,reg]
      }else{
        me_reg_effect_list[[me]][[reg]]$p=c(me_reg_effect_list[[me]][[reg]]$p,assoc_table_p[me,reg])
        me_reg_effect_list[[me]][[reg]]$beta=c(me_reg_effect_list[[me]][[reg]]$beta,assoc_table_beta[me,reg])
      }
    }
  }
}

save(me_reg_effect_list,file="data_user/04_WGCNA/04_01_A_04_ME_Reg_Effect_Permutations.RData")

### (3) Test Severity

rm(list=ls())
wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))
  
load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
load("data_user/04_WGCNA/04_01_A_04_ME_Reg_Effect_Permutations.RData")
load("data_user/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")
  
MEs = MEs[,-1]

mes = colnames(MEs)

# Make tables for association of MEs with all traits

### just get region specific effects of ASD

assoc_table_p <- assoc_table_beta <-  data.frame(matrix(NA,ncol(MEs),11))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-   c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                                                             "BA39.40","BA17","BA41.42.22","BA20.37","BA38")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-   colnames(MEs)

for(j in c(1:ncol(assoc_table_p))){
  reg = colnames(assoc_table_p)[j]
  for(me in colnames(MEs)){
    assoc_table_p[me,reg] = signif(me_reg_effect_list[[me]][[reg]]$p[1],3)
    assoc_table_beta[me,reg] = signif(me_reg_effect_list[[me]][[reg]]$beta[1],3)
  }
}

### now make plots for the individual regions

### ASD

ME_names=factor(rownames(assoc_table_p),levels = mes)

assoc_table_p_tmp = data.frame(MEs=ME_names,assoc_table_p)
melt_assoc_table_p <- melt(assoc_table_p_tmp, id.vars = "MEs")
colnames(melt_assoc_table_p) <- c("ME","Covariate","Value")

ME_assoc_table=melt_assoc_table_p
ME_assoc_table$Value = ME_assoc_table$Value

assoc_table_beta_tmp = data.frame(MEs=ME_names,assoc_table_beta)
melt_assoc_table_beta <- melt(assoc_table_beta_tmp, id.vars = "MEs")
colnames(melt_assoc_table_beta) <- c("ME","Covariate","Value")

ME_assoc_table$Beta <- melt_assoc_table_beta$Value

ME_assoc_table$Signed_log10_PValue <- signif((sign(ME_assoc_table$Beta)*-log10(ME_assoc_table$Value)),3)
colnames(ME_assoc_table)[3] <- "PValue"

idx = which(signif(ME_assoc_table$PValue,2) <= 0.05)
ME_assoc_table$Sig <- rep(NA,dim(ME_assoc_table)[1])
ME_assoc_table$Sig[idx] = "*"

ME_assoc_table_test = ME_assoc_table[idx,]

ME_assoc_table_test$Perm_P = rep(NA,dim(ME_assoc_table_test)[1])

### make plots by region (so 11 PDFs)

for(reg in unique(ME_assoc_table_test$Covariate)){
  
  print(reg)
  pdf(file=paste("plots/04_WGCNA/04_01_A_04_ME_Permutation_",reg,".pdf",sep=""))
  idx = which(ME_assoc_table_test$Covariate==reg)
  
  for(i in idx){
    me=ME_assoc_table_test$ME[i]
    reg=ME_assoc_table_test$Covariate[i]
    tmp=me_reg_effect_list[[as.character(me)]][[as.character(reg)]]
    perm_p=-log10(tmp$p[-1])*sign(tmp$beta[-1])
    real_p=-log10(ME_assoc_table_test$PValue[i])*sign(ME_assoc_table_test$Beta[i])
    
    cent=scale(perm_p,scale=FALSE)
    mn=mean(perm_p)
    diff=real_p-mn
    if(sign(real_p) > 0){
      more=length(which(cent >= diff))
    }else{
      more=length(which(cent <= diff)) 
    }
    
    p = (more)/10001
    ME_assoc_table_test$Perm_P[i] = signif(p,5)
    
    plot.perm=data.frame(perm_p)
    
    print(ggplot(plot.perm,aes(x=perm_p))+
            geom_density()+
            ggtitle(paste(me,", ",reg,"; p=",signif(p,3),sep=""))+
            xlab("Signed log10 Permutation P-Value")+
            theme(plot.title = element_text(hjust = 0.5))+
            geom_vline(xintercept=real_p,col="red"))
  }
  
  dev.off()
  
}

ME_assoc_table_test$Perm_sig = rep(NA,dim(ME_assoc_table_test)[1])
ME_assoc_table_test$Perm_sig[which(round(ME_assoc_table_test$Perm_P,2) <= 0.05)] = "*"

ME_assoc_table_test$FDR_Perm_P = p.adjust(ME_assoc_table_test$Perm_P,method = "fdr")

ME_assoc_table_perm_test = ME_assoc_table_test

save(ME_assoc_table_perm_test,file="data_user/04_WGCNA/04_01_A_04_ME_Permutation_Results.RData")
  
##### (5) Gene biotype permutation #####

rm(list=ls())

allowWGCNAThreads()

wkdir="C:/Users/jillh/Dropbox/GitHub/"
load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

load("data_user/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")
### or: load("data_provided/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")

### Count gene biotypes in modules

RNA_types = list()
all_genes = read.csv("main_datasets/A_GeneLevel/AllGenesByModule_wDGE_wAnnotation.csv")

for(mod in module_key$Module_Name[-1]){
  print(mod)
  tmp = all_genes[which(all_genes$WGCNA_module==mod),]
  RNA_types[[mod]]$length = nrow(tmp)
  RNA_types[[mod]]$types = table(tmp$gene_biotype)
}

uniq_types = NA

for(i in c(1:35)){
  tmp_types = names(RNA_types[[i]]$types)
  uniq_types = c(uniq_types,tmp_types)
  uniq_types = uniq_types[!duplicated(uniq_types)]
}

uniq_types = uniq_types[-1]

RNA_type_table = matrix(NA,nrow=35,ncol=length(uniq_types))
rownames(RNA_type_table) = module_key$Module_Name[-1]
colnames(RNA_type_table) = uniq_types

for(i in c(1:length(uniq_types))){
  for(j in c(1:35)){
    mod = rownames(RNA_type_table)[j]
    type = colnames(RNA_type_table)[i]
    if( type %in% names(RNA_types[[mod]]$types) ){
      RNA_type_table[j,i] = RNA_types[[mod]]$types[which(names(RNA_types[[mod]]$types) == type)]
    }else{
      RNA_type_table[j,i] = 0
    }
  }
}

RNA_type_table_true = RNA_type_table
save(uniq_types,RNA_type_table_true,file="data_user/04_WGCNA/04_01_A_05_RNA_biotypes_trueCount.RData")
  
### Run permutation analysis

### make random module assignments
### It is recommended to run this next section on a high-performance computing cluster.
  
for(iter in c(1:10000)){
  
  set.seed(as.numeric(iter))

  RNA_types = list()
  all_genes = read.csv("main_datasets/A_GeneLevel/AllGenesByModule_wDGE_wAnnotation.csv")
  all_genes$WGCNA_module = sample(all_genes$WGCNA_module,replace=FALSE,size=nrow(all_genes))
  
  for(mod in module_key$Module_Name[-1]){
    print(mod)
    tmp = all_genes[which(all_genes$WGCNA_module==mod),]
    RNA_types[[mod]]$length = nrow(tmp)
    RNA_types[[mod]]$types = table(tmp$gene_biotype)
  }
  
  RNA_type_table = matrix(NA,nrow=35,ncol=length(uniq_types))
  rownames(RNA_type_table) = module_key$Module_Name[-1]
  colnames(RNA_type_table) = uniq_types
  
  for(i in c(1:length(uniq_types))){
    for(j in c(1:35)){
      mod = rownames(RNA_type_table)[j]
      type = colnames(RNA_type_table)[i]
      if( type %in% names(RNA_types[[mod]]$types) ){
        RNA_type_table[j,i] = RNA_types[[mod]]$types[which(names(RNA_types[[mod]]$types) == type)]
      }else{
        RNA_type_table[j,i] = 0
      }
    }
  }
  
  RNA_type_table_perm = RNA_type_table
  save(RNA_type_table_perm,file=paste("data_user/04_WGCNA/04_01_A_05_permutations/RNA_biotypes_perm",iter,".RData",sep=""))
  
}

### Compile output
  
perm_dist = list()

for(i in c(1:10000)){
  load(paste("data_user/04_WGCNA/04_01_A_05_permutations/RNA_biotypes_perm",i,".RData",sep=""))
  for(mod in module_key$Module_Name[-1]){
    for(type in uniq_types){
      if(i==1){
        perm_dist[[mod]][[type]] <- rep(NA,5)
        perm_dist[[mod]][[type]] <- c(perm_dist[[mod]][[type]],RNA_type_table_perm[mod,type])
      }else if(i!=1 & i!=10000){
        perm_dist[[mod]][[type]] <- c(perm_dist[[mod]][[type]],RNA_type_table_perm[mod,type])
      }else if(i==10000){
        perm_dist[[mod]][[type]] <- c(perm_dist[[mod]][[type]],RNA_type_table_perm[mod,type])
        perm_dist[[mod]][[type]] <- perm_dist[[mod]][[type]][-c(1:5)]
      }
    }
  }
}

save(perm_dist,file="data_user/04_WGCNA/04_01_A_05_BiotypePermutationResults.RData")

### Examine results
  
biotype_perm_p = matrix(NA,nrow=35,ncol=length(uniq_types))
rownames(biotype_perm_p) = module_key$Module_Name[-1]
colnames(biotype_perm_p) = uniq_types

for(i in c(1:35)){
  for(j in c(1:length(uniq_types))){
    mod = rownames(biotype_perm_p)[i]
    type = colnames(biotype_perm_p)[j]
    perm_dist_iter = perm_dist[[mod]][[type]]
    true_count = RNA_type_table_true[mod,type]
    more=length(which(perm_dist_iter >= true_count))
    p = (more)/10001
    biotype_perm_p[i,j] = p
  }
}

biotype_perm_fdr <- biotype_perm_p

for(i in c(1:35)){
  idx = which(biotype_perm_p[i,-which(colnames(biotype_perm_p)=="protein_coding")] < 0.05)
  if(length(idx) > 0){
    print(rownames(biotype_perm_p)[i])
    cat("\n")
    print("P")
    print(colnames(biotype_perm_p)[-which(colnames(biotype_perm_p)=="protein_coding")][idx])
    cat("\n")
  }
  biotype_perm_fdr[i,] = p.adjust(biotype_perm_p[i,],method="fdr")
  idx = which(biotype_perm_fdr[i,-which(colnames(biotype_perm_fdr)=="protein_coding")] < 0.05)
  if(length(idx) > 0){
    print("FDR")
    print(colnames(biotype_perm_fdr)[-which(colnames(biotype_perm_fdr)=="protein_coding")][idx])
    cat("\n")
  }
}

library(reshape2)
biotype_perm_p_melt = melt(biotype_perm_p)
biotype_perm_fdr_melt = melt(biotype_perm_fdr)
colnames(biotype_perm_p_melt) = c("Module","Gene_Biotype","P")
biotype_perm_results = data.frame(cbind(biotype_perm_p_melt,biotype_perm_fdr_melt[,3]))
colnames(biotype_perm_results)[4] = "FDR"

### also test explicitly for less than expected by chance

biotype_perm_p_lessPC = matrix(NA,nrow=35,ncol=length(uniq_types))
rownames(biotype_perm_p_lessPC) = module_key$Module_Name[-1]
colnames(biotype_perm_p_lessPC) = uniq_types

for(i in c(1:35)){
  for(j in c(1:length(uniq_types))){
    mod = rownames(biotype_perm_p_lessPC)[i]
    type = colnames(biotype_perm_p_lessPC)[j]
    perm_dist_iter = perm_dist[[mod]][[type]]
    true_count = RNA_type_table_true[mod,type]
    less=length(which(perm_dist_iter <= true_count))
    p = (less)/10001
    biotype_perm_p_lessPC[i,j] = p
  }
}

biotype_perm_fdr_lessPC <- biotype_perm_p_lessPC

for(i in c(1:35)){
  idx = which(biotype_perm_p_lessPC[i,] < 0.05)
  if(length(idx) > 0){
    print(rownames(biotype_perm_p_lessPC)[i])
    cat("\n")
    print("P")
    print(colnames(biotype_perm_p_lessPC)[idx])
    cat("\n")
  }
  biotype_perm_fdr_lessPC[i,] = p.adjust(biotype_perm_p_lessPC[i,],method="fdr")
  idx = which(biotype_perm_fdr_lessPC[i,] < 0.05)
  if(length(idx) > 0){
    print("FDR")
    print(colnames(biotype_perm_fdr_lessPC)[idx])
    cat("\n")
  }
}

PC = biotype_perm_p_lessPC[,which(colnames(biotype_perm_p_lessPC)=="protein_coding")]
names(PC) = rownames(biotype_perm_p_lessPC)
PC[which(PC < 0.05)]
#   M9_magenta M18_lightgreen 
#0.00779922     0.03329667 

biotype_perm_p_melt = melt(biotype_perm_p_lessPC)
biotype_perm_fdr_melt = melt(biotype_perm_fdr_lessPC)
colnames(biotype_perm_p_melt) = c("Module","Gene_Biotype","P")
biotype_perm_results_less = data.frame(cbind(biotype_perm_p_melt,biotype_perm_fdr_melt[,3]))
colnames(biotype_perm_results_less)[4] = "FDR"

biotype_perm_results_more = biotype_perm_results

save(biotype_perm_results_more,biotype_perm_results_less,file="data_user/04_WGCNA/04_01_A_05_BiotypePermutation_SigResults.RData")

## plot with ggplot2 heatmap

pdf(file="plots/04_WGCNA/04_01_A_05_Biotype_Permutation.pdf",width=10,height=14)  

datPlot = biotype_perm_results_more

datPlot$log10p = -log10(as.numeric(datPlot$P) + 1e-5)
datPlot$star = rep(NA,nrow(datPlot))
datPlot$star[which(datPlot$P < 0.05)] = "*"
datPlot$Module = factor(datPlot$Module,levels=rev(unique(datPlot$Module)))

lab_size=22
text_size=22

biotype_p <- ggplot(datPlot, aes(Gene_Biotype, Module,fill=log10p)) + 
  geom_tile(aes(fill = log10p),colour = "grey50") + 
  scale_fill_gradient2(low = "dodgerblue",mid="white",high = "firebrick",limits=c(-5,5)) +
  geom_text(aes(label = star), color = "black", size = 10,vjust=0.75) +
  xlab("") +
  ylab("") +
  labs(fill="Signed\nlog10(P)") +
  theme_bw() +
  theme(legend.title = element_text(size = lab_size),
        legend.text = element_text(size = text_size),
        axis.text.x = element_text(size = text_size,angle=45,hjust=1),
        axis.text.y = element_text(size = text_size),
        legend.position = "left",
        legend.direction = "vertical")

print(biotype_p)

datPlot = biotype_perm_results_more

datPlot$log10p = -log10(as.numeric(datPlot$FDR) + 1e-5)
datPlot$star = rep(NA,nrow(datPlot))
datPlot$star[which(datPlot$FDR < 0.05)] = "*"
datPlot$Module = factor(datPlot$Module,levels=rev(unique(datPlot$Module)))

lab_size=22
text_size=22

biotype_fdr <- ggplot(datPlot, aes(Gene_Biotype, Module,fill=log10p)) + 
  geom_tile(aes(fill = log10p),colour = "grey50") + 
  scale_fill_gradient2(low = "dodgerblue",mid="white",high = "firebrick",limits=c(-5,5)) +
  geom_text(aes(label = star), color = "black", size = 10,vjust=0.75) +
  xlab("") +
  ylab("") +
  labs(fill="Signed\nlog10(FDR)") +
  theme_bw() +
  theme(legend.title = element_text(size = lab_size),
        legend.text = element_text(size = text_size),
        axis.text.x = element_text(size = text_size,angle=45,hjust=1),
        axis.text.y = element_text(size = text_size),
        legend.position = "left",
        legend.direction = "vertical")

print(biotype_fdr)

datPlot = biotype_perm_results_less

datPlot$log10p = log10(as.numeric(datPlot$P) + 1e-5)
datPlot$star = rep(NA,nrow(datPlot))
datPlot$star[which(datPlot$P < 0.05)] = "*"
datPlot$Module = factor(datPlot$Module,levels=rev(unique(datPlot$Module)))

lab_size=22
text_size=22

biotype_p <- ggplot(datPlot, aes(Gene_Biotype, Module,fill=log10p)) + 
  geom_tile(aes(fill = log10p),colour = "grey50") + 
  scale_fill_gradient2(low = "dodgerblue",mid="white",high = "firebrick",limits=c(-5,5)) +
  geom_text(aes(label = star), color = "black", size = 10,vjust=0.75) +
  xlab("") +
  ylab("") +
  labs(fill="Signed\nlog10(P)") +
  theme_bw() +
  theme(legend.title = element_text(size = lab_size),
        legend.text = element_text(size = text_size),
        axis.text.x = element_text(size = text_size,angle=45,hjust=1),
        axis.text.y = element_text(size = text_size),
        legend.position = "left",
        legend.direction = "vertical")

print(biotype_p)

datPlot = biotype_perm_results_less

datPlot$log10p = log10(as.numeric(datPlot$FDR) + 1e-5)
datPlot$star = rep(NA,nrow(datPlot))
datPlot$star[which(datPlot$FDR < 0.05)] = "*"
datPlot$Module = factor(datPlot$Module,levels=rev(unique(datPlot$Module)))

lab_size=22
text_size=22

biotype_fdr <- ggplot(datPlot, aes(Gene_Biotype, Module,fill=log10p)) + 
  geom_tile(aes(fill = log10p),colour = "grey50") + 
  scale_fill_gradient2(low = "dodgerblue",mid="white",high = "firebrick",limits=c(-5,5)) +
  geom_text(aes(label = star), color = "black", size = 10,vjust=0.75) +
  xlab("") +
  ylab("") +
  labs(fill="Signed\nlog10(FDR)") +
  theme_bw() +
  theme(legend.title = element_text(size = lab_size),
        legend.text = element_text(size = text_size),
        axis.text.x = element_text(size = text_size,angle=45,hjust=1),
        axis.text.y = element_text(size = text_size),
        legend.position = "left",
        legend.direction = "vertical")

print(biotype_fdr)

dev.off()


