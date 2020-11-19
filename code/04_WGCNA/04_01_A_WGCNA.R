##### 04_01_A_WGCNA.R
##### Perform WGCNA (Gene-Level)
##### May 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(WGCNA); library(gridExtra); library(ggplot2); library(lmtest)

load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8
rm(datExpr)
tdatExpr=t(datExpr.reg)

##### Pipeline #####

###(1) Choose parameters for rWGCNA
###(2) Perform rWGCNA
###(3) Identify module associations with covariates

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

