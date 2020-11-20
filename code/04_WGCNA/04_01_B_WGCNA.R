##### 04_01_B_WGCNA.R
##### Perform WGCNA (Isoform-Level)
##### November 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="/path/to/my/directory"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(WGCNA); library(gridExtra); library(ggplot2); library(lmtest); library(limma)

load("data_provided/04_WGCNA/04_01_B_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_B_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_B_RegressedExpression.RData") ### produced by 01_02_B_CountsProcessing.R, section 8
rm(datExpr)
tdatExpr=t(datExpr.reg)

##### Pipeline #####

###(1) Perform WGCNA
###(2) Compare isoform modules to gene modules
###(3) Identify module associations with covariates
###(4) Gene biotype permutation 

##### (1) Perform WGCNA #####

### Calculate soft power threshold and get WGCNA Modules

powers = c(seq(1,9,by=1),seq(10,30,by=2))
powerTable = pickSoftThreshold(data= tdatExpr, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = 5000)

### output:
### Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k.  max.k.
###1      1  0.00109  -4.27          0.812 5.11e+04  5.12e+04 52300.0
###2      2  0.06450 -17.60          0.833 2.59e+04  2.59e+04 27200.0
###3      3  0.19100 -16.50          0.926 1.32e+04  1.32e+04 14500.0
###4      4  0.44700 -16.90          0.971 6.74e+03  6.75e+03  8060.0
###5      5  0.62400 -14.30          0.969 3.48e+03  3.46e+03  4620.0
###6      6  0.75300 -11.60          0.961 1.81e+03  1.78e+03  2740.0
###7      7  0.80500  -9.18          0.944 9.45e+02  9.20e+02  1690.0
###8      8  0.86500  -7.14          0.937 4.99e+02  4.77e+02  1080.0
###9      9  0.91800  -6.31          0.956 2.66e+02  2.48e+02   751.0
###10    10  0.93500  -5.32          0.966 1.44e+02  1.29e+02   546.0
###11    12  0.96500  -3.76          0.986 4.39e+01  3.56e+01   324.0
###12    14  0.96300  -2.97          0.991 1.45e+01  9.95e+00   222.0
###13    16  0.96200  -2.43          0.995 5.36e+00  2.82e+00   165.0
###14    18  0.95400  -2.09          0.988 2.26e+00  8.15e-01   127.0
###15    20  0.95800  -1.87          0.990 1.10e+00  2.40e-01   102.0
###16    22  0.96100  -1.72          0.985 6.03e-01  7.16e-02    82.7
###17    24  0.96600  -1.63          0.988 3.65e-01  2.18e-02    68.3
###18    26  0.97700  -1.57          0.992 2.37e-01  6.74e-03    57.1
###19    28  0.98300  -1.52          0.993 1.61e-01  2.12e-03    48.2
###20    30  0.98500  -1.48          0.992 1.15e-01  6.79e-04    41.0

pdf("plots/04_WGCNA/04_01_B_01_CalculateSftPowerThresh_techRegressed.pdf",width=12)
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

### choose 10

### get TOM

### It is recommended to run this next section on a high-performance computing cluster.

pow=10

### use the same module cutting parameters derived from the gene level analysis

enableWGCNAThreads()
net = blockwiseModules(datExpr=tdatExpr, maxBlockSize=26000,networkType="signed",corType="bicor",  
                       power = pow, mergeCutHeight= 0.1, minModuleSize= 50, pamStage=FALSE, reassignThreshold=1e-6, 
                       saveTOMFileBase = "data_user/04_WGCNA/04_01_B_01_WGCNA_TOM", saveTOMs = TRUE, 
                       verbose = Inf, deepSplit=4)
save(net,file="data_user/04_WGCNA/04_01_B_01_WGCNA_network.RData")

### Cut Tree (get modules from TOM)
  
### For each block: load TOM, visualize dendrogram w/ modules

load("data_user/04_WGCNA/04_01_B_01_WGCNA_network.RData")
ds = 4; minModSize = 50; dthresh = 0.1; pam = FALSE

for(block in c(1:4)){
  
  load(paste("data_user/04_WGCNA/04_01_B_01_WGCNA_TOM-block.",block,".RData",sep=""))
  
  networks=list()
  networks$datExpr=t(tdatExpr[,net$blockGenes[[block]]])
  networks$tree = hclust(1-TOM, method="average")
  networks$cut = cutreeHybrid(dendro = networks$tree, pamStage=pam, minClusterSize= minModSize, cutHeight = 0.99999, deepSplit=ds, distM=as.matrix(1-TOM))
  networks$merged = mergeCloseModules(exprData= t(networks$datExpr), colors = networks$cut$labels, cutHeight=dthresh)
  networks$membership=data.frame(isoform=rownames(networks$datExpr),membership=networks$merged$colors)
  pdf(file=paste("plots/04_WGCNA/04_01_B_01_WGCNA_block",block,"_dendrogram.pdf",sep=""),height=8,width=11)
  plotDendroAndColors(networks$tree, colors=labels2colors(networks$merged$colors), dendroLabels = F)
  dev.off()
  save(file=paste("data_user/04_WGCNA/04_01_B_01_WGCNA_block",block,"_recut.RData",sep=""), networks) 
  
}


### Permute modules for robustness
  
enableWGCNAThreads()
allowWGCNAThreads()

nBlocks=4
maxUsedLabel=0
moduleSignif=data.frame()

for (i in 1:nBlocks){
  print(paste("working on block",i))
  lname=paste("data_user/04_WGCNA/04_01_B_01_WGCNA_block",i,"_recut.RData",sep="")
  load(lname)
  load(paste("data_user/04_WGCNA/04_01_B_01_WGCNA_TOM-block.",i,".RData",sep=""))
  TOM.mat=as.matrix(TOM)
  rm(TOM);gc()
  blockLength=dim(TOM.mat)[1]
  modules=networks$merged$colors
  for (j in 1:max(modules)){
    print(paste("working on module ",j))
    currentTOM=TOM.mat[modules==j,modules==j]
    currentSize=dim(currentTOM)[1]
    currentMean=mean(currentTOM)
    currentPermu=c()
    for (k in 1:5000){
      if ((k %% 100) == 0) {print(paste("iteration",k))}
      perm=sample(1:blockLength,currentSize,replace = F)
      permTOM=TOM.mat[perm,perm]
      permMean=mean(permTOM)
      currentPermu=c(currentPermu,permMean)
    }
    moduleSignif=rbind(moduleSignif,data.frame(module=j+maxUsedLabel,moduleMeanTOM=currentMean,permutMean=mean(currentPermu),permutSd=sd(currentPermu),p.value=1-pnorm(currentMean,mean = mean(currentPermu),sd = sd(currentPermu))))
  }
  maxUsedLabel=maxUsedLabel+max(modules)
  rm(networks,TOM.mat);gc()
}
write.table(moduleSignif,file = "data_user/04_WGCNA/04_01_B_01_module_significance_permutation.txt",quote = F,sep = "\t",row.names = F)
  
### final tree cut

enableWGCNAThreads()
allowWGCNAThreads()

nBlocks=4
maxUsedLabel=0
finalmember=data.frame()

for (block in 1:nBlocks){
  print(paste("Working on block ",block,".",sep=""))
  lname=paste("data_user/04_WGCNA/04_01_B_01_WGCNA_block",block,"_recut.RData",sep="")
  load(lname)
  current.member=networks$membership
  current.member$membership[current.member$membership > 0]=current.member$membership[current.member$membership > 0]+maxUsedLabel
  finalmember=rbind(finalmember,current.member)
  maxUsedLabel=max(finalmember$membership)
  rm(networks)
  gc()
}

print(dim(finalmember))

### Modules that did not pass the sig test (p < 0.01) in the previous step: 11,28,51,62,73,75
### Remove them

make_grey = which(finalmember[,2] %in% c(11,28,51,62,73,75))
finalmember[make_grey,2] = rep(0,length(make_grey))

# 77330 grey: 99819 - 77330 = 22489 transcripts in modules

#dthresh = 0.1 too many modules - try dthresh = 0.2
dthresh = 0.2
datExpr=tdatExpr
finalmember1=finalmember[match(colnames(datExpr),finalmember$isoform),]
merged=mergeCloseModules(exprData= datExpr, colors = as.vector(finalmember1$membership), cutHeight=dthresh)
# down to 62 from 76 with the higher merge threshold (61 only 5 more than PECDC). happy to continue with this one for now.
MEs = moduleEigengenes(datExpr, colors=merged$colors, softPower=10)
kMEtable = signedKME(datExpr, datME = MEs$eigengenes,corFnc = "bicor")

networks=list()
networks$merged=merged
networks$colors=merged$colors
networks$MEs=MEs
networks$kMEtable=kMEtable
save(file = "data_user/04_WGCNA/04_01_B_01_Full_Isoform_WGCNA_recut.RData",networks)
print(table(merged$colors))

##### (2) Identify module associations with covariates #####
  
load("data_user/04_WGCNA/04_01_B_01_Full_Isoform_WGCNA_recut.RData")  

### make a key - will name modules accordingly. important since some modules were merged, confusing to have mod 40 when only 35 mods

module_key_iso = data.frame("Color"=labels2colors(networks$merged$colors),
                            "Number_Original"=networks$merged$colors)

module_key_iso = module_key_iso[-which(duplicated(module_key_iso$Color)==TRUE),]
module_key_iso = module_key_iso[order(module_key_iso$Number_Original),]
module_key_iso$Number_Finalized = c(0:(nrow(module_key_iso)-1))
module_key_iso$Module_Name = paste("M",module_key_iso$Number_Finalized,"_",module_key_iso$Color,sep="")
module_key_iso$ME_Name = paste("ME",module_key_iso$Number_Finalized,"_",module_key_iso$Color,sep="")

labs <- as.numeric(gsub("ME","",colnames(networks$MEs$eigengenes)))
cols <- labels2colors(labs)
MEs = networks$MEs$eigengenes
colnames(MEs) <- paste("ME",module_key_iso$Number_Finalized[match(labs,module_key_iso$Number_Original)],"_",cols,sep = "")
rownames(MEs) <- colnames(datExpr.reg)

kME = networks$kMEtable
colors = labels2colors(networks$merged$colors)

modules = list()

for(col in module_key_iso$Color){
  modules[[col]]=rownames(datExpr.reg)[which(colors==col)]
}

MEs_iso = MEs
kME_iso = kME
colors_iso = colors
modules_iso = modules

save(module_key_iso,MEs_iso,kME_iso,colors_iso,modules_iso,file="data_user/04_WGCNA/04_01_B_02_isoformModules_filtered.RData")

### cluster iso MEs w/ the gene level MEs

load("data_user/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")

## OR: load("data_provided/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData")

# Calculate dissimilarity of module eigengenes
MEs_iso = MEs_iso[,-1]
MEs = MEs[,-1]
MEs_genes = MEs

colnames(MEs_iso) = paste("Isoform",colnames(MEs_iso),sep="_")
colnames(MEs_genes) = paste("Gene",colnames(MEs_genes),sep="_")

MEs = data.frame(cbind(MEs_genes,MEs_iso))

MEDiss = 1-bicor(MEs);

# Cluster module eigengenes

METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result

pdf(file="plots/04_WGCNA/04_01_B_02_Joint-ME-Tree.pdf",width=20,height=10)
plot(METree, main = "Clustering of module eigengenes",
   xlab = "", sub = "")
abline(h=0.5,col="blue")
abline(h=0.75,col="red")
dev.off()

ME_order = METree$order

MEs_all = MEs
MEs = MEs_iso

## 96 modules total: 35 gene, 61 isoform.

### Define 'iso-specific' module clusters

ME_groups = cutreeHybrid(dendro = METree, pamStage=FALSE,
                       minClusterSize= 1, cutHeight = 0.99999, deepSplit=4, distM=MEDiss)

pdf(file="plots/04_WGCNA/04_01_B_02_Joint-ME-Tree_wColors.pdf",width=20,height=10)

plotDendroAndColors(METree, colors=labels2colors(ME_groups$labels))

dev.off()

ME_groups_colors = ME_groups$labels
names(ME_groups_colors) = colnames(MEDiss)

iso_specific_modules = names(ME_groups_colors)[which(ME_groups_colors==0)]
iso_specific_modules = iso_specific_modules[-grep("Gene",iso_specific_modules)]
iso_specific_modules = c(iso_specific_modules,"Isoform_ME26_skyblue3") # this one is deserving too

for(i in c(1:max(ME_groups_colors))){
idx = grep("Gene",names(ME_groups_colors)[which(ME_groups_colors==i)])
if(length(idx) < 1){
  iso_specific_modules = c(iso_specific_modules,names(ME_groups_colors)[which(ME_groups_colors==i)])
}
}

## 37 isoform specific modules out of the original 61

iso_specific_modules

## Validate that the iso-conserved modules are conserved with respective gene level modules: fisher exact test

iso_conserved_modules = colnames(MEs)[-which(colnames(MEs) %in% iso_specific_modules)]

corr_gene_mods = list()

for(mod in iso_conserved_modules){
group = ME_groups_colors[which(names(ME_groups_colors)==mod)]
all_in_group = names(ME_groups_colors)[which(ME_groups_colors==group)]
gene_in_group = all_in_group[grep("Gene",all_in_group)]
corr_gene_mods[[mod]]$Corr_Gene = gene_in_group
}

source("code/04_WGCNA/fisher_overlap.R")
load("data_provided/04_WGCNA/gencode.v25lift37.annotation.gene_to_transcript.key.RData")

load("data_provided/04_WGCNA/04_01_B_02_WholeCortex_FullDGE.RData") ### extracted from 02_01_A_DEGenes.R, section 1

### get Fisher's Exact Test p-values

for(mod1 in iso_conserved_modules){
for(mod2 in corr_gene_mods[[mod1]]$Corr_Gene){
  f_out=ORA(substr(unique(key$gene_id[match(modules_iso[[which(names(modules_iso)==module_key_iso$Color[which(module_key_iso$ME_Name==gsub("Isoform_","",mod1))])]],key$transcript_id)]),1,15),
            substr(modules[[module_key$Color[which(module_key$ME_Name==gsub("Gene_","",mod2))]]],1,15),
            substr(unique(key$gene_id[match(rownames(datExpr.reg),key$transcript_id)]),1,15),
            substr(rownames(tt_ASD_Region),1,15))
  corr_gene_mods[[mod1]]$OR[which(corr_gene_mods[[mod1]]$Corr_Gene==mod2)]=as.numeric(f_out[1])
  corr_gene_mods[[mod1]]$P[which(corr_gene_mods[[mod1]]$Corr_Gene==mod2)]=as.numeric(f_out[2])
}
}

### compile list into table

iso_gene_conserved_mods = data.frame("IsoMod"=NA,"GeneMod"=NA,"OR"=NA,"P"=NA)

for(mod1 in iso_conserved_modules){
for(mod2 in corr_gene_mods[[mod1]]$Corr_Gene){
  iso_gene_conserved_mods = data.frame(rbind(iso_gene_conserved_mods,
                                             data.frame("IsoMod"=mod1,
                                                        "GeneMod"=mod2,
                                                        "OR"=corr_gene_mods[[mod1]]$OR[which(corr_gene_mods[[mod1]]$Corr_Gene==mod2)],
                                                        "P"=corr_gene_mods[[mod1]]$P[which(corr_gene_mods[[mod1]]$Corr_Gene==mod2)])))
}
}

iso_gene_conserved_mods = iso_gene_conserved_mods[-1,]
idx = which(iso_gene_conserved_mods$P > 0.01)
iso_gene_conserved_mods[idx,]

#                     IsoMod        GeneMod       OR         P
#34        Isoform_ME55_blue2 Gene_ME7_black 1.316680 0.5360104
#36 Isoform_ME61_navajowhite1 Gene_ME14_cyan 2.912641 0.2957798

length(modules_iso[["blue2"]]) #74
length(modules[["black"]]) # 253
length(modules_iso[["navajowhite1"]]) # 54
length(modules[["cyan"]]) #160

length(intersect(modules[["black"]],key$gene_id[match(modules_iso[["blue2"]],key$transcript_id)])) # 1
length(intersect(modules[["cyan"]],key$gene_id[match(modules_iso[["navajowhite1"]],key$transcript_id)])) # 1

### based on this, add these two iso mods to iso distinct list and remove from iso conserved list

iso_specific_modules = c(iso_specific_modules,"Isoform_ME55_blue2","Isoform_ME61_navajowhite1")
iso_conserved_modules = iso_conserved_modules[-which(iso_conserved_modules == "Isoform_ME61_navajowhite1" | iso_conserved_modules == "Isoform_ME55_blue2")]
iso_gene_conserved_mods = iso_gene_conserved_mods[-which(iso_gene_conserved_mods$IsoMod== "Isoform_ME61_navajowhite1" | iso_gene_conserved_mods$IsoMod== "Isoform_ME55_blue2")]

### finally: overlap with DTE/DGE intersects for the remaining isoform conserved modules?

load("data_provided/04_WGCNA/04_01_B_02_intersect_DE.RData") ### intersect/outersect of DE genes found in section 1 of 02_01_A_DEGenes.R and 02_01_B_DEGenes.R 

for(reg in names(gene_diff_intersect_logFC)){
overlap_p <- overlap_or <- rep(NA,nrow(iso_gene_conserved_mods))
for(i in c(1:nrow(iso_gene_conserved_mods))){
  f_out=ORA(substr(unique(key$gene_id[match(modules_iso[[which(names(modules_iso)==module_key_iso$Color[which(module_key_iso$ME_Name==gsub("Isoform_","",iso_gene_conserved_mods$IsoMod[i]))])]],key$transcript_id)]),1,15),
            substr(gene_diff_intersect_logFC[[reg]]$gene,1,15),
            substr(unique(key$gene_id[match(rownames(datExpr.reg),key$transcript_id)]),1,15),
            substr(rownames(tt_ASD_Region),1,15))
  overlap_or[i]=as.numeric(f_out[1])
  overlap_p[i]=as.numeric(f_out[2])
}
iso_gene_conserved_mods = data.frame(cbind(iso_gene_conserved_mods,overlap_p,overlap_or))
colnames(iso_gene_conserved_mods)= c(colnames(iso_gene_conserved_mods)[-c((ncol(iso_gene_conserved_mods)-1),ncol(iso_gene_conserved_mods))],paste(reg,"_int_p",sep=""),paste(reg,"_int_or",sep=""))
}

# not all, but a lot of these do overlap with these isoform modules.

colnames(iso_gene_conserved_mods)[c(3,4)] = c("Gene_v_Iso_OR","Gene_v_Iso_P")
colnames(iso_gene_conserved_mods)[c(5:ncol(iso_gene_conserved_mods))] = paste("DE_Intersect_",colnames(iso_gene_conserved_mods)[c(5:ncol(iso_gene_conserved_mods))],sep="")
colnames(iso_gene_conserved_mods) = gsub("int_","",colnames(iso_gene_conserved_mods))

save(iso_gene_conserved_mods,file="data_user/04_WGCNA/04_01_B_02_conserved_isoform_WGCNA_intDE_overlap.RData")

## Check and see if the iso-distinct modules are enriched in any gene-level modules and evaluate

iso_specific_modules # 39 now

## same thing, except now we are looking for an unsig p value in ALL gene modules

fisher_iso_distinct_p <- fisher_iso_distinct_or <- matrix(NA,nrow=length(iso_specific_modules),ncol=length(names(modules)))
rownames(fisher_iso_distinct_p) <- rownames(fisher_iso_distinct_or) <- iso_specific_modules
colnames(fisher_iso_distinct_p) <- colnames(fisher_iso_distinct_or) <- names(modules)

for(mod1 in iso_specific_modules){
for(mod2 in names(modules)){
  f_out=ORA(substr(unique(key$gene_id[match(modules_iso[[which(names(modules_iso)==module_key_iso$Color[which(module_key_iso$ME_Name==gsub("Isoform_","",mod1))])]],key$transcript_id)]),1,15),
            substr(modules[[mod2]],1,15),
            substr(unique(key$gene_id[match(rownames(datExpr.reg),key$transcript_id)]),1,15),
            substr(rownames(tt_ASD_Region),1,15))
  fisher_iso_distinct_or[mod1,mod2]=as.numeric(f_out[1])
  if(f_out[1] < 1){
    fisher_iso_distinct_p[mod1,mod2] = 1
  }else{
    fisher_iso_distinct_p[mod1,mod2]=as.numeric(f_out[2])
  }
}
}

## save any iso mod and gene mod pairs where p < 0.05

iso_distinct_mods_gene_overlaps = list()

for(mod1 in iso_specific_modules){
idx = which(fisher_iso_distinct_p[mod1,] < 0.05)
overlap_mods = module_key$Module_Name[match(colnames(fisher_iso_distinct_p)[idx],module_key$Color)]
overlap_mods = paste("Gene_",overlap_mods,sep="")
rem = grep("Gene_M0_grey",overlap_mods)
if(length(rem)==1){
  idx = idx[-rem]
  overlap_mods = overlap_mods[-rem]
}
if(length(idx) >= 1){
  iso_distinct_mods_gene_overlaps[[mod1]] = data.frame("Gene_Mod"=overlap_mods,
                                                       "P"=fisher_iso_distinct_p[mod1,idx],
                                                       "OR"=fisher_iso_distinct_or[mod1,idx])
}
}

save(iso_distinct_mods_gene_overlaps,file="data_user/04_WGCNA/04_01_B_02_distinct_isoMod_all_geneMod_overlap.RData")

### Now, save the isoform specific/distinct modules.

save(iso_specific_modules,file="data_user/04_WGCNA/04_01_B_02_isoMods_isoSpecific_noGeneCluster.RData")

##### (3) Identify module associations with covariates #####

### get whole cortex and region specific asd effects for isoform MEs, as well as effects of age and sex 

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

assoc_table_list_iso = list("WC_P"=assoc_table_p_wc,"WC_B"=assoc_table_beta_wc,
                            "ASDReg_P"=assoc_table_p_asdReg,"ASDReg_B"=assoc_table_beta_asdReg,
                            "Dup15Reg_P"=assoc_table_p_dup15Reg,"Dup15Reg_B"=assoc_table_beta_dup15Reg)

save(assoc_table_list_iso,file="data_user/04_WGCNA/04_01_B_03_isoME_assoc_table.RData")

##### (4) Gene biotype permutation #####

rm(list=ls())

allowWGCNAThreads()

wkdir="/path/to/my/directory"
load("data_provided/04_WGCNA/04_01_B_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_B_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_B_RegressedExpression.RData") ### produced by 01_02_B_CountsProcessing.R, section 8
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

load("data_user/04_WGCNA/04_01_B_02_isoformModules_filtered.RData")
### or: load("data_provided/04_WGCNA/04_01_B_02_isoformModules_filtered.RData")

### Count gene biotypes in modules
  
RNA_types = list()
all_genes = read.csv("main_datasets/B_IsoformLevel/AllIsoformsByModule_wDTE_wAnnotation.csv")

for(mod in module_key_iso$Module_Name[-1]){
  print(mod)
  tmp = all_genes[which(all_genes$WGCNA_module==mod),]
  RNA_types[[mod]]$length = nrow(tmp)
  RNA_types[[mod]]$types = table(tmp$gene_biotype)
}

uniq_types = NA

for(i in c(1:length(names(RNA_types)))){
  tmp_types = names(RNA_types[[i]]$types)
  uniq_types = c(uniq_types,tmp_types)
  uniq_types = uniq_types[!duplicated(uniq_types)]
}

uniq_types = uniq_types[-1]

RNA_type_table = matrix(NA,nrow=length(names(RNA_types)),ncol=length(uniq_types))
rownames(RNA_type_table) = module_key_iso$Module_Name[-1]
colnames(RNA_type_table) = uniq_types

for(i in c(1:length(uniq_types))){
  for(j in c(1:length(names(RNA_types)))){
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
save(uniq_types,RNA_type_table_true,file="data_user/04_WGCNA/04_01_B_04_RNA_biotypes_trueCount.RData")
  
### Run permutation analysis
### make random module assignments
  
for(iter in c(1:10000)){
  
  set.seed(as.numeric(iter))
  
  RNA_types = list()
  all_genes = read.csv("main_datasets/B_IsoformLevel/AllIsoformsByModule_wDTE_wAnnotation.csv")
  all_genes$WGCNA_module = sample(all_genes$WGCNA_module,replace=FALSE,size=nrow(all_genes))
  
  for(mod in module_key_iso$Module_Name[-1]){
    print(mod)
    tmp = all_genes[which(all_genes$WGCNA_module==mod),]
    RNA_types[[mod]]$length = nrow(tmp)
    RNA_types[[mod]]$types = table(tmp$gene_biotype)
  }
  
  RNA_type_table = matrix(NA,nrow=length(module_key_iso$Module_Name[-1]),ncol=length(uniq_types))
  rownames(RNA_type_table) = module_key_iso$Module_Name[-1]
  colnames(RNA_type_table) = uniq_types
  
  for(i in c(1:length(uniq_types))){
    for(j in c(1:length(module_key_iso$Module_Name[-1]))){
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
  save(RNA_type_table_perm,file=paste("data_user/04_WGCNA/04_01_B_04_permutations/RNA_biotypes_perm",iter,".RData",sep=""))
  
}

### Compile results
  
perm_dist = list()

for(i in c(1:10000)){
  load(paste("data_user/04_WGCNA/04_01_B_04_permutations/RNA_biotypes_perm",i,".RData",sep=""))
  for(mod in module_key_iso$Module_Name[-1]){
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

save(perm_dist,file="data_user/04_WGCNA/04_01_B_04_BiotypePermutationResults.RData")

### Examine results
  
biotype_perm_p = matrix(NA,nrow=length(module_key_iso$Module_Name[-1]),ncol=length(uniq_types))
rownames(biotype_perm_p) = module_key_iso$Module_Name[-1]
colnames(biotype_perm_p) = uniq_types

for(i in c(1:length(module_key_iso$Module_Name[-1]))){
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

for(i in c(1:length(module_key_iso$Module_Name[-1]))){
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

biotype_perm_p_lessPC = matrix(NA,nrow=length(module_key_iso$Module_Name[-1]),ncol=length(uniq_types))
rownames(biotype_perm_p_lessPC) = module_key_iso$Module_Name[-1]
colnames(biotype_perm_p_lessPC) = uniq_types

for(i in c(1:length(module_key_iso$Module_Name[-1]))){
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

for(i in c(1:length(module_key_iso$Module_Name[-1]))){
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

biotype_perm_p_melt = melt(biotype_perm_p_lessPC)
biotype_perm_fdr_melt = melt(biotype_perm_fdr_lessPC)
colnames(biotype_perm_p_melt) = c("Module","Gene_Biotype","P")
biotype_perm_results_less = data.frame(cbind(biotype_perm_p_melt,biotype_perm_fdr_melt[,3]))
colnames(biotype_perm_results_less)[4] = "FDR"

biotype_perm_results_more = biotype_perm_results

save(biotype_perm_results_more,biotype_perm_results_less,file="data_user/04_WGCNA/04_01_B_04_BiotypePermutation_SigResults.RData")

## plot with ggplot2 heatmap
 
pdf(file="plots/04_WGCNA/04_01_B_04_Biotype_Permutation.pdf",width=10,height=14)  

  datPlot = biotype_perm_results_more
  
  datPlot$log10p = -log10(as.numeric(datPlot$P) + 1e-5)
  datPlot$star = rep(NA,nrow(datPlot))
  datPlot$star[which(datPlot$P < 0.05)] = "*"
  datPlot$Module = factor(datPlot$Module,levels=rev(unique(datPlot$Module)))
  
  lab_size=16
  text_size=16
  
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
  
  #lab_size=22
  #text_size=22
  
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
  
  #lab_size=22
  #text_size=22
  
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
  
  #lab_size=22
  #text_size=22
  
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

