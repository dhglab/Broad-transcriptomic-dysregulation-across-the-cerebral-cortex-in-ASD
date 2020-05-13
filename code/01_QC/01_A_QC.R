##### 01_A_QC.R
##### RNA-seq Data Processing (Gene-Level)
##### May 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Transcriptomic-changes-are-pervasive-across-11-cortical-regions-in-ASD//",sep=""))

library(limma); library(edgeR); library(cqn); library(biomaRt)
library(WGCNA)

##### Pipeline #####

### (1) Gene Filtering
### (2) Normalization
### (3) Outlier Removal
### (4) Build Linear Model
### (5) Obtain Regressed Gene Expression Data

##### (1) Gene Filtering #####

load("data/01_QC/01_A_01_RawData.RData")

### 12 NCTL samples (Angelman Syndrome, CNV chr3p, and Epilepsy patients)
### Keep these samples only for the gene filtering step

### Keep genes with cpm > 0.1 in 30% of samples

cpm <- apply(rsem_gene,2, function(x) (x/sum(x))*1000000) 
keep = apply(cpm>0.1,1,sum) 
idx = which(keep > 0.3*dim(cpm)[2]) ## cpm > 0.1 in 30% of samples
cpm = cpm[idx,]

### 26,475 genes genes remain - about 44% of genes pass filter

### Remove any remaining genes with an effective gene length <= 15 bp

length=rsem_gene_effLen[match(rownames(cpm),names(rsem_gene_effLen))] 

idx = which(length <= 15)
short_remove = rownames(cpm)[idx]
idx_remove = match(short_remove,rownames(cpm))
cpm = cpm[-idx_remove,]
rsem_gene_effLen=rsem_gene_effLen[match(rownames(cpm),names(rsem_gene_effLen))]
rsem_gene = rsem_gene[match(rownames(cpm),rownames(rsem_gene)),]

rm(cpm,cpm_unfilt,length,keep,short_remove,idx,idx_remove)

dim(rsem_gene)
### 26364 genes, 854 samples

### Remove the NCTLs from the analysis

idx_remove=which(datMeta$Diagnosis=="NCTL")
datMeta=datMeta[-idx_remove,]
rsem_gene=rsem_gene[,-idx_remove]
datSeq=datSeq[-idx_remove,]
datSeq_numeric=datSeq_numeric[-idx_remove,]
datMeta$Diagnosis=droplevels(datMeta$Diagnosis)

##### (2) Normalization #####

counts = rsem_gene
rm(rsem_gene)

### Get CQN GC Content and Read Length Adjustment

getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gene_gc_content")

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="grch37.ensembl.org") 

geneAnno1 <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values=substr(rownames(counts),1,15),mart=mart)

geneAnno <- geneAnno1[match(substr(rownames(counts),1,15),geneAnno1[,1]),]
counts_cqn <- counts[!is.na(geneAnno[,1]),]
geneAnno <- geneAnno[!is.na(geneAnno[,1]),]

geneAnno$effective_length = rsem_gene_effLen[match(geneAnno[,1],substr(names(rsem_gene_effLen),1,15))] 
rownames(geneAnno) = geneAnno[,1]

cqn.dat <- cqn(counts_cqn,lengths = as.numeric(geneAnno$effective_length), x = as.numeric(geneAnno$percentage_gene_gc_content),
               lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with no quantile normalization

### limma-trend normalization, using the cqn offset values

dge2 <- DGEList(counts=counts_cqn)
dge2 <- calcNormFactors(dge2)
logCPM_offset<- cpm(dge2, log=TRUE, offset=cqn.dat$offset, prior.count=0.5) # prior count helps to dampen variance of low count genes

datExpr = logCPM_offset
rsem_gene_effLen=rsem_gene_effLen[match(rownames(datExpr),names(rsem_gene_effLen))]

save(datExpr,datMeta,datSeq,datSeq_numeric,rsem_gene_effLen,file="data/01_QC/01_A_02_NormalizedGeneFilteredRNAseq.RData")

rm(list=ls())

##### (3) Outlier Removal #####

load("data/01_QC/01_A_02_NormalizedGeneFilteredRNAseq.RData")

mod_batch=model.matrix(~0+datMeta$batch)
mod_reg=model.matrix(~0+datMeta$lobe)
colnames(mod_batch) = gsub("datMeta[$]batch","",colnames(mod_batch))
colnames(mod_reg)= gsub("datMeta[$]lobe","",colnames(mod_reg))

### PCA-based Outlier Removal
### Which samples have a z-score greater than 3 on the first 10 PCs? by batch.

outliers = list()

for(i in c(1:dim(mod_batch)[2])){
  for(j in c(1:dim(mod_reg)[2])){
    idx = which(mod_batch[,i]==1 & mod_reg[,j]==1)
    if(length(idx)==0){
      next
    }else{
      expr = datExpr[,idx]
      norm <- t(scale(t(expr),scale=F))
      PC <- prcomp(norm,center=FALSE)
      varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
      if(dim(PC$rotation)[2] >=10){
        topPC <- PC$rotation[,1:10] 
      }else{
        topPC <- PC$rotation[,1:dim(PC$rotation)[2]] 
      }
      for(z in c(1:dim(topPC)[2])){
        zscore=(topPC[,z]-mean(topPC[,z]))/sqrt(var(topPC[,z]))
        out = which(abs(zscore) > 3)
        if(length(out)!=0){
          out.flag = colnames(expr)[out]
          outliers[[paste(colnames(topPC)[z],
                          colnames(mod_reg)[j],
                          colnames(mod_batch)[i],sep="_")]] = out.flag
        }
      }
    }
  }
}

for(i in c(1:length(names(outliers)))){
  if(i==1){
    out.remove=outliers[[i]]
  }else{
    out.remove=c(out.remove,outliers[[i]])
  }
}

out.remove=out.remove[!duplicated(out.remove)]
out.remove_z3_pca=out.remove

### Connectivity-based Outlier Removal

outliers=list()

for(i in c(1:dim(mod_batch)[2])){
  for(j in c(1:dim(mod_reg)[2])){
    idx = which(mod_batch[,i]==1 & mod_reg[,j]==1)
    if(length(idx)==0){
      next
    }else{
      expr = datExpr[,idx]
      normadj <- (0.5+0.5*bicor(expr))^2 ## Calculate signed adjacency
      netsummary <- fundamentalNetworkConcepts(normadj) ## Calculate connectivity
      ku <- netsummary$Connectivity
      z.ku <- (ku-mean(ku))/sqrt(var(ku))
      idx <- z.ku < -2
      out.flag = colnames(expr)[idx]
      outliers[[paste(colnames(mod_batch)[i],colnames(mod_reg)[j],sep=" ")]] = out.flag
    }
  }
}

for(i in c(1:length(names(outliers)))){
  if(i==1){
    out.remove.conn=outliers[[i]]
  }else{
    out.remove.conn=c(out.remove.conn,outliers[[i]])
  }
}

### This method finds 50 outliers

int_outliers = intersect(out.remove_z3_pca,out.remove.conn)
length(int_outliers) ## 34 intersect

### Remove intersecting outliers

idx=which(rownames(datMeta) %in% int_outliers)
datExpr=datExpr[,-idx]
datMeta=datMeta[-idx,]
datSeq=datSeq[-idx,]
datSeq_numeric=datSeq_numeric[-idx,]

### Re-scale datSeq (now that outliers are removed)

idx_no_scale=grep("N",datSeq[1,])
datSeq_scaled = apply(datSeq_numeric[,-idx_no_scale],2,scale)
datSeq = as.data.frame(datSeq_scaled)
rownames(datSeq) = rownames(datSeq_numeric)
datSeq_model = datSeq

save(datExpr,datMeta,datSeq_model,datSeq_numeric,topPC,rsem_gene_effLen,file="data/01_QC/01_A_03_NormalizedOutliersRemovedRNAseq.RData")

##### (4) Build Linear Model #####

rm(list=ls())
load(file="data/01_QC/01_A_03_NormalizedOutliersRemovedRNAseq.RData")

### Format meta data for input into EARTH/MARS
datMeta$Brain_Bank_Source=gsub("London Brain","Oxford",datMeta$Brain_Bank_Source)
datMeta$Brain_Bank_Source=as.factor(datMeta$Brain_Bank_Source)

datMeta_model = datMeta[,c(2,3,5,6,8,9,14,33,34,35,36,38,40)]
datMeta_model = apply(datMeta_model,2,factor)
datMeta_model=data.frame(datMeta_model)

### For the 30 samples with no PMI, use the average of all other samples
idx_replace=which(is.na(datMeta_model$PMI)==TRUE)
PMI_average=mean(as.numeric(datMeta_model$PMI)[-idx_replace])
datMeta_model$PMI[idx_replace]=PMI_average

datMeta_model$PMI=as.numeric(datMeta_model$PMI)
datMeta_model$Age=as.numeric(datMeta_model$Age)
datMeta_model$RIN=as.numeric(datMeta_model$RIN)

datMeta_model$subject=as.factor(datMeta_model$subject)
datMeta_model$region=factor(datMeta_model$region, levels = c("BA17","BA20-37","BA24",
                                                             "BA3-1-2-5","BA38","BA39-40",
                                                             "BA4-6","BA41-42-22","BA44-45","BA7","BA9"))
datMeta_model$Diagnosis=factor(datMeta_model$Diagnosis,levels=c("CTL","ASD","Dup15q"))
datMeta_model$Ancestry_Genotype=factor(datMeta_model$Ancestry_Genotype,levels=c("EUR","AFR","ASN","HISPANIC","MIXED"))
datMeta_model$Sex=factor(datMeta_model$Sex,levels=c("M","F"))

colnames(datMeta_model)[c(2,8)]=c("Region","Batch")

datMeta_model$Brain_Bank_Source=as.factor(datMeta_model$Brain_Bank_Source)
datMeta_model$Batch=as.factor(datMeta_model$Batch)
datMeta_model$Read_Length=as.factor(datMeta$Read_Length)

### Scale the continuous datMeta_model variables

datMeta_model$Age=scale(datMeta_model$Age)
datMeta_model$PMI=scale(datMeta_model$PMI)
datMeta_model$RIN=scale(datMeta_model$RIN)

datMeta_model$Age=as.numeric(datMeta_model$Age)
datMeta_model$PMI=as.numeric(datMeta_model$PMI)
datMeta_model$RIN=as.numeric(datMeta_model$RIN)

datMeta_model$SeqMethod = as.factor(datMeta_model$SeqMethod)
datMeta_model$lobe = as.factor(datMeta_model$lobe)

### Now determine which datMeta_model covariates to remove going forward (obvious confounds): Read_Length, lobe, SeqMethod

datMeta_model = datMeta_model[-which(colnames(datMeta_model) %in% c("Read_Length","lobe","SeqMethod"))]

### Filter datSeq_model so that I am only evaluating terms that are not collinear with other covariates.
### For each covariate identify other covariates that have an adjusted R2 > 0.95.
### First, identify seq covariates (continuous) that are confounds for any meta covariates. 

allmat_meta = list()

for(i in c(1:dim(datMeta_model)[2])){
  tmp1 = datMeta_model[,i]
  allmat_meta[[colnames(datMeta_model)[i]]]=NA
  for(j in c(1:dim(datSeq_model)[2])){
    tmp2 = datSeq_model[,j]
    if(is.factor(tmp1)==TRUE & is.factor(tmp2)==FALSE){
      mod=summary(lm(tmp2~tmp1))
      r2=mod$adj.r.squared
      if(r2 >= 0.95){
        allmat_meta[[colnames(datMeta_model)[i]]] = c(allmat_meta[[colnames(datMeta_model)[i]]],colnames(datSeq_model)[j])
      }
    }else{
      mod=summary(lm(tmp1~tmp2))
      r2=mod$adj.r.squared
      if(r2 >= 0.95){
        allmat_meta[[colnames(datMeta_model)[i]]] = c(allmat_meta[[colnames(datMeta_model)[i]]],colnames(datSeq_model)[j])
      }
    }
  }
  allmat_meta[[colnames(datMeta_model)[i]]] = allmat_meta[[colnames(datMeta_model)[i]]][-1]
}

### Remove the seq covariates that are confounds with the meta covariates

set.seed(101788)

allmat_meta_org = allmat_meta

covs_keep = colnames(datSeq_model)
covs_iterate = names(allmat_meta_org)

for(cov in covs_iterate){
  print(cov)
  if(length(allmat_meta[[cov]])==0){
    next
  }else{
    keep0 = which(allmat_meta[[cov]] %in% covs_keep)
    remove = allmat_meta[[cov]][keep0]
    covs_keep = covs_keep[-which(covs_keep %in% remove)]
    print(length(covs_keep))
  }
}

### Remove 4 covariates

datSeq_model = datSeq_model[,which(colnames(datSeq_model) %in% covs_keep)]

### Now identify confounds among the remaining Seq covariates

allmat_seq = list()

for(i in c(1:dim(datSeq_model)[2])){
  tmp1 = datSeq_model[,i]
  allmat_seq[[colnames(datSeq_model)[i]]]=NA
  for(j in c(1:dim(datSeq_model)[2])){
    tmp2 = datSeq_model[,j]
    if(i==j){
      next
    }else{
      mod=summary(lm(tmp1~tmp2))
      r2=mod$adj.r.squared
      if(r2 >= 0.95){
        allmat_seq[[colnames(datSeq_model)[i]]] = c(allmat_seq[[colnames(datSeq_model)[i]]],colnames(datSeq_model)[j])
      }
    }
  }
  allmat_seq[[colnames(datSeq_model)[i]]] = allmat_seq[[colnames(datSeq_model)[i]]][-1]
}

## Iterate through each covariate, determine confounds. and if there are confounds
## For confounds, determine which of the covariates has the highest overall association with the top 5 expr PCs (adj R2).
## Then keep this covariate and remove the remaining from the loop.

## Get top expression PCs.

norm <- t(scale(t(datExpr),scale=F))
PC <- prcomp(norm,center = FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topPC <- PC$rotation[,1:15] ## these first 10 explain ~52% of the variance
colnames(topPC) <- paste("PC",c(1:15),"_",(signif(varexp[c(1:15)],2)*100),"%",sep="")

## Now begin filtering algorithm.

allmat_seq_org = allmat_seq
covs_keep = names(allmat_seq_org)
covs_iterate = names(allmat_seq_org)

for(cov in covs_iterate){
  print(cov)
  if((cov %in% covs_keep)==FALSE){
    next
  }
  if(length(allmat_seq[[cov]])==0){
    next
  }else{
    keep0 = which(allmat_seq[[cov]] %in% covs_keep)
    covs_test = c(cov,allmat_seq[[cov]][keep0])
    r2_mat = matrix(NA,nrow=5,ncol=length(covs_test))
    rownames(r2_mat)=colnames(topPC)[1:5]
    colnames(r2_mat)= covs_test
    for(i in c(1:5)){
      for(j in c(1:length(covs_test))){
        mod = summary(lm(topPC[,i] ~ datSeq_model[,which(colnames(datSeq_model)==covs_test[j])]))
        r2_mat[i,j]=mod$adj.r.squared
      }
    }
    r2_sum = apply(r2_mat,2,sum)
    max_val = max(r2_sum)
    max_cov = names(r2_sum)[which(r2_sum==max_val)]
    if(length(max_cov) > 1){
      keep = sample(max_cov,1)
    }else{
      keep = max_cov
    }
    allmat_seq[[keep]] = allmat_seq[[keep]][-c(1:length(allmat_seq[[keep]]))]
    remove = covs_test[-which(covs_test==keep)]
    covs_keep = covs_keep[-which(covs_keep %in% remove)]
    print(length(covs_keep))
  }
}

## Keeping 55 covariates (removing 11 covariates).
### Filter datSeq_model.

datSeq_model = datSeq_model[,which(colnames(datSeq_model) %in% covs_keep)]

### Run EARTH/MARS.
### Use earth-inifite R package installation.


##### (5) Obtain Regressed Gene Expression Data #####
