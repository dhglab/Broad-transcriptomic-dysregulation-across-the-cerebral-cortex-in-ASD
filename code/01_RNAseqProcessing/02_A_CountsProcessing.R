##### 02_A_CountsProcessing.R
##### RNA-seq Counts Processing (Gene-Level)
##### May 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Transcriptomic-changes-are-pervasive-across-11-cortical-regions-in-ASD/",sep=""))

library(limma); library(edgeR); library(cqn); library(biomaRt)
library(WGCNA); library(devtools); library(paralell)
library(ggplot2); library(gridExtra); library(variancePartition)
library(doParallel); library(lmerTest)

## Install custom 'earth' R package to allow infinite genes as input
devtools::install(paste(wkdir,"code/earth-infGenes",sep=""))
library(earth)

##### Pipeline #####

### (1) Gene Filtering
### (2) Normalization
### (3) Outlier Removal
### (4) Prepare Metadata for Linear Model Building
### (5) Run CV MARS
### (6) Build Linear Model
### (7) Partition Variance Across Covariates
### (8) Obtain Regressed Gene Expression Data

##### (1) Gene Filtering #####

load("data_provided/01_RNAseqProcessing/02_A_01_RawData.RData")

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

save(datExpr,datMeta,datSeq,datSeq_numeric,rsem_gene_effLen,file="data_user/01_RNAseqProcessing/02_A_02_NormalizedGeneFilteredRNAseq.RData")

rm(list=ls())

##### (3) Outlier Removal #####

load("data_user/01_RNAseqProcessing/02_A_02_NormalizedGeneFilteredRNAseq.RData")

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

save(datExpr,datMeta,datSeq_model,datSeq_numeric,topPC,rsem_gene_effLen,file="data_user/01_RNAseqProcessing/02_A_03_NormalizedOutliersRemovedRNAseq.RData")

##### (4) Prepare Metadata for Linear Model Building #####

rm(list=ls())
load(file="data_user/01_RNAseqProcessing/02_A_03_NormalizedOutliersRemovedRNAseq.RData")

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
datMeta_model_org = datMeta_model
datMeta_model = datMeta_model[,-which(colnames(datMeta_model)=="subject")]
# remove subject, since we already know that we need to include subject in the model

save(datMeta_model,datSeq_model,datExpr,file="data_user/01_RNAseqProcessing/02_A_04_dat4MARS.RData")

datMeta_model = datMeta_model_org

save(datExpr,datMeta,datSeq_numeric,datMeta_model,datSeq_model,topPC,rsem_gene_effLen,
     file="data_user/01_RNAseqProcessing/02_A_04_AllDataPreMARS.RData")

##### (5) Run CV MARS #####

### Run EARTH/MARS.
### Use earth-inifite R package installation.
### It is recommended to run this section on a high-performance computing cluster.

runMARS <- function(gene.expr, covars, n.predictors=NULL, 
                    n.replicates=1, n.cores=1, allow.interaction=FALSE, 
                    batch.interactions=FALSE,allow.nonlinear=FALSE,cr=10) {
  
  # uses the packages `earth` to determine an appropriate linear model for expression, given
  # technical covariates.
  # Inputs
  #  gene.expr     - The gene expression matrix (genes x samples)
  #  covars        - The sample covariate matrix (samples x covariates)
  #  n.predictors  - The (maximum) number of predictors to use in the forward earth  model
  #  no.cross      - vector of covariates to exclude from cross terms
  #
  # Returns:
  #  earth         - the fitted `earth` object
  #  formula       - formula giving the right-hand-size of the covariate-correction LM
  #  terms         - the terms included
  #  term.imp      - importance for each term: this is the 80th percentile of max(beta)
  #  model.matrix  - the model matrix resulting from model.matrix(formula, covars)
  
  gene.expr = t(gene.expr)
  
  binary.terms = apply(covars,2,function(x) { 
    return(length(levels(as.factor(x)))==2 && min(x)==0 && max(x)==1)})
  square.terms = covars[,!binary.terms]^2
  colnames(square.terms) = paste0( colnames(square.terms), "^2")
  covars = cbind(covars, square.terms)
  
  if(allow.interaction==TRUE) {
    degree=2;
  } else {
    degree=1;
  }
  
  allowed.fx <- function(degree, pred, parents, namesx, first) {
    if(batch.interactions==FALSE & degree > 1) {
      bad.idx <- grep("BATCH", toupper(namesx))
      
      return(!(pred %in% bad.idx || which(parents!=0) %in% bad.idx))
    } else {
      return(TRUE)
    }
    
  }
  
  if(allow.nonlinear==TRUE) {
    linear=FALSE;
  } else {
    linear=TRUE;
  }
  
    out = mclapply(1:n.replicates, function(r) {
    e = earth(x=covars, y=gene.expr, trace=2, degree=degree,linpreds=linear, allowed=allowed.fx,nfold=cr)
    med= apply(abs(e$coefficients),1,median)
    beta = data.frame(row.names=1:length(med), var=as.character(names(med)), abs.beta=med, replicate=r, 
                      modfit = paste0("RSS=",e$rss, " RSQ=",e$rsq, " GCV=",e$gcv, " GRSQ=",e$grsq))
    all = list("beta"=beta,"e"=e)
    }, mc.cores=n.cores)
  
  return(out)
}

load("data_user/01_RNAseqProcessing/02_A_04_dat4MARS.RData")

datSeq=datSeq_model
rm(datSeq_model)

datMeta=datMeta_model
rm(datMeta_model)

meta_dat=data.frame(datMeta)
all_dat=data.frame(datSeq)

form_meta_dat=paste0(colnames(meta_dat)[c(1:ncol(meta_dat))],collapse=" + ")
form_meta_dat=paste0("~",form_meta_dat,collapse="")

form_all_dat=paste0(colnames(all_dat)[c(1:ncol(all_dat))],collapse=" + ")
form_all_dat=paste0("~0+",form_all_dat,collapse="")

Y = datExpr
rm(datExpr)

X.meta = data.frame(model.matrix(as.formula(form_meta_dat),data=meta_dat))

X.all = data.frame(model.matrix(as.formula(form_all_dat),data=all_dat))

Y = data.frame(Y[,match(rownames(X.meta), colnames(Y))])  ### I know that everything is matched already

## add meta to the two seq model matrices

X.all=data.frame(X.meta,X.all)

#Remove covariates with 0 variance

to_keep = !(apply(X.all,2,sd)==0) 
X.all = X.all[,to_keep]

output="data_user/01_RNAseqProcessing/"

if(TRUE | !file.exists(paste0(output,"/02_A_05_MARSOutput.RData")))  {
  
  e.all = runMARS(gene.expr = Y, covars = X.all, n.replicates = 1, allow.interaction = F,batch.interactions=F,allow.nonlinear=F)
  save(e.all, file=paste0(output,"/02_A_05_MARSOutput.RData"))
  
  beta_e_all=e.all[[1]]$beta
  coef_e_all=e.all[[1]]$e$coefficients
  
  save(beta_e_all,coef_e_all,file=paste0(output,"/02_A_05_MARSOutput_subset.RData"))
}

### Now, compile cross-validated MARS output

cv_coef_list=list()
cv_r2_list=list()

for(i in c(1:10)){
  
  cv_coef_list[[i]]=e.all[[1]]$e$cv.list[[i]]$coefficients
  cv_r2_list[[i]]=e.all[[1]]$e$cv.list[[i]]$rsq
}

cv_rsq_all=e.all[[1]]$e$cv.rsq.tab[,24837]

save(cv_coef_list,cv_r2_list,cv_rsq_all,file==paste0(output,"/02_A_05_MARSOutput_CVCompiled.RData"))

##### (6) Build Linear Model #####

rm(list=ls())

### It is recommended to use the provided CV MARS data to be consistent with the publication.
load("data_user/01_RNAseqProcessing/02_A_04_AllDataPreMARS.RData")
load("data_provided/01_RNAseqProcessing/02_A_05_MARSOutput_subset.RData")
load("data_provided/01_RNAseqProcessing/02_A_05_MARSOutput_CVCompiled.RData")

pdf(file="output/01_RNAseqProcessing/02_A_06_MARS_CV_plots.pdf",width=16,height=12)
for(i in c(1:10)){
  
  dat = rbind(cv_coef_list[[i]])
  med_dat=abs(apply(dat,1,median))
  plot_dat=data.frame("Coef"=rownames(dat),"Median_Beta"=med_dat)
  plot_dat=plot_dat[-grep("(Intercept)",plot_dat[,1]),]
  
  unq_rsq=signif(cv_r2_list[[i]],3)
  cv_rsq=signif(cv_rsq_all[[i]],3)
  
  g1 <- ggplot(plot_dat, aes(x=reorder(Coef, Median_Beta), y=Median_Beta)) + 
    geom_bar(stat="identity")  +
    coord_flip() + 
    geom_label(aes(y=Inf,x=-Inf,hjust=1,vjust=-0.1,label=paste0("R^2=",unq_rsq,"; CV R^2=",cv_rsq)),
               size=8) +
    xlab("Absolute Median Beta") +
    ylab("Coefficient") + 
    ggtitle(paste("Median Beta Across All Genes: Iter ",i)) +
    theme(plot.title = element_text(hjust = 0.5,size=28),
          axis.text.x =element_text(size=22),
          axis.title=element_text(size=26))
  
  grid.arrange(g1)
  
}

dev.off()

### Keep the terms that are in the best cross validation (iteration 6; CV R2 of 0.391, regular R2 of 0.401).

final_coefs=rownames(cv_coef_list[[6]])

### Add RIN, Subject, and Region to the MARS-selected model, since we want to evaluate these terms.

final_coefs=c(final_coefs,"RIN","Subject","Region")
final_coefs=final_coefs[-1]

### Assemble final model matrix.
### This includes making the DxReg term, to enable region-specific diagnosis effect contrasts.

datMeta_model$Region = gsub("-","_",datMeta_model$Region)
datMeta_model$DxReg = as.factor(paste(datMeta_model$Diagnosis,datMeta_model$Region,sep="_"))

### Evaluate selected coefficients with variancePartition.
## mod_meta_cat: for the model
## mod_meta_cat_vp: for VariancePartition

mod_meta_cat=data.frame("Subject"=datMeta_model$subject,"DxReg"=datMeta_model$DxReg,
                        "SeqBatch"=datMeta_model$Batch,
                        "Sex"=datMeta_model$Sex,
                        "Ancestry"=datMeta_model$Ancestry_Genotype)
mod_meta_cat_vp=data.frame("Subject"=datMeta_model$subject,
                           "Diagnosis"=datMeta_model$Diagnosis,
                           "Region"=datMeta_model$Region,
                           "SeqBatch"=datMeta_model$Batch,
                           "Sex"=datMeta_model$Sex,
                           "Ancestry"=datMeta_model$Ancestry_Genotype)
mod_meta_cont=data.frame("Age"=datMeta_model$Age,"Age_sqd"=(datMeta_model$Age)^2,"PMI"=datMeta_model$PMI,
                         "RIN"=datMeta_model$RIN,
                         "picard_gcbias.AT_DROPOUT"=datSeq_model$picard_gcbias.AT_DROPOUT,
                         "star.deletion_length"=datSeq_model$star.deletion_length,
                         "picard_rnaseq.PCT_INTERGENIC_BASES"=datSeq_model$picard_rnaseq.PCT_INTERGENIC_BASES,
                         "picard_insert.MEDIAN_INSERT_SIZE"=datSeq_model$picard_insert.MEDIAN_INSERT_SIZE,
                         "picard_alignment.PCT_CHIMERAS"=datSeq_model$picard_alignment.PCT_CHIMERAS,
                         "picard_alignment.PCT_PF_READS_ALIGNED"=datSeq_model$picard_alignment.PCT_PF_READS_ALIGNED,
                         "star.multimapped_percent"=datSeq_model$star.multimapped_percent,
                         "picard_rnaseq.MEDIAN_5PRIME_BIAS"=datSeq_model$picard_rnaseq.MEDIAN_5PRIME_BIAS,
                         "star.unmapped_other_percent"=datSeq_model$star.unmapped_other_percent,
                         "picard_rnaseq.PCT_USABLE_BASES"=datSeq_model$picard_rnaseq.PCT_USABLE_BASES,
                         "picard_alignment.PCT_CHIMERAS_sqd"=(datSeq_model$picard_alignment.PCT_CHIMERAS^2),
                         "star.uniquely_mapped_percent_sqd"=(datSeq_model$star.uniquely_mapped_percent^2))

cl <- makeCluster(4)
registerDoParallel(cl)

covs_cat=paste0("(1|",colnames(mod_meta_cat_vp),")",collapse=" + ")
covs=c(covs_cat,colnames(mod_meta_cont))
form = paste0(covs,collapse=" + ")
form = paste0("~",form,collapse=" ")

datMeta_model=data.frame(mod_meta_cat,mod_meta_cont)
# 21 total model covariates

datMeta_model_vp = data.frame(mod_meta_cat_vp,mod_meta_cont)

datMeta_model_keep = datMeta_model
datMeta_model = datMeta_model_vp

save(datExpr,form,datMeta_model_vp,file="data_user/01_RNAseqProcessing/02_A_06_VarPartData.RData")

datSeq = datSeq_model

save(datMeta_model,datMeta,datSeq,datSeq_numeric,
     topPC,datExpr,rsem_gene_effLen,file="data_user/01_RNAseqProcessing/02_A_06_AllProcessedData_wModelMatrix.RData")

##### (7) Partition Variance Across Covariates #####

### It is recommended to run this section on a high-performance computing cluster.

load("data_user/01_RNAseqProcessing/02_A_06_VarPartData.RData")

varPart <- fitExtractVarPartModel( datExpr, form, datMeta_model )

save(varPart,file="data_user/01_RNAseqProcessing/02_A_07_VarPartResults.RData")

vp <- sortCols( varPart )

pdf(file="output/01_RNAseqProcessing/02_A_07_varPartPlot.pdf",width=18,height=12)
plotPercentBars( vp[1:10,] )
plotVarPart(vp)
dev.off()

median(vp$Residuals) # 0.3054697
median(vp$Diagnosis) # 0.004216617

##### (8) Obtain Regressed Gene Expression Data #####

### It is recommended to run this section on a high-performance computing cluster.
### Genes can also be split up into blocks to run in parallel to generate the regressed expression data.

rm(list=ls())

load("data_user/01_RNAseqProcessing/02_A_06_AllProcessedData_wModelMatrix.RData")

### Make regressed gene expression data for visualization and WGCNA.
### Keep: Subject, Region, Sex, Diagnosis, Ancestry, Age, Age_sqd.

datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

form_dat=paste0(colnames(datMeta_model)[c(2:ncol(datMeta_model))],collapse=" + ")
form_dat=paste0("y ~ ",form_dat," + (1| Subject)",collapse="")

mod_remove=model.matrix(~SeqBatch,data=datMeta_model)[,-1]

d=c(1:nrow(datExpr))
max=2000
d1 <- split(d, ceiling(d/max))

for(j in c(1:length(names(d1)))){
  
  lmmod <- apply(as.matrix(datExpr[d1[[j]],]),1,function(y) lmer(as.formula(form_dat),data=datMeta_model))
  
  for (i in c(1:length(d1[[j]]))) {
    if (i%%100 == 0) {print(paste(i,"; Set ",j,sep=""))}
    
    summary = summary(lmmod[[i]])
    
    row_idx_batch <- grep("BatchASD",rownames(coef(summary)))
    
    datExpr.reg[d1[[j]][i],] <- datExpr[d1[[j]][i],] - mod_remove %*% coef(summary)[row_idx_batch,1] - as.matrix(datMeta_model[,c(8:dim(datMeta_model)[2])]) %*% coef(summary)[c(46:dim(coef(summary))[1]),1]
    
  }
}

### Get PCs for regressed gene expression.

norm <- t(scale(t(datExpr.reg),scale=F))
PC <- prcomp(norm,center = FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topPC_reg <- PC$rotation[,1:10] 
colnames(topPC_reg) <- paste("PC",c(1:10),"_",(signif(varexp[c(1:10)],2)*100),"%",sep="")

save(datExpr.reg,topPC_reg,file="data_user/01_RNAseqProcessing/02_A_08_RegressedExpression.RData")

### Make regressed gene expression data for TRI.
### Keep: Subject, Region, Diagnosis.

rm(list=ls())

load("data_user/01_RNAseqProcessing/02_A_06_AllProcessedData_wModelMatrix.RData")

datExpr.reg <- matrix(NA,nrow=nrow(datExpr),ncol=ncol(datExpr))
rownames(datExpr.reg) <- rownames(datExpr)
colnames(datExpr.reg) <- colnames(datExpr)

form_dat=paste0(colnames(datMeta_model)[c(2:ncol(datMeta_model))],collapse=" + ")
form_dat=paste0("y ~ ",form_dat," + (1| Subject)",collapse="")

mod_remove=model.matrix(~SeqBatch+Sex+Ancestry,data=datMeta_model)[,-1]

## do it in blocks to save memory

d=c(1:nrow(datExpr))
max=2000
d1 <- split(d, ceiling(d/max))

for(j in c(1:length(names(d1)))){
  
  lmmod <- apply(as.matrix(datExpr[d1[[j]],]),1,function(y) lmer(as.formula(form_dat),data=datMeta_model))
  
  for (i in c(1:length(d1[[j]]))) {
    if (i%%100 == 0) {print(paste(i,"; Set ",j,sep=""))}
    
    summary = summary(lmmod[[i]])
    
    row_idx_batch <- grep("BatchASD",rownames(coef(summary)))
    row_idx_sex <- grep("Sex",rownames(coef(summary)))
    row_idx_anc <- grep("Ancestry",rownames(coef(summary)))
    
    datExpr.reg[d1[[j]][i],] <- datExpr[d1[[j]][i],] - mod_remove %*% coef(summary)[c(row_idx_batch,row_idx_sex,row_idx_anc),1] - as.matrix(datMeta_model[,c(6:dim(datMeta_model)[2])]) %*% coef(summary)[c(44:dim(coef(summary))[1]),1]
    
  }
}

norm <- t(scale(t(datExpr.reg),scale=F)) ## center genes, NOT samples
PC <- prcomp(norm,center = FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topPC_reg <- PC$rotation[,1:10] 
colnames(topPC_reg) <- paste("PC",c(1:10),"_",(signif(varexp[c(1:10)],2)*100),"%",sep="")

save(datExpr.reg,topPC_reg,file="data_user/01_RNAseqProcessing/02_A_08_RegressedExpression_TRI.RData")
