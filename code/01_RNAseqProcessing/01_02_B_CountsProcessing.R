##### 01_02_B_CountsProcessing.R
##### RNA-seq Counts Processing (Isoform-Level)
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

### Keep intersect of:
### (1) isoforms with cpm > 0.1 in 30% of samples
### (2) isoforms which intersect with gene-level analysis

### Get gene-level dataset genes

load("data_provided/02_DEGenesIsoforms/02_01_A_AllProcessedData_wModelMatrix.RData")
genes = substr(rownames(datExpr),1,15)
rm(datExpr,datMeta,datMeta_model,datSeq,datSeq_numeric,topPC,rsem_gene_effLen)

### Filter isoforms (cpm > 0.1 in 30% of samples)

load("data_provided/01_RNAseqProcessing/01_02_B_01_RawData.RData")

### 12 NCTL samples (Angelman Syndrome, CNV chr3p, and Epilepsy patients)
### Keep these samples only for the gene filtering step

cpm <- apply(rsem_tx,2, function(x) (x/sum(x))*1000000) 
keep = apply(cpm>0.1,1,sum) 
idx = which(keep > 0.3*dim(cpm)[2]) ## cpm > 0.1 in 30% of samples
cpm = cpm[idx,]

### Remove any remaining genes with an effective gene length <= 15 bp

length=rsem_transcript_effLen[match(rownames(cpm),names(rsem_transcript_effLen))] 

idx = which(length <= 15)
short_remove = rownames(cpm)[idx]
idx_remove = match(short_remove,rownames(cpm))
cpm = cpm[-idx_remove,]

### Filter current CPM list to match gene-level dataset

getinfo <- c("ensembl_transcript_id","ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gene_gc_content")

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="http://apr2018.archive.ensembl.org") 

## Use this ensembl version since, at the time the analysis was run,
## this was the most recent version.
## Using this version maximized the number of isoform matches with genes.

geneAnno <- getBM(attributes = getinfo,filters=c("ensembl_transcript_id"),
                   values=substr(rownames(cpm),1,15),mart=mart)

isos_keep = rownames(cpm)[which(substr(rownames(cpm),1,15) %in% geneAnno$ensembl_transcript_id)]
match_genes = geneAnno$ensembl_gene_id[match(substr(isos_keep,1,15),geneAnno$ensembl_transcript_id)]
isos_keep = isos_keep[which(match_genes %in% genes)]
length(isos_keep)

rsem_tx = rsem_tx[match(isos_keep,rownames(rsem_tx)),]
rsem_transcript_effLen=rsem_transcript_effLen[match(isos_keep,names(rsem_transcript_effLen))]
geneAnno = geneAnno[match(substr(isos_keep,1,15),geneAnno$ensembl_transcript_id),]

rm(cpm,length,keep,short_remove,idx,idx_remove,genes,isos_keep,match_genes,mart,getinfo)

dim(rsem_tx)
## 99819 isoforms,  854 samples

### Remove the NCTLs from the analysis

idx_remove=which(datMeta$Diagnosis=="NCTL")
datMeta=datMeta[-idx_remove,]
rsem_tx=rsem_tx[,-idx_remove]
datSeq=datSeq[-idx_remove,]
datSeq_numeric=datSeq_numeric[-idx_remove,]
datMeta$Diagnosis=droplevels(datMeta$Diagnosis)

##### (2) Normalization #####

counts = rsem_tx
rm(rsem_tx)

### Get CQN GC Content and Read Length Adjustment

geneAnno$effective_length = rsem_transcript_effLen[match(geneAnno[,1],substr(names(rsem_transcript_effLen),1,15))] 
rownames(geneAnno) = geneAnno[,1]

cqn.dat <- cqn(counts_cqn,lengths = as.numeric(geneAnno$effective_length), 
               x = as.numeric(geneAnno$percentage_gene_gc_content),
               lengthMethod=c("smooth"),sqn=FALSE) ## Run cqn with no quantile normalization

### limma-trend normalization, using the cqn offset values

dge2 <- DGEList(counts=counts_cqn)
dge2 <- calcNormFactors(dge2)
logCPM_offset<- cpm(dge2, log=TRUE, offset=cqn.dat$offset, prior.count=0.5) # prior count helps to dampen variance of low count genes

datExpr = logCPM_offset
rsem_transcript_effLen=rsem_transcript_effLen[match(rownames(datExpr),names(rsem_transcript_effLen))]

save(datExpr,datMeta,datSeq,datSeq_numeric,rsem_transcript_effLen,file="data_user/01_RNAseqProcessing/01_02_B_02_NormalizedGeneFilteredRNAseq.RData")

rm(list=ls())

##### (3) Outlier Removal #####

### Use the same outliers as the gene-level analysis to be consistent.

load("data_provided/02_DEGenesIsoforms/02_01_A_01_AllProcessedData_wModelMatrix.RData")
samples = colnames(datExpr)
rm(datExpr,datMeta,datMeta_model,datSeq,datSeq_numeric,topPC,rsem_gene_effLen)

load("data_user/01_RNAseqProcessing/01_02_B_02_NormalizedGeneFilteredRNAseq.RData")

idx_keep = match(samples,rownames(datMeta))
datExpr = datExpr[,idx_keep]
datMeta = datMeta[idx_keep,]
datSeq = datSeq[idx_keep,]
datSeq_numeric = datSeq_numeric[idx_keep,]

### Regenerate top expression PCs.

norm <- t(scale(t(datExpr),scale=F))
PC <- prcomp(norm,center = FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topPC <- PC$rotation[,1:15] ## these first 10 explain ~52% of the variance
colnames(topPC) <- paste("PC",c(1:15),"_",(signif(varexp[c(1:15)],2)*100),"%",sep="")

### Re-scale datSeq (now that outliers are removed)

idx_no_scale=grep("N",datSeq[1,])
datSeq_scaled = apply(datSeq_numeric[,-idx_no_scale],2,scale)
datSeq = as.data.frame(datSeq_scaled)
rownames(datSeq) = rownames(datSeq_numeric)
datSeq_model = datSeq

save(datExpr,datMeta,datSeq_model,datSeq_numeric,topPC,rsem_transcript_effLen,file="data_user/01_RNAseqProcessing/01_02_B_03_NormalizedOutliersRemovedRNAseq.RData")

##### (4) Prepare Metadata for Linear Model Building #####

rm(list=ls())
load(file="data_user/01_RNAseqProcessing/01_02_B_03_NormalizedOutliersRemovedRNAseq.RData")

### At this stage, use the datMeta and datMeta_model already processed for the gene-level analysis.

load("data_provided/01_RNAseqProcessing/01_02_B_04_InputCovariatesMARS.RData")

save(datMeta_model,datSeq_model,datExpr,file="data_user/01_RNAseqProcessing/01_02_B_04_dat4MARS.RData")

datMeta_model = datMeta_model_org
rm(datMeta_model_org)

save(datExpr,datMeta,datSeq_numeric,datMeta_model,datSeq_model,topPC,rsem_transcript_effLen,
     file="data_user/01_RNAseqProcessing/01_02_B_04_AllDataPreMARS.RData")

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

load("data_user/01_RNAseqProcessing/01_02_B_04_dat4MARS.RData")

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

if(TRUE | !file.exists(paste0(output,"/01_02_B_05_MARSOutput.RData")))  {
  
  e.all = runMARS(gene.expr = Y, covars = X.all, n.replicates = 1, allow.interaction = F,batch.interactions=F,allow.nonlinear=F)
  save(e.all, file=paste0(output,"/01_02_B_05_MARSOutput.RData"))
  
  beta_e_all=e.all[[1]]$beta
  coef_e_all=e.all[[1]]$e$coefficients
  
  save(beta_e_all,coef_e_all,file=paste0(output,"/01_02_B_05_MARSOutput_subset.RData"))
}

### Now, compile cross-validated MARS output

cv_coef_list=list()
cv_r2_list=list()

for(i in c(1:10)){
  
  cv_coef_list[[i]]=e.all[[1]]$e$cv.list[[i]]$coefficients
  cv_r2_list[[i]]=e.all[[1]]$e$cv.list[[i]]$rsq
}

cv_rsq_all=e.all[[1]]$e$cv.rsq.tab[,24837]

save(cv_coef_list,cv_r2_list,cv_rsq_all,file==paste0(output,"/01_02_B_05_MARSOutput_CVCompiled.RData"))

##### (6) Build Linear Model #####

rm(list=ls())

### It is recommended to use the provided CV MARS data to be consistent with the publication.
load("data_user/01_RNAseqProcessing/01_02_B_04_AllDataPreMARS.RData")
load("data_provided/01_RNAseqProcessing/01_02_B_05_MARSOutput_subset.RData")
load("data_provided/01_RNAseqProcessing/01_02_B_05_MARSOutput_CVCompiled.RData")

pdf(file="plots/01_RNAseqProcessing/01_02_B_06_MARS_CV_plots.pdf",width=16,height=12)
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
    ggtitle(paste("Median Beta Across All Isoforms: Iter ",i)) +
    theme(plot.title = element_text(hjust = 0.5,size=28),
          axis.text.x =element_text(size=22),
          axis.title=element_text(size=26))
  
  grid.arrange(g1)
  
}

dev.off()

### Keep the terms that are in the best cross validation (iteration 1; CV R2 of 0.218, regular R2 of 0.192).

final_coefs=rownames(cv_coef_list[[1]])

### Add some excluded terms to the MARS-selected model, since we want to evaluate these terms.

final_coefs=c(final_coefs,"RIN","Subject","Region","Diagnosis","Sex","Ancestry_Genotype","PMI")
final_coefs=final_coefs[-1]

## Use previous gene-level dataset to make final model matrix.

datExpr_tx = datExpr
topPC_tx = topPC

load("data_provided/02_DEGenesIsoforms/02_01_A_01_AllProcessedData_wModelMatrix.RData")

rm(datExpr,topPC)

datExpr = datExpr_tx
topPC = topPC_tx

rm(datExpr_tx,topPC_tx,rsem_gene_effLen)

## Get isoform datMeta_model.

colnames(datMeta_model)

datMeta_model = datMeta_model[,c(1:9)]

mod_meta_tech_cont=data.frame("picard_rnaseq.PCT_MRNA_BASES"=datSeq_model$picard_rnaseq.PCT_MRNA_BASES,
                              "picard_gcbias.AT_DROPOUT"=datSeq_model$picard_gcbias.AT_DROPOUT,
                              "picard_rnaseq.PCT_UTR_BASES"=datSeq_model$picard_rnaseq.PCT_UTR_BASES,
                              "star.multimapped_toomany_percent"=datSeq_model$star.multimapped_toomany_percent,                            
                              "picard_rnaseq.MEDIAN_CV_COVERAGE"=datSeq_model$picard_rnaseq.MEDIAN_CV_COVERAGE,
                              "picard_insert.MEDIAN_INSERT_SIZE"=datSeq_model$picard_insert.MEDIAN_INSERT_SIZE,
                              "picard_rnaseq.PCT_INTERGENIC_BASES"=datSeq_model$picard_rnaseq.PCT_INTERGENIC_BASES,
                              "picard_rnaseq.PF_BASES"=datSeq_model$picard_rnaseq.PF_BASES)

datMeta_model = data.frame(cbind(datMeta_model,mod_meta_tech_cont)) # 17 covariates

### Evaluate selected coefficients with variancePartition.
## mod_meta_cat_vp: for VariancePartition

mod_meta_cat_vp=data.frame("Subject"=datMeta_model$Subject,
                           "Diagnosis"=datMeta$Diagnosis,
                           "Region"=datMeta$region,
                           "SeqBatch"=datMeta_model$SeqBatch,
                           "Sex"=datMeta_model$Sex,
                           "Ancestry"=datMeta_model$Ancestry)
mod_meta_cont=data.frame("Age"=datMeta_model$Age,"Age_sqd"=(datMeta_model$Age)^2,"PMI"=datMeta_model$PMI,
                            "RIN"=datMeta_model$RIN,
                            "picard_rnaseq.PCT_MRNA_BASES"=datSeq_model$picard_rnaseq.PCT_MRNA_BASES,
                            "picard_gcbias.AT_DROPOUT"=datSeq_model$picard_gcbias.AT_DROPOUT,
                            "picard_rnaseq.PCT_UTR_BASES"=datSeq_model$picard_rnaseq.PCT_UTR_BASES,
                            "star.multimapped_toomany_percent"=datSeq_model$star.multimapped_toomany_percent,                            
                            "picard_rnaseq.MEDIAN_CV_COVERAGE"=datSeq_model$picard_rnaseq.MEDIAN_CV_COVERAGE,
                            "picard_insert.MEDIAN_INSERT_SIZE"=datSeq_model$picard_insert.MEDIAN_INSERT_SIZE,
                            "picard_rnaseq.PCT_INTERGENIC_BASES"=datSeq_model$picard_rnaseq.PCT_INTERGENIC_BASES,
                            "picard_rnaseq.PF_BASES"=datSeq_model$picard_rnaseq.PF_BASES)

cl <- makeCluster(4)
registerDoParallel(cl)

covs_cat=paste0("(1|",colnames(mod_meta_cat_vp),")",collapse=" + ")
covs=c(covs_cat,colnames(mod_meta_cont))
form = paste0(covs,collapse=" + ")
form = paste0("~",form,collapse=" ")

datMeta_model_vp = data.frame(mod_meta_cat_vp,mod_meta_cont)

save(datExpr,form,datMeta_model_vp,file="data_user/01_RNAseqProcessing/01_02_B_06_VarPartData.RData")

## Save final output for downstream analysis.

datSeq = datSeq_model

save(datMeta_model,datMeta,datSeq,datSeq_numeric,
     topPC,datExpr,rsem_transcript_effLen,file="data_user/01_RNAseqProcessing/01_02_B_06_AllProcessedData_wModelMatrix.RData")

##### (7) Partition Variance Across Covariates #####

### It is recommended to run this section on a high-performance computing cluster.

load("data_user/01_RNAseqProcessing/01_02_B_06_VarPartData.RData")

varPart <- fitExtractVarPartModel( datExpr, form, datMeta_model )

save(varPart,file="data_user/01_RNAseqProcessing/01_02_B_07_VarPartResults.RData")

vp <- sortCols( varPart )

pdf(file="plots/01_RNAseqProcessing/01_02_B_07_varPartPlot.pdf",width=18,height=12)
plotPercentBars( vp[1:10,] )
plotVarPart(vp)
dev.off()

median(vp$Residuals) # 0.6770287
median(vp$Diagnosis) # 0.0001032235

##### (8) Obtain Regressed Gene Expression Data #####

### It is recommended to run this section on a high-performance computing cluster.
### Genes can also be split up into blocks to run in parallel to generate the regressed expression data.

rm(list=ls())

load("data_user/01_RNAseqProcessing/01_02_B_06_AllProcessedData_wModelMatrix.RData")

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

save(datExpr.reg,topPC_reg,file="data_user/01_RNAseqProcessing/01_02_B_08_RegressedExpression.RData")

### Make regressed gene expression data for TRI.
### Keep: Subject, Region, Diagnosis.

rm(list=ls())

load("data_user/01_RNAseqProcessing/01_02_B_06_AllProcessedData_wModelMatrix.RData")

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

save(datExpr.reg,topPC_reg,file="data_user/01_RNAseqProcessing/01_02_B_08_RegressedExpression_TRI.RData")
