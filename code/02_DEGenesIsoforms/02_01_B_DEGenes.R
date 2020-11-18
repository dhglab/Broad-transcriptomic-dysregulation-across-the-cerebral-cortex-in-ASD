##### 02_01_B_DEGenes.R
##### Identify DE Isoforms (Isoform-Level)
##### November 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(limma); library(statmod); library(WGCNA)

##### Pipeline #####

###(1) Identify DE genes for covariates

##### (1) Identify DE isoforms for covariates #####

load("data_provided/02_DEGenesIsoforms/02_01_B_AllProcessedData_wModelMatrix.RData")

mod=model.matrix(~ 0 + DxReg+SeqBatch+Sex+Ancestry,data=datMeta_model)
design=data.frame(mod,datMeta_model[,c(6:(dim(datMeta_model)[2]))])

### It is recommended to run this next section on a high-performance computing cluster.

corfit <- duplicateCorrelation(datExpr,design,block=datMeta_model$Subject)
fit <- lmFit(datExpr,design,block=datMeta_model$Subject,correlation=corfit$consensus)

### The remaining steps do not require extensive computing resources.

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

this_contrast_asd <- data.frame(this_contrast_asd[match(colnames(fit$coefficients),rownames(this_contrast_asd)),])
this_contrast_dup15q <- data.frame(this_contrast_dup15q[match(colnames(fit$coefficients),rownames(this_contrast_dup15q)),])

### Get the whole cortex and region-specific contrasts

## Make tables for association of genes with biological traits
## First: whole cortex

assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,nrow(datExpr),5))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-  c("SexF","Age","Age_sqd","Diagnosis_ASD","Diagnosis_Dup15q")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  rownames(datExpr)

coefs = colnames(fit$coefficients)

for(j in c(1:ncol(assoc_table_p))){
  print(j)
  if(length(grep("Diagnosis",colnames(assoc_table_p)[j])) > 0){
    if(j==4){
      fit_asd = contrasts.fit(fit, this_contrast_asd)
      fit2_asd= eBayes(fit_asd,trend = T, robust = T)
      tt_ASD_Region= topTable(fit2_asd, coef=1, number=Inf, sort.by = 'none')
      write.csv(tt_ASD_Region,file="data_user/02_DEGenesIsoforms/02_01_A_01_ttable_ASD_WholeCortex.csv")
      for(gene in rownames(datExpr)){
        assoc_table_p[gene,j] = signif(tt_ASD_Region$adj.P.Val[which(rownames(tt_ASD_Region)==gene)],3)
        assoc_table_beta[gene,j] = signif(tt_ASD_Region$logFC[which(rownames(tt_ASD_Region)==gene)],3)
      }
    }else if(j==5){
      fit_dup = contrasts.fit(fit, this_contrast_dup15q)
      fit2_dup= eBayes(fit_dup,trend = T, robust = T)
      tt_Dup15q_Region= topTable(fit2_dup, coef=1, number=Inf, sort.by = 'none')
      write.csv(tt_Dup15q_Region,file="data_user/02_DEGenesIsoforms/02_01_A_01_ttable_dup15q_WholeCortex.csv")
      for(gene in rownames(datExpr)){
        assoc_table_p[gene,j] = signif(tt_Dup15q_Region$adj.P.Val[which(rownames(tt_Dup15q_Region)==gene)],3)
        assoc_table_beta[gene,j] = signif(tt_Dup15q_Region$logFC[which(rownames(tt_Dup15q_Region)==gene)],3)
      }    
    }
  }else{
    idx=grep(paste("\\b",colnames(assoc_table_p)[j],"\\b",sep=""),coefs)
    fit2= eBayes(fit,trend = T, robust = T)
    tt_tmp = topTable(fit2,coef=idx,number=Inf, sort.by = 'none')
    for(gene in rownames(datExpr)){
      assoc_table_p[gene,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_tmp)==gene)],3)
      assoc_table_beta[gene,j] = signif(tt_tmp$logFC[which(rownames(tt_tmp)==gene)],3)
    }
  }
}

assoc_table_p_wc = assoc_table_p  # contains whole cortex p-values for: ASD, sex, age, age_sqd
assoc_table_beta_wc = assoc_table_beta  # contains whole cortex effects/betas for: ASD, sex, age, age_sqd

save(assoc_table_p_wc,assoc_table_beta_wc,
     tt_ASD_Region,tt_Dup15q_Region,file="data_user/02_DEGenesIsoforms/02_01_B_01_WholeCortexFullDGE.RData")

## Next: get region-specific effects (one table for ASD, and one for dup15q)

## ASD

assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,nrow(datExpr),11))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-  c("BA9","BA24","BA44_45","BA4_6","BA3_1_2_5","BA7",
                                                            "BA39_40","BA17","BA41_42_22","BA20_37","BA38")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  rownames(datExpr)

ASD_Region_ttables=list()

for(j in c(1:ncol(assoc_table_p))){
  print(j)
  reg = colnames(assoc_table_p)[j]
  form = paste("DxRegASD_",reg," - DxRegCTL_",reg,sep="")
  this_contrast = makeContrasts(contrast=form,levels=design)
  this_contrast <- data.frame(this_contrast[match(colnames(fit$coefficients),rownames(this_contrast)),])
  fit_reg = contrasts.fit(fit, this_contrast)
  fit2= eBayes(fit_reg,trend = T, robust = T)
  tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
  for(gene in rownames(datExpr)){
    assoc_table_p[gene,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_tmp)==gene)],3)
    assoc_table_beta[gene,j] = signif(tt_tmp$logFC[which(rownames(tt_tmp)==gene)],3)
  }
  ASD_Region_ttables[[reg]]=tt_tmp
}

assoc_table_p_asdReg = assoc_table_p   # contains whole cortex p-values for all ASD region-specific contrasts
assoc_table_beta_asdReg = assoc_table_beta  # contains whole cortex effects/betas for all ASD region-specific contrasts

save(ASD_Region_ttables,assoc_table_p_asdReg,
     assoc_table_beta_asdReg,file="data_user/02_DEGenesIsoforms/02_01_B_01_ASD_RegionSpecific_FullDGE.RData")

for(reg in names(ASD_Region_ttables)){
  reg_name=gsub("[.]","-",reg)
  write.csv(ASD_Region_ttables[[reg]],paste("data_user/02_DEGenesIsoforms/02_01_B_01_ASD_",reg_name,"_DGE.csv",sep=""))
}

## dup15q

assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,nrow(datExpr),11))
colnames(assoc_table_p) <- colnames(assoc_table_beta) <-  c("BA9","BA24","BA44_45","BA4_6","BA3_1_2_5","BA7",
                                                            "BA39_40","BA17","BA41_42_22","BA20_37","BA38")
rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  rownames(datExpr)

Dup15q_Region_ttables=list()

for(j in c(1:ncol(assoc_table_p))){
  print(j)
  reg = colnames(assoc_table_p)[j]
  form = paste("DxRegDup15q_",reg," - DxRegCTL_",reg,sep="")
  this_contrast = makeContrasts(contrast=form,levels=design)
  this_contrast <- data.frame(this_contrast[match(colnames(fit$coefficients),rownames(this_contrast)),])
  fit_reg = contrasts.fit(fit, this_contrast)
  fit2= eBayes(fit_reg,trend = T, robust = T)
  tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
  for(gene in rownames(datExpr)){
    assoc_table_p[gene,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_tmp)==gene)],3)
    assoc_table_beta[gene,j] = signif(tt_tmp$logFC[which(rownames(tt_tmp)==gene)],3)
  }
  Dup15q_Region_ttables[[reg]]=tt_tmp
}

assoc_table_p_dup15Reg = assoc_table_p   # contains whole cortex p-values for all dup15q region-specific contrasts
assoc_table_beta_dup15Reg = assoc_table_beta  # contains whole cortex effects/betas for all dup15q region-specific contrasts

save(Dup15q_Region_ttables,assoc_table_p_dup15Reg,
     assoc_table_beta_dup15Reg,file="data_user/02_DEGenesIsoforms/02_01_B_01_dup15q_RegionSpecific_FullDGE.RData")

for(reg in names(Dup15q_Region_ttables)){
  reg_name=gsub("[.]","-",reg)
  write.csv(Dup15q_Region_ttables[[reg]],paste("data_user/02_DEGenesIsoforms/02_01_B_01_dup15q_",reg_name,"_DGE.csv",reg_name,".csv",sep=""))
}