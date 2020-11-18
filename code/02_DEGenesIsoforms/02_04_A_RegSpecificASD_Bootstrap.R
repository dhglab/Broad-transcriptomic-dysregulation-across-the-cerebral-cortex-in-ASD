##### 02_04_A_RegSpecificASD_Bootstrap.R
##### ASD region-specific DGE bootstrap (Gene-Level)
##### May 2020, Jillian Haney

options(stringsAsFactors = FALSE)

library(ggplot2); library(reshape2); library(gridExtra); library(grid); library(limma); library(statmod)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

load("data_provided/02_DEGenesIsoforms/02_01_A_AllProcessedData_wModelMatrix.RData")
load("data_provided/02_DEGenesIsoforms/02_01_A_RegressedExpression.RData")
load("data_provided/02_DEGenesIsoforms/02_04_A_WholeCortex_DGE_logFC.RData")

boot=TRUE
compile=FALSE

if(boot==TRUE){
  
  ### Run 1,000 bootstraps (bootstrap=i) as below:
  ### Rscript 02_04_A_RegSpecificASD_Bootstrap.R i
  
  args = commandArgs(trailingOnly=TRUE)
  iter=as.numeric(args[1])
  
  if(iter==1){
    datMeta_model_perm = datMeta_model
    dir.create("data_user/02_DEGenesIsoforms/02_04_A_01_RegSpecificASD_bootstrap")
  }
  
  datMeta_model$Diagnosis = factor(datMeta$Diagnosis,levels=c("CTL","ASD","Dup15q"))
  datMeta_model$Region = as.factor(datMeta$region)
  dx_regs = unique(datMeta_model$DxReg)
  
  ### sample with replacement
  
  set.seed(as.numeric(iter))
  datMeta_model_perm = datMeta_model
  datMeta_model_perm$Sample_ID = rep(NA,nrow(datMeta_model_perm))
  rownames(datMeta_model_perm) <- rownames(datMeta_model) <- rownames(datMeta)
  
  for(group in dx_regs){
    idx = which(datMeta_model$DxReg==group)
    new_idx = sample(idx,replace=TRUE,size=length(idx))
    datMeta_model_perm[idx,]=datMeta_model[new_idx,]
    datMeta_model_perm$Sample_ID[idx] = rownames(datMeta_model)[new_idx]
  }
  
  rownames(datMeta_model_perm) = c(1:nrow(datMeta_model_perm))
  datExpr_match = datExpr[,match(datMeta_model_perm$Sample_ID,colnames(datExpr))]
  datExpr = datExpr_match
  rm(datMeta_model,datExpr_match,datExpr.reg)
  
  datMeta_model_perm$Region=gsub("-",".",datMeta_model_perm$Region)
  datMeta_model_perm$DxRegion = factor(paste0(datMeta_model_perm$Diagnosis, ".", datMeta_model_perm$Region))
  
  mod=model.matrix(~ 0 + DxRegion+SeqBatch+Sex+Ancestry,data=datMeta_model_perm)
  design=data.frame(mod,datMeta_model_perm[,c(6:21)])
  
  corfit <- duplicateCorrelation(datExpr,design,block=datMeta_model_perm$Subject)
  fit <- lmFit(datExpr, design,block=datMeta_model_perm$Subject,correlation=corfit$consensus)
  
  ### get region specific effects of ASD
  
  regs <- c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7","BA39.40","BA17","BA41.42.22","BA20.37","BA38")
  
  ASD_Region_ttables=list()
  
  for(j in c(1:length(regs))){
    print(j)
    reg = regs[j]
    reg_name = gsub("[.]","_",reg)
    form = paste("DxRegionASD.",reg," - DxRegionCTL.",reg,sep="")
    this_contrast = makeContrasts(contrast=form,levels=design)
    this_contrast <- data.frame(this_contrast[match(colnames(fit$coefficients),rownames(this_contrast)),])
    fit_reg = contrasts.fit(fit, this_contrast)
    fit2= eBayes(fit_reg,trend = T, robust = T)
    tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
    ASD_Region_ttables[[reg_name]]=tt_tmp
  }
  
  ### and dup15q
  
  this_contrast_dup15q = makeContrasts(contrast=(DxRegionDup15q.BA17 + DxRegionDup15q.BA20.37 + DxRegionDup15q.BA24 + DxRegionDup15q.BA3.1.2.5 +    
                                                   DxRegionDup15q.BA38 + DxRegionDup15q.BA39.40 + DxRegionDup15q.BA4.6 + DxRegionDup15q.BA41.42.22 +   
                                                   DxRegionDup15q.BA44.45 + DxRegionDup15q.BA7 + DxRegionDup15q.BA9 - DxRegionCTL.BA17 -         
                                                   DxRegionCTL.BA20.37 - DxRegionCTL.BA24 - DxRegionCTL.BA3.1.2.5 - DxRegionCTL.BA38 -          
                                                   DxRegionCTL.BA39.40 - DxRegionCTL.BA4.6 - DxRegionCTL.BA41.42.22 - DxRegionCTL.BA44.45 -      
                                                   DxRegionCTL.BA7 - DxRegionCTL.BA9)/11, levels=design)
  
  this_contrast_dup15q <- data.frame(this_contrast_dup15q[match(colnames(fit$coefficients),rownames(this_contrast_dup15q)),])
  
  fit_dup = contrasts.fit(fit, this_contrast_dup15q)
  fit2_dup= eBayes(fit_dup,trend = T, robust = T)
  tt_Dup15q= topTable(fit2_dup, coef=1, number=Inf, sort.by = 'none')
  
  ### get PC 1 slope
  
  wc_dge_genes= names(whole_cortex_logFC)
  
  pairsdat <- data.frame("Whole_Cortex"=whole_cortex_logFC,
                         "BA9"=ASD_Region_ttables[["BA9"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA44_45"=ASD_Region_ttables[["BA44_45"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA24"=ASD_Region_ttables[["BA24"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA4_6"=ASD_Region_ttables[["BA4_6"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA3_1_2_5"=ASD_Region_ttables[["BA3_1_2_5"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA7"=ASD_Region_ttables[["BA7"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA39_40"=ASD_Region_ttables[["BA39_40"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA38"=ASD_Region_ttables[["BA38"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA20_37"=ASD_Region_ttables[["BA20_37"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA41_42_22"=ASD_Region_ttables[["BA41_42_22"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "BA17"=ASD_Region_ttables[["BA17"]]$logFC[match(wc_dge_genes,rownames(ASD_Region_ttables[["BA9"]]))],
                         "Dup15q"=tt_Dup15q$logFC[match(wc_dge_genes,rownames(tt_Dup15q))],
                         "ID"=wc_dge_genes)
  plotDat= melt(pairsdat,id=c("ID","Whole_Cortex"))
  
  ##Calculate the slope of transcriptome overlap using principal components regression
  pcreg = function(ds1, ds2) {
    #Principal components regression to calculate slope 
    r = prcomp(~ds1+ds2)
    slope <- r$rotation[2,1] / r$rotation[1,1]
    intercept <- r$center[2] - slope*r$center[1]
    return(list(slope,intercept))
  }
  
  slope_boot = rep(NA,12)
  names(slope_boot) = c("BA9","BA44_45","BA24","BA4_6","BA3_1_2_5","BA7","BA39_40",
                        "BA38","BA20_37","BA41_42_22","BA17","Dup15q")
  for(reg in names(slope_boot)){
    slope_boot[reg] = pcreg(pairsdat$Whole_Cortex,pairsdat[[which(names(pairsdat)==reg)]])[[1]]
  }
  
  save(ASD_Region_ttables,tt_Dup15q,slope_boot,file=paste("data_user/02_DEGenesIsoforms/02_04_A_01_RegSpecificASD_bootstrap/ASDPan_DGE_bootstrap_",iter,".RData",sep=""))
  
}

if(compile==TRUE){
  
  ## get distribution for each region and dup15q
  
  slope_dist = list()
  regs = c("BA9","BA44_45","BA24","BA4_6","BA3_1_2_5","BA7","BA39_40",
           "BA38","BA20_37","BA41_42_22","BA17","Dup15q")
  
  for(i in c(1:1000)){
    if (i%%100 == 0) {print(paste(i,"/",length(list.files("data_user/02_DEGenesIsoforms/02_04_A_01_RegSpecificASD_bootstrap/")),sep=""))}
    load(paste("data_user/02_DEGenesIsoforms/02_04_A_01_RegSpecificASD_bootstrap/ASDPan_DGE_bootstrap_",i,".RData",sep=""))
    for(reg in regs){
      if(i==1){
        slope_dist[[reg]] = slope_boot[reg]
      }else{
        slope_dist[[reg]] = c(slope_boot[reg],slope_dist[[reg]])
      }
    }
  }
  
  save(slope_dist,file="data_user/02_DEGenesIsoforms/02_04_A_02_RegSpecificASD_bootstrap_results.RData")
  
}
