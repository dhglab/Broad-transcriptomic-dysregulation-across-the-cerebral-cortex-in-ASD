##### 02_02_A_Parikshak_RegSpecificASD_Permutation.R
##### Parikshak et al. module region-specific permutation (Gene-Level)
##### Determine if region-specific ASD effect is significantly different from the whole cortex ASD effect.
##### May 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Transcriptomic-changes-are-pervasive-across-11-cortical-regions-in-ASD/",sep=""))

library(WGCNA); library(ggplot2); library(reshape2); library(limma)

PERM=TRUE
COMPILE=FALSE
SIGTEST=FALSE

##### Pipeline #####

### (1) Run Permutation
### (2) Compile Permutation Results
### (3) Test Severity of Regional Dysregulation

##### (1) Run Permutation #####

if(PERM==TRUE){
  
  load("data_provided/02_DEGenesIsoforms/02_01_A_AllProcessedData_wModelMatrix.RData")
  load("data_provided/02_DEGenesIsoforms/02_01_A_RegressedExpression.RData")
  
  ### Run 10,001 permutations (permutation=i) as below:
  ### Rscript 02_02_A_Parikshak_RegSpecificASD_Permutation.R i
  
  args = commandArgs(trailingOnly=TRUE)
  iter=as.numeric(args[1])
  
  datMeta_model$Diagnosis = factor(datMeta$Diagnosis,levels=c("CTL","ASD","Dup15q"))
  datMeta_model$Region = as.factor(datMeta$region)
  
  ### permute regional identity
  
  if(iter==1){
    datMeta_model_perm = datMeta_model
    dir.create("data_user/02_DEGenesIsoforms/02_02_A_01_parikshak_permutations")
  }else{
    set.seed(as.numeric(iter-1))
    datMeta_model_perm = datMeta_model
    for(sub in unique(datMeta_model$Subject)){
      idx = which(datMeta_model$Subject==sub)
      datMeta_model_perm$Region[idx]=sample(datMeta_model_perm$Region[idx])
    }
  }
  
  rm(datMeta_model)
  
  ### Parikshak et al. modules
  
  parik_mods = read.csv("data_provided/02_DEGenesIsoforms/02_01_A_ParikshakModules.csv")
  
  # extract sig_mods
  
  sig_mods = c(4,9,10,16,19,20)
  keep = which(parik_mods$WGCNA.Module.Label %in% sig_mods)
  parik_mods = parik_mods[keep,]
  
  int_genes = intersect(parik_mods$ENSEMBL.ID,substr(rownames(datExpr.reg),1,15))
  parik_mods = parik_mods[match(int_genes,parik_mods$ENSEMBL.ID),]
  datExpr.reg_match = datExpr.reg[match(int_genes,substr(rownames(datExpr.reg),1,15)),]
  
  parik_mes = moduleEigengenes(t(datExpr.reg_match), colors=parik_mods$WGCNA.Module.Label, softPower=8)
  parik_mes = parik_mes$eigengenes
  colnames(parik_mes) = gsub("E","",colnames(parik_mes))
  
  datMeta_model_perm$Region=gsub("-",".",datMeta_model_perm$Region)
  datMeta_model_perm$DxRegion = factor(paste0(datMeta_model_perm$Diagnosis, ".", datMeta_model_perm$Region))
  design = model.matrix(~ 0 + DxRegion + Sex + Ancestry + Age + Age_sqd, data = datMeta_model_perm)
  
  corfit= duplicateCorrelation(t(parik_mes),design,block=datMeta_model_perm$Subject)
  lm = lmFit(t(parik_mes), design,block=datMeta_model_perm$Subject, correlation = corfit$consensus)
  
  this_contrast_asd = makeContrasts(contrast=(DxRegionASD.BA17 + DxRegionASD.BA20.37 + DxRegionASD.BA24 + DxRegionASD.BA3.1.2.5 +
                                                DxRegionASD.BA38 + DxRegionASD.BA39.40 + DxRegionASD.BA4.6 + DxRegionASD.BA41.42.22 +
                                                DxRegionASD.BA44.45 + DxRegionASD.BA7 + DxRegionASD.BA9 - DxRegionCTL.BA17 -
                                                DxRegionCTL.BA20.37 - DxRegionCTL.BA24 - DxRegionCTL.BA3.1.2.5 - DxRegionCTL.BA38 -
                                                DxRegionCTL.BA39.40 - DxRegionCTL.BA4.6 - DxRegionCTL.BA41.42.22 - DxRegionCTL.BA44.45 -
                                                DxRegionCTL.BA7 - DxRegionCTL.BA9)/11, levels=design)
  
  ### Make tables for association of MEs with all traits
  
  ### just get region specific effects of ASD
  
  assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,ncol(parik_mes),11))
  colnames(assoc_table_p) <- colnames(assoc_table_beta) <- c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                                                             "BA39.40","BA17","BA41.42.22","BA20.37","BA38")
  rownames(assoc_table_p) <- rownames(assoc_table_beta) <-  colnames(parik_mes)
  
  for(j in c(1:ncol(assoc_table_p))){
    reg = colnames(assoc_table_p)[j]
    form = paste("DxRegionASD.",reg," - DxRegionCTL.",reg,sep="")
    this_contrast = makeContrasts(contrast=form,levels=design)
    fit = contrasts.fit(lm, this_contrast)
    fit2= eBayes(fit,trend = T, robust = T)
    tt_tmp= topTable(fit2, coef=1, number=Inf, sort.by = 'none')
    for(me in colnames(parik_mes)){
      assoc_table_p[me,j] = signif(tt_tmp$adj.P.Val[which(rownames(tt_tmp)==me)],3)
      assoc_table_beta[me,j] = signif(tt_tmp$logFC[which(rownames(tt_tmp)==me)],3)
    }
  }
  
  save(assoc_table_p,assoc_table_beta,
       file=paste("data_user/02_DEGenesIsoforms/02_02_A_01_parikshak_permutations/Parikshak_ME_permutation_",iter,".RData",sep=""))
  
}


##### (2) Compile Permutation Results #####

if(COMPILE==TRUE){
  
  ## iter=1 is the real result
  
  me_reg_effect_list = list()
  
  for(i in c(1:10001)){
    if (i%%100 == 0) {print(paste(i,"/",length(list.files("data_user/02_DEGenesIsoforms/02_02_A_01_parikshak_permutations")),sep=""))}
    load(paste("data_user/02_DEGenesIsoforms/02_02_A_01_parikshak_permutations/Parikshak_ME_permutation_",i,".RData",sep=""))
    for(me in c("M4","M9","M10","M16","M19","M20")){
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
  
  save(me_reg_effect_list,file="data_user/02_DEGenesIsoforms/02_02_A_02_ParikshakME_RegEffect_Permutations.RData")
  
}

##### (3) Test Severity #####

if(SIGTEST==TRUE){
  
  load("data_provided/02_DEGenesIsoforms/02_01_A_AllProcessedData_wModelMatrix.RData")
  load("data_provided/02_DEGenesIsoforms/02_01_A_RegressedExpression.RData")
  load("data_user/02_DEGenesIsoforms/02_02_A_02_ParikshakME_RegEffect_Permutations.RData")
  
  ### Parikshak modules
  
  parik_mods = read.csv("data_provided/02_DEGenesIsoforms/02_01_A_ParikshakModules.csv")
  
  # extract sig_mods
  
  sig_mods = c(4,9,10,16,19,20)
  keep = which(parik_mods$WGCNA.Module.Label %in% sig_mods)
  parik_mods = parik_mods[keep,]
  
  int_genes = intersect(parik_mods$ENSEMBL.ID,substr(rownames(datExpr.reg),1,15))
  parik_mods = parik_mods[match(int_genes,parik_mods$ENSEMBL.ID),]
  datExpr.reg_match = datExpr.reg[match(int_genes,substr(rownames(datExpr.reg),1,15)),]
  
  parik_mes = moduleEigengenes(t(datExpr.reg_match), colors=parik_mods$WGCNA.Module.Label, softPower=8)
  parik_mes = parik_mes$eigengenes
  colnames(parik_mes) = gsub("E","",colnames(parik_mes))
  
  # Make tables for association of MEs with region-specific ASD
  
  assoc_table_p <- assoc_table_beta <-  data.frame(matrix(NA,ncol(parik_mes),11))
  colnames(assoc_table_p) <- colnames(assoc_table_beta) <-   c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                                                               "BA39.40","BA17","BA41.42.22","BA20.37","BA38")
  rownames(assoc_table_p) <- rownames(assoc_table_beta) <-   colnames(parik_mes)
  
  for(j in c(1:ncol(assoc_table_p))){
    reg = colnames(assoc_table_p)[j]
    for(me in colnames(parik_mes)){
      assoc_table_p[me,reg] = signif(me_reg_effect_list[[me]][[reg]]$p[1],3)
      assoc_table_beta[me,reg] = signif(me_reg_effect_list[[me]][[reg]]$beta[1],3)
    }
  }
  
  ### now get permuation p-values and plots for the individual regions
  
  ME_names=factor(rownames(assoc_table_p),levels = c("M4","M10","M16","M9","M19","M20"))
  
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
  
  pdf(file="plots/02_DEGenesIsoforms/02_02_A_03_ParikshakME_Permutation.pdf")
  
  for(i in c(1:dim(ME_assoc_table_test)[1])){
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
  
  ME_assoc_table_test$Perm_sig = rep(NA,dim(ME_assoc_table_test)[1])
  ME_assoc_table_test$Perm_sig[which(round(ME_assoc_table_test$Perm_P,2) <= 0.05)] = "*"
  
  ME_assoc_table_test$FDR_Perm_P = p.adjust(ME_assoc_table_test$Perm_P,method = "fdr")
  
  ME_assoc_table_ParikshakME_test = ME_assoc_table_test
  
  save(ME_assoc_table_ParikshakME_test,file="data_user/02_DEGenesIsoforms/02_02_A_03_ParikshakME_Permutation_Results.RData")
  
}
