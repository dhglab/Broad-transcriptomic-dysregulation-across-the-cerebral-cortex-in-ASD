##### 02_03_A_RegSpecificASD_Permutation.R
##### ASD region-specific DGE permutation (Gene-Level)
##### Determine if region-specific ASD effect is significantly different from the whole cortex ASD effect.
##### May 2020, Jillian Haney

options(stringsAsFactors = FALSE)

library(WGCNA); library(limma); library(statmod); library(ggplot2); library(reshape2)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

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
  
  ### Run 10,001 permutations (permutation=i) as below:
  ### Rscript 02_03_A_RegSpecificASD_Permutation.R i
  
  args = commandArgs(trailingOnly=TRUE)
  iter=as.numeric(args[1])
  
  datMeta_model$Diagnosis = factor(datMeta$Diagnosis,levels=c("CTL","ASD","Dup15q"))
  datMeta_model$Region = as.factor(datMeta$region)
  
  ### permute regional identity
  
  if(iter==1){
    datMeta_model_perm = datMeta_model
    dir.create("data_user/02_DEGenesIsoforms/02_03_A_01_RegSpecificASD_permutations")
  }else{
    set.seed(as.numeric(iter-1))
    datMeta_model_perm = datMeta_model
    for(sub in unique(datMeta_model$Subject)){
      idx = which(datMeta_model$Subject==sub)
      datMeta_model_perm$Region[idx]=sample(datMeta_model_perm$Region[idx])
    }
  }
  
  rm(datMeta_model)
  
  datMeta_model_perm$Region=gsub("-",".",datMeta_model_perm$Region)
  datMeta_model_perm$DxRegion = factor(paste0(datMeta_model_perm$Diagnosis, ".", datMeta_model_perm$Region))
  
  mod=model.matrix(~ 0 + DxRegion+SeqBatch+Sex+Ancestry,data=datMeta_model_perm)
  design=data.frame(mod,datMeta_model_perm[,c(6:21)])
  
  corfit <- duplicateCorrelation(datExpr,design,block=datMeta_model_perm$Subject)
  fit <- lmFit(datExpr, design,block=datMeta_model_perm$Subject,correlation=corfit$consensus)
  
  ### just get region specific effects of ASD
  
  assoc_table_p <- assoc_table_beta <- data.frame(matrix(NA,nrow(datExpr),11))
  colnames(assoc_table_p) <- colnames(assoc_table_beta) <- c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                                                             "BA39.40","BA17","BA41.42.22","BA20.37","BA38")
  rownames(assoc_table_p) <- rownames(assoc_table_beta) <- rownames(datExpr)
  
  ASD_Region_ttables=list()
  
  for(j in c(1:ncol(assoc_table_p))){
    print(j)
    reg = colnames(assoc_table_p)[j]
    form = paste("DxRegionASD.",reg," - DxRegionCTL.",reg,sep="")
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
  
  assoc_table_p_asdReg = assoc_table_p
  assoc_table_beta_asdReg = assoc_table_beta
  
  save(ASD_Region_ttables,assoc_table_p,assoc_table_beta,file=paste("data_user/02_DEGenesIsoforms/02_03_A_01_RegSpecificASD_permutations/ASD_DGE_permutation_",iter,".RData",sep=""))
  
}


##### (2) Compile Permutation Results #####

if(COMPILE==TRUE){
  
  ## iter=1 is the real result
  
  dge_reg_effect_list = list()
  
  for(i in c(1:10001)){
    if (i%%100 == 0) {print(paste(i,"/",length(list.files("data_user/02_DEGenesIsoforms/02_03_A_01_RegSpecificASD_permutations/")),sep=""))}
    load(paste("data_user/02_DEGenesIsoforms/02_03_A_01_RegSpecificASD_permutations/ASD_DGE_permutation_",i,".RData",sep=""))
    for(reg in c("BA9","BA24","BA44.45","BA4.6","BA3.1.2.5","BA7",
                 "BA39.40","BA17","BA41.42.22","BA20.37","BA38")){
      if(i==1){
        dge_reg_effect_list[[reg]]=length(which(assoc_table_p[,reg] < 0.05))
      }else{
        dge_reg_effect_list[[reg]]=c(dge_reg_effect_list[[reg]],length(which(assoc_table_p[,reg] < 0.05)))
      }
    }
  }
  
  save(dge_reg_effect_list,file="data_user/02_DEGenesIsoforms/02_03_A_02_ASDDGE_RegSpecific_Permutations.RData")
  
}

##### (3) Test Severity #####

if(SIGTEST==TRUE){
  
  load("data_provided/02_DEGenesIsoforms/02_01_A_AllProcessedData_wModelMatrix.RData")
  load("data_provided/02_DEGenesIsoforms/02_01_A_RegressedExpression.RData")
  load("data_user/02_DEGenesIsoforms/02_03_A_02_ASDDGE_RegSpecific_Permutations.RData")
  
  ### now get permuation p-values and plots for the individual regions
  
  reg_names=names(dge_reg_effect_list)
  
  dge_assoc_table = data.frame("Region"=reg_names,
                               "Real_NumDGE"=rep(NA,length(reg_names)),
                               "Permutation_P"=rep(NA,length(reg_names)),
                               "Permutation_Sig"=rep(NA,length(reg_names)))
  
  pdf(file="plots/02_DEGenesIsoforms/02_03_A_02_RegSpecificASD_DEGenes_Permutation.pdf")
  
  for(i in c(1:11)){
    
    reg=dge_assoc_table$Region[i]
    real_num = dge_reg_effect_list[[reg]][1]
    tmp=dge_reg_effect_list[[reg]][-1]
    
    cent=scale(tmp,scale=FALSE)
    mn=mean(tmp)
    diff=real_num-mn
    more=length(which(abs(cent) >= abs(diff)))
    
    p = (more)/10001
    dge_assoc_table$Permutation_P[i] = signif(p,5)
    dge_assoc_table$Real_NumDGE[i] = real_num
    
    if(p < 0.05){
      dge_assoc_table$Permutation_Sig[i] = "*"
    }
    
    plot.perm=data.frame(tmp)
    
    print(ggplot(plot.perm,aes(x=tmp))+
            geom_density()+
            ggtitle(paste(reg,"; p=",signif(p,3),sep=""))+
            xlab("Signed log10 Permutation P-Value")+
            theme(plot.title = element_text(hjust = 0.5))+
            geom_vline(xintercept=real_num,col="red"))
  }
  
  dev.off()
  
  save(dge_assoc_table,file="data_user/02_DEGenesIsoforms/02_03_A_02_RegSpecificASD_DEGenes.RData")

}
