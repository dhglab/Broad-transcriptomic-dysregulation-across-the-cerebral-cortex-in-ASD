##### 03_01_A_TRI.R
##### Identify any changes in Transcriptomic Regional Identity (TRI) in ASD (Gene-Level)
##### November 2020, Jillian Haney
##### Note: many variables are named 'cp_xyz'; cp stands for 'cortical patterning', another definition for TRI

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(ggplot2); library(plyr); library(gridExtra); library(limma)
library(ggpubr); library(statmod); library(biomaRt)

##### Pipeline #####

### (1) TRI Permutation Analysis
### (2) TRI Bootstrap Analysis
### (3) Identify ARI Genes

##### (1) TRI Permutation Analysis #####

load("data_provided/03_TRI/03_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/03_TRI/03_01_A_RegressedExpression_TRI.RData") ### produced by 01_02_A_CountsProcessing.R, section 8

### First: establish subjects/samples that I will use for each of the 55 comparisons
### (paired, so that each subject has a sample for each region in the comparison)

datMeta_model = data.frame(datMeta_model,"Region"=as.factor(datMeta$region),"Diagnosis"=as.factor(datMeta$Diagnosis))

cont_mat=t(combn(paste(unique(datMeta_model$Region),sep=""), 2))
colnames(cont_mat)=c("Region1","Region2")
idx_comp_all=list()

for(i in c(1:nrow(cont_mat))){
  idx=which(datMeta_model$Region==cont_mat[i,1] | datMeta_model$Region==cont_mat[i,2])
  datMeta_comp=datMeta_model[idx,]
  for(sub in unique(datMeta_comp$Subject)){
    length = length(unique(datMeta_comp$Region[which(datMeta_comp$Subject==sub)]))
    if(length < 2 ){
      datMeta_comp=datMeta_comp[-which(datMeta_comp$Subject==sub),]
    }
  }
  idx_comp=list()
  for(dx in c("ASD","CTL")){
    idx2=which(datMeta_comp$Diagnosis==dx)
    datMeta_comp_dx=datMeta_comp[idx2,]
    idx_keep=which(duplicated(datMeta_comp_dx$Subject)==TRUE)
    subs_keep=datMeta_comp_dx$Subject[idx_keep]
    idx_comp[[dx]]=which(datMeta_model$Subject %in% subs_keep & (datMeta_model$Region==cont_mat[i,1] | datMeta_model$Region==cont_mat[i,2]) )
    sub_xtra_remove=names(table(datMeta_model$Subject[idx_comp[[dx]]]))[which(table(datMeta_model$Subject[idx_comp[[dx]]]) > 2)]
    if(length(sub_xtra_remove)==0){
      next
    }
    tmp0=datMeta_model[idx_comp[[dx]],]
    idx_remove=which(tmp0$Subject %in% sub_xtra_remove)
    tmp=tmp0[idx_remove,]
    tmp=tmp[order(tmp$Subject),]
    
    dup_samps=NA
    for(samp in unique(tmp$Subject)){
      dup=which(duplicated(tmp$Region[which(tmp$Subject==samp)])==TRUE)
      dup_samps=c(dup_samps,paste(tmp$Subject[which(tmp$Subject==samp)][dup],
                                  tmp$Region[which(tmp$Subject==samp)][dup],
                                  tmp$picard_rnaseq.PCT_INTERGENIC_BASES[which(tmp$Subject==samp)][dup],sep="_"))
    }
    dup_samps=dup_samps[-1]
    
    idx_remove2=NA
    for(dup in dup_samps){
      dup_string=strsplit(dup,split="_")
      idx_remove2=c(idx_remove2,
                    which(tmp0$Subject==dup_string[[1]][1] & 
                            tmp0$Region==dup_string[[1]][2] &
                            tmp0$picard_rnaseq.PCT_INTERGENIC_BASES==dup_string[[1]][3]))
    }
    idx_remove2=idx_remove2[-1]
    
    idx_comp[[dx]]=idx_comp[[dx]][-idx_remove2]
  }
  idx_comp_all[[paste(cont_mat[i,1],cont_mat[i,2],sep=" - ")]]=c(idx_comp[["CTL"]],idx_comp[["ASD"]])
}

n_vec_asd <- n_vec_ctl <- rep(NA,55)
names(n_vec_asd) <- names(n_vec_ctl) <- names(idx_comp_all)

for(comp in names(idx_comp_all)){
  datMeta_comp=datMeta_model[idx_comp_all[[comp]],]
  datMeta_unq=datMeta_comp[!duplicated(datMeta_comp$Subject),]
  asd_subs=length(datMeta_unq$Subject[which(datMeta_unq$Diagnosis=="ASD")])
  ctl_subs=length(datMeta_unq$Subject[which(datMeta_unq$Diagnosis=="CTL")])
  n_vec_asd[comp]=asd_subs
  n_vec_ctl[comp]=ctl_subs
}

save(idx_comp_all,n_vec_asd,n_vec_ctl,file="data_user/03_TRI/03_01_A_01_Permutation_Index.RData")
save(datMeta_model,file="data_user/03_TRI/03_01_A_01_datMetaModel_TRI.RData")

### Run Wilcoxon test with extracted sample pairs (get 'true' values of #DGE between regions)
### Get numbers of DE genes between pairs of regions (stratified by ASD and CTL)

cp_list=list()

for(comp in names(idx_comp_all)){
  print(comp)
  for(dx in c("ASD","CTL")){
    print(dx)
    idx_dx=which(datMeta_model$Diagnosis[idx_comp_all[[comp]]]==dx)
    
    comp_meta_mat=datMeta_model[idx_comp_all[[comp]][idx_dx],]
    comp_expr_mat=datExpr.reg[,idx_comp_all[[comp]][idx_dx]]
    names=unique(comp_meta_mat$Subject)
    
    reg1_name=strsplit(comp,split=" - ")[[1]][1]
    reg2_name=strsplit(comp,split=" - ")[[1]][2]
    
    reg1=comp_expr_mat[,which(comp_meta_mat$Region==reg1_name)]
    reg2=comp_expr_mat[,which(comp_meta_mat$Region==reg2_name)]
    colnames(reg1)=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg1_name)]
    colnames(reg2)=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg2_name)]
    
    reg1=reg1[,match(names,colnames(reg1))]
    reg2=reg2[,match(names,colnames(reg2))]
    
    ttable=data.frame("V"=rep(NA,nrow(comp_expr_mat)),"PVal"=rep(NA,nrow(comp_expr_mat)))
    rownames(ttable)=rownames(comp_expr_mat)
    
    for(i in c(1:nrow(comp_expr_mat))){
      if (i%%1000 == 0) {print(i)}
      tmp=wilcox.test(reg1[i,],reg2[i,],paired=TRUE)
      ttable[i,1]=tmp$statistic
      ttable[i,2]=tmp$p.value
    }
    
    cp_list[[comp]][[dx]]$data=ttable
    cp_list[[comp]][[dx]]$n=ncol(reg1)
  }
}

save(cp_list,file=paste("data_user/03_TRI/03_01_A_01_True_TRI_DGE_CompleteList.RData",sep=""))

### Compile results into a matrix

cp_mat <- data.frame(matrix(NA,nrow=length(names(cp_list)),ncol=5))
rownames(cp_mat) <- names(cp_list)
colnames(cp_mat) <- c("Comparison","Num_DGE_ASD","Num_DGE_CTL","Num_Samps_ASD","Num_Samps_CTL")
cp_mat$Comparison <- names(cp_list)
for(comp in cp_mat$Comparison){
  tmp=cp_list[[which(names(cp_list)==comp)]]$ASD$data
  tmp$FDR=p.adjust(tmp$PVal,method="fdr")
  cp_mat$Num_DGE_ASD[which(cp_mat$Comparison==comp)]=length(which(tmp$FDR < 0.05))
  tmp=cp_list[[which(names(cp_list)==comp)]]$CTL$data
  tmp$FDR=p.adjust(tmp$PVal,method="fdr")
  cp_mat$Num_DGE_CTL[which(cp_mat$Comparison==comp)]=length(which(tmp$FDR < 0.05))
  cp_mat$Num_Samps_ASD[which(cp_mat$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$ASD$n
  cp_mat$Num_Samps_CTL[which(cp_mat$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$CTL$n
}

save(cp_mat,file="data_user/03_TRI/03_01_A_01_True_TRI_DGE_SummarizedMatrix.RData")

### Get the difference in the # of DE genes between regions for ASD and CTL

diff_list=list()
### take difference of num dge CTL - num DGE ASD. Positive diff = CTL greater, Neg diff = ASD greater.

for(comp in names(idx_comp_all)){
  diff_list[[comp]]$Difference=cp_mat[comp,3]-cp_mat[comp,2]
  diff_list[[comp]]$ASD_N=cp_mat[comp,4]
  diff_list[[comp]]$CTL_N=cp_mat[comp,5]
}

save(diff_list,file="data_user/03_TRI/03_01_A_01_True_TRI_DGE_ASD_v_Control.RData")

### then do this again, but permute the ASD/Control status of the subjects (get difference in DGE each time)

rm(list=ls())

load("data_provided/03_TRI/03_01_A_AllProcessedData_wModelMatrix.RData")
load("data_provided/03_TRI/03_01_A_RegressedExpression_TRI.RData")
rm(datMeta_model)
load("data_user/03_TRI/03_01_A_01_datMetaModel_TRI.RData")
load("data_user/03_TRI/03_01_A_01_Permutation_Index.RData")

### 10,000 permutations were run for the original analysis.
### It is recommended to run individual permutations in parallel on a high-performance computing cluster.

for(iter in c(1:10000)){
  
  set.seed(iter)
  
  datMeta_model_perm=list()
  
  for(comp in names(idx_comp_all)){
    datMeta_tmp=datMeta_model[idx_comp_all[[comp]],]
    datMeta_tmp_unq=datMeta_tmp[!duplicated(datMeta_tmp$Subject),]
    dx_permute=c(rep("CTL",length(which(datMeta_tmp_unq$Diagnosis=="CTL"))),
                 rep("ASD",length(which(datMeta_tmp_unq$Diagnosis=="ASD"))))
    dx_permute=sample(dx_permute)
    datMeta_tmp_unq$Diagnosis=dx_permute
    datMeta_tmp$Diagnosis=datMeta_tmp_unq$Diagnosis[match(datMeta_tmp$Subject,datMeta_tmp_unq$Subject)]
    datMeta_model_perm[[comp]]=datMeta_tmp
  }
  
  datExpr.reg_perm=datExpr.reg
  colnames(datExpr.reg_perm)=paste(datMeta_model$Subject,datMeta_model$Region,sep="_")
  
  cp_list=list()
  
  for(comp in names(idx_comp_all)){
    print(comp)
    for(dx in c("ASD","CTL")){
      print(dx)
      idx_dx=which(datMeta_model_perm[[comp]]$Diagnosis==dx)
      comp_meta_mat=datMeta_model_perm[[comp]][idx_dx,]
      names_expr=paste(comp_meta_mat$Subject,comp_meta_mat$Region,sep="_")
      
      comp_expr_mat=datExpr.reg_perm[,match(names_expr,colnames(datExpr.reg_perm))]
      
      names=unique(comp_meta_mat$Subject)
      
      reg1_name=strsplit(comp,split=" - ")[[1]][1]
      reg2_name=strsplit(comp,split=" - ")[[1]][2]
      
      reg1=comp_expr_mat[,which(comp_meta_mat$Region==reg1_name)]
      reg2=comp_expr_mat[,which(comp_meta_mat$Region==reg2_name)]
      colnames(reg1)=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg1_name)]
      colnames(reg2)=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg2_name)]
      
      reg1=reg1[,match(names,colnames(reg1))]
      reg2=reg2[,match(names,colnames(reg2))]
      
      ttable=data.frame("V"=rep(NA,nrow(comp_expr_mat)),"PVal"=rep(NA,nrow(comp_expr_mat)))
      rownames(ttable)=rownames(comp_expr_mat)
      
      for(i in c(1:nrow(comp_expr_mat))){
        if (i%%1000 == 0) {print(i)}
        tmp=wilcox.test(reg1[i,],reg2[i,],paired=TRUE)
        ttable[i,1]=tmp$statistic
        ttable[i,2]=tmp$p.value
      }
      
      cp_list[[comp]][[dx]]$data=ttable
      cp_list[[comp]][[dx]]$n=ncol(reg1)
    }
  }
  
  save(cp_list,file=paste("data_user/03_01_A_01_permutation/TRI_list_permutation",iter,".RData",sep=""))

}

### Compile permutations

cp_mat_list=list()

for(i in c(1:length(list.files("data_user/03_01_A_01_permutation/")))){
  if (i%%50 == 0) {print(paste(i,"/",length(list.files("data_user/03_01_A_01_permutation/")),sep=""))}
  load(paste("data_user/03_01_A_01_permutation/TRI_list_permutation",i,".RData",sep=""))
  cp_mat_list[[i]] <- data.frame(matrix(NA,nrow=length(names(cp_list)),ncol=5))
  rownames(cp_mat_list[[i]]) <- names(cp_list)
  colnames(cp_mat_list[[i]]) <- c("Comparison","Num_DGE_ASD","Num_DGE_CTL","Num_Samps_ASD","Num_Samps_CTL")
  cp_mat_list[[i]]$Comparison <- names(cp_list)
  for(comp in cp_mat_list[[i]]$Comparison){
    tmp=cp_list[[which(names(cp_list)==comp)]]$ASD$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    cp_mat_list[[i]]$Num_DGE_ASD[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
    tmp=cp_list[[which(names(cp_list)==comp)]]$CTL$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    cp_mat_list[[i]]$Num_DGE_CTL[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
    cp_mat_list[[i]]$Num_Samps_ASD[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$ASD$n
    cp_mat_list[[i]]$Num_Samps_CTL[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$CTL$n
  }
}

save(cp_mat_list,file="data_user/03_TRI/03_01_A_01_Permuted_TRI_DGE_SummarizedList.RData")

### Now, use permutations to create a null distribution to test if 
### the # of DE genes between ASD and Controls is significant for each regional pair.

rm(list=ls())

load("data_user/03_TRI/03_01_A_01_Permuted_TRI_DGE_SummarizedList.RData")
load("data_user/03_TRI/03_01_A_01_Permutation_Index.RData")

perm_diff_list=list()

for(comp in names(idx_comp_all)){
  perm_diff_list[[comp]]$data=data.frame(matrix(NA,nrow=10000,ncol=4))
  colnames(perm_diff_list[[comp]]$data)=c("Iteration","Difference","ASD_Num_DGE","CTL_Num_DGE")
  rownames(perm_diff_list[[comp]]$data)=seq(1:10000)
}

for(i in c(1:10000)){
  if (i%%50 == 0) {print(paste(i,"/10000",sep=""))}
  cp_mat_tmp=cp_mat_list[[i]]
  for(comp in names(idx_comp_all)){
    perm_diff_list[[comp]]$data$Iteration[i]=i
    perm_diff_list[[comp]]$data$Difference[i]=cp_mat_tmp[comp,3]-cp_mat_tmp[comp,2]
    perm_diff_list[[comp]]$data$ASD_Num_DGE[i]=cp_mat_tmp[comp,2]
    perm_diff_list[[comp]]$data$CTL_Num_DGE[i]=cp_mat_tmp[comp,3]
  }
}

save(perm_diff_list,file="data_user/03_TRI/03_01_A_01_Permuted_TRI_DGE_ASD_v_Control.RData")

load("data_user/03_TRI/03_01_A_01_True_TRI_DGE_ASD_v_Control.RData")

## to get p-value: # of observations greater than true value / 10001
## use absolute value of the difference for significance testing

ttable=data.frame(matrix(NA,nrow=length(cp_mat_list[[1]]$Comparison),ncol=6))
colnames(ttable) = c("Comparison","CTL_ASD_NumDGE_MeanDiff","Perm_PVal","FDR_adj_Perm_PVal","ASD_N","CTL_N")
rownames(ttable) <- ttable$Comparison <- cp_mat_list[[1]]$Comparison

pdf(file="plots/03_TRI/03_01_A_01_TRI_Permutation.pdf")

for( comp in names(idx_comp_all)){
  cent=round(scale(perm_diff_list[[comp]]$data$Difference,scale=FALSE),0)
  mn=round(mean(perm_diff_list[[comp]]$data$Difference),0)
  diff=diff_list[[comp]]$Difference-mn
  ttable$CTL_ASD_NumDGE_MeanDiff[which(rownames(ttable)==comp)]=diff
  
  more=length(which(abs(cent) >= abs(diff)))
  p = (more)/10001
  ttable$Perm_PVal[which(rownames(ttable)==comp)]=signif(p,5)
  ttable$ASD_N[which(rownames(ttable)==comp)]=diff_list[[comp]]$ASD_N
  ttable$CTL_N[which(rownames(ttable)==comp)]=diff_list[[comp]]$CTL_N
  
  plot.cent=data.frame(cent)
  
  print(ggplot(plot.cent,aes(x=cent))+
          geom_density()+
          ggtitle(paste(comp,"; Real_Diff=",diff,"; p=",signif(p,3),"; N_CTL=",diff_list[[comp]]$CTL_N,
                        "; N_ASD=",diff_list[[comp]]$ASD_N,sep=""))+
          xlab("CTL_NDiff-ASD_NDiff (Mean Centered)")+
          theme(plot.title = element_text(hjust = 0.5))+
          geom_vline(xintercept=diff,col="red"))
}

dev.off()

ttable$FDR_adj_Perm_PVal=p.adjust(ttable$Perm_PVal,method="fdr")

write.csv(ttable,file="data_user/03_TRI/03_01_A_01_Permutation_RegComp-TTable.csv")

##### (2) TRI Bootstrap Analysis #####

### First, subsample down to 10 matched ASD cases and controls
  
set.seed(101788)

load("data_provided/03_TRI/03_01_A_AllProcessedData_wModelMatrix.RData")
load("data_user/03_TRI/03_01_A_01_Permutation_Index.RData")
  
idx_cp_bs <- list()

for(comp in names(idx_comp_all)){
  print(comp)
  idx = idx_comp_all[[comp]]
  datMeta_tmp = datMeta[idx,]
  datMeta_tmp_ctl = datMeta_tmp[which(datMeta_tmp$Diagnosis=="CTL"),]
  datMeta_tmp_asd = datMeta_tmp[which(datMeta_tmp$Diagnosis=="ASD"),]
  ctl_subs = unique(datMeta_tmp_ctl$subject)
  asd_subs = unique(datMeta_tmp_asd$subject)
  ctl_subs_age = datMeta_tmp_ctl$Age[match(ctl_subs,datMeta_tmp_ctl$subject)]
  asd_subs_age = datMeta_tmp_asd$Age[match(asd_subs,datMeta_tmp_asd$subject)]
  names(ctl_subs_age) <- ctl_subs
  names(asd_subs_age) <- asd_subs
  
  ### match for age
  
  med_age = median(datMeta_tmp$Age)
  
  # Controls
  
  print("Control")  
  mod=0
  idx_keep=0
  while(length(idx_keep) < 10 ){
    idx_keep = which(ctl_subs_age >= (med_age-mod) & ctl_subs_age <= (med_age+mod))
    mod=mod+1
    if(mod%%10==0){print(mod)}
  }
  ctl_subs = ctl_subs[idx_keep]
  ctl_subs_age = ctl_subs_age[idx_keep]
  if(length(ctl_subs) > 10){
    while(length(ctl_subs) > 10){
      ctl_subs_age = ctl_subs_age[order(ctl_subs_age)]
      extremes = c(names(ctl_subs_age)[1],names(ctl_subs_age)[length(ctl_subs_age)])
      remove_samp = sample(extremes,1)
      ctl_subs = ctl_subs[-which(ctl_subs==remove_samp)]
      ctl_subs_age = ctl_subs_age[-which(names(ctl_subs_age)==remove_samp)]
    }
  }
  
  # ASD
  
  print("ASD")  
  mod=0
  idx_keep=0
  while(length(idx_keep) < 10 ){
    idx_keep = which(asd_subs_age >= (med_age-mod) & asd_subs_age <= (med_age+mod))
    mod=mod+1
    if(mod%%10==0){print(mod)}
  }
  asd_subs = asd_subs[idx_keep]
  asd_subs_age = asd_subs_age[idx_keep]
  if(length(asd_subs) > 10){
    while(length(asd_subs) > 10){
      asd_subs_age = asd_subs_age[order(asd_subs_age)]
      extremes = c(names(asd_subs_age)[1],names(asd_subs_age)[length(asd_subs_age)])
      remove_samp = sample(extremes,1)
      asd_subs = asd_subs[-which(asd_subs==remove_samp)]
      asd_subs_age = asd_subs_age[-which(names(asd_subs_age)==remove_samp)]
    }
  }
  
  # Add both new indeces to list
  
  subjects_keep = c(ctl_subs,asd_subs)   
  sample_names = rownames(datMeta_tmp)[which(datMeta_tmp$subject %in% subjects_keep)]
  idx_cp_bs[[comp]] <- match(sample_names,rownames(datMeta))
}

save(idx_cp_bs,file="data_user/03_TRI/03_01_A_02_TRI_Bootstrap_Index_N10.RData")

### Run 10,000 bootstraps
### It is recommended to run individual bootstraps in parallel on a high-performance computing cluster.

for(iter in c(1:10000)){
  
  set.seed(iter)
  
  load("data_provided/03_TRI/03_01_A_AllProcessedData_wModelMatrix.RData")
  load("data_provided/03_TRI/03_01_A_RegressedExpression_TRI.RData")
  load("data_user/03_TRI/03_01_A_02_TRI_Bootstrap_Index_N10.RData")
  
  datMeta_model = data.frame(datMeta_model,"Region"=as.factor(datMeta$region),"Diagnosis"=as.factor(datMeta$Diagnosis))
  datMeta_model$Sample_ID = rownames(datMeta)
  
  datMeta_model_boot=list()
  
  for(comp in names(idx_cp_bs)){
    datMeta_tmp=datMeta_model[idx_cp_bs[[comp]],]
    datMeta_tmp_unq=datMeta_tmp[!duplicated(datMeta_tmp$Subject),]
    ctl_subs = datMeta_tmp_unq$Subject[which(datMeta_tmp_unq$Diagnosis=="CTL")]
    asd_subs = datMeta_tmp_unq$Subject[which(datMeta_tmp_unq$Diagnosis=="ASD")]
    
    samples_boot = c(as.character(sample(ctl_subs,10,replace=TRUE)),
                     as.character(sample(asd_subs,10,replace=TRUE)))
    
    i=1
    for(samp in samples_boot){
      if(i==1){
        datMeta_boot_add = datMeta_tmp[which(datMeta_tmp$Subject==samp),]
      }else{
        datMeta_boot_add = data.frame(rbind(datMeta_boot_add,datMeta_tmp[which(datMeta_tmp$Subject==samp),]))
      }
      i=i+1  
    }
    
    datMeta_model_boot[[comp]]=datMeta_boot_add
  }
  
  cp_list=list()
  
  for(comp in names(idx_cp_bs)){
    print(comp)
    for(dx in c("ASD","CTL")){
      print(dx)
      
      comp_meta_mat = datMeta_model_boot[[comp]]
      comp_expr_mat = datExpr.reg[,match(comp_meta_mat$Sample_ID,colnames(datExpr.reg))]
      
      idx_dx = which(comp_meta_mat$Diagnosis==dx)
      comp_meta_mat = comp_meta_mat[idx_dx,]
      comp_expr_mat = comp_expr_mat[,idx_dx]
      
      reg1_name=strsplit(comp,split=" - ")[[1]][1]
      reg2_name=strsplit(comp,split=" - ")[[1]][2]
      
      names=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg1_name)]
      
      reg1=comp_expr_mat[,which(comp_meta_mat$Region==reg1_name)]
      reg2=comp_expr_mat[,which(comp_meta_mat$Region==reg2_name)]
      colnames(reg1)=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg1_name)]
      colnames(reg2)=comp_meta_mat$Subject[which(comp_meta_mat$Region==reg2_name)]
      
      reg1=reg1[,match(names,colnames(reg1))]
      reg2=reg2[,match(names,colnames(reg2))]
      
      ttable=data.frame("V"=rep(NA,nrow(comp_expr_mat)),"PVal"=rep(NA,nrow(comp_expr_mat)))
      rownames(ttable)=rownames(comp_expr_mat)
      
      for(i in c(1:nrow(comp_expr_mat))){
        if (i%%1000 == 0) {print(i)}
        tmp=wilcox.test(reg1[i,],reg2[i,],paired=TRUE)
        ttable[i,1]=tmp$statistic
        ttable[i,2]=tmp$p.value
      }
      
      cp_list[[comp]][[dx]]$data=ttable
      cp_list[[comp]][[dx]]$n=ncol(reg1)
    }
  }
  
  save(cp_list,file=paste("data_user/03_TRI/03_01_A_02_bootstrap/TRI_list_bootstrap",iter,".RData",sep=""))
}  

##### Compile bootstrap output

cp_mat_list=list()

for(i in c(1:length(list.files("data_user/03_TRI/03_01_A_02_bootstrap/")))){
  if (i%%50 == 0) {print(paste(i,"/",length(list.files("data_user/03_TRI/03_01_A_02_bootstrap/")),sep=""))}
  load(paste("data_user/03_TRI/03_01_A_02_bootstrap/TRI_list_bootstrap",i,".RData",sep=""))
  cp_mat_list[[i]] <- data.frame(matrix(NA,nrow=length(names(cp_list)),ncol=7))
  rownames(cp_mat_list[[i]]) <- names(cp_list)
  colnames(cp_mat_list[[i]]) <- c("Comparison","Num_DGE_ASD","Num_DGE_CTL","Num_DGE_ASD_P","Num_DGE_CTL_P","Num_Samps_ASD","Num_Samps_CTL")
  cp_mat_list[[i]]$Comparison <- names(cp_list)
  if(i==1){
    pdf(file="Bootstrap1_PValue_Histograms.pdf")
    par(mfrow=c(1,2))
    for(comp in cp_mat_list[[i]]$Comparison){
      tmp=cp_list[[which(names(cp_list)==comp)]]$ASD$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      cp_mat_list[[i]]$Num_DGE_ASD[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
      cp_mat_list[[i]]$Num_DGE_ASD_P[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$PVal < 0.05))
      hist(tmp$PVal,main=paste("ASD",comp,sep=": "))
      tmp=cp_list[[which(names(cp_list)==comp)]]$CTL$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      cp_mat_list[[i]]$Num_DGE_CTL[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
      cp_mat_list[[i]]$Num_DGE_CTL_P[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$PVal < 0.05))
      hist(tmp$PVal,main=paste("CTL",comp,sep=": "))
      cp_mat_list[[i]]$Num_Samps_ASD[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$ASD$n
      cp_mat_list[[i]]$Num_Samps_CTL[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$CTL$n
    }
    dev.off()
  }else{
    for(comp in cp_mat_list[[i]]$Comparison){
      tmp=cp_list[[which(names(cp_list)==comp)]]$ASD$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      cp_mat_list[[i]]$Num_DGE_ASD[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
      cp_mat_list[[i]]$Num_DGE_ASD_P[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$PVal < 0.05))
      tmp=cp_list[[which(names(cp_list)==comp)]]$CTL$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      cp_mat_list[[i]]$Num_DGE_CTL[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$FDR < 0.05))
      cp_mat_list[[i]]$Num_DGE_CTL_P[which(cp_mat_list[[i]]$Comparison==comp)]=length(which(tmp$PVal < 0.05))
      cp_mat_list[[i]]$Num_Samps_ASD[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$ASD$n
      cp_mat_list[[i]]$Num_Samps_CTL[which(cp_mat_list[[i]]$Comparison==comp)]=cp_list[[which(names(cp_list)==comp)]]$CTL$n
    }
  }
}

save(cp_mat_list,file="data_user/03_TRI/03_01_A_02_TRI_Bootstrap_FullResultsList.RData")
  
### Get mean of each bootstrap for each comparison, for CTL and ASD separately

load("data_user/03_TRI/03_01_A_02_TRI_Bootstrap_FullResultsList.RData")

cp_list = cp_mat_list

cp_mean_table <- data.frame(matrix(NA,nrow=length(cp_list[[1]]$Comparison),ncol=7))
colnames(cp_mean_table) <- c("Comparison","NumDiff_ASD_P","NumDiff_CTL_P","%_ASD/CTL_P",
                             "NumDiff_ASD_FDR","NumDiff_CTL_FDR","%_ASD/CTL_FDR")
cp_mean_table$Comparison = cp_list[[1]]$Comparison

for(comp in cp_list[[1]]$Comparison){
  print(comp)
  mean_vec_ctl_p <- mean_vec_asd_p <- mean_vec_ctl_fdr <- mean_vec_asd_fdr <- rep(NA,10000)
  for(i in c(1:10000)){
    mean_vec_ctl_p[i] = cp_list[[i]]$Num_DGE_CTL_P[which(cp_list[[i]]$Comparison==comp)]
    mean_vec_asd_p[i] = cp_list[[i]]$Num_DGE_ASD_P[which(cp_list[[i]]$Comparison==comp)]
    mean_vec_ctl_fdr[i] = cp_list[[i]]$Num_DGE_CTL[which(cp_list[[i]]$Comparison==comp)]
    mean_vec_asd_fdr[i] = cp_list[[i]]$Num_DGE_ASD[which(cp_list[[i]]$Comparison==comp)]
  }
  cp_mean_table$NumDiff_ASD_P[which(cp_mean_table$Comparison==comp)]=mean(mean_vec_asd_p)
  cp_mean_table$NumDiff_CTL_P[which(cp_mean_table$Comparison==comp)]=mean(mean_vec_ctl_p)
  cp_mean_table[which(cp_mean_table$Comparison==comp),4]=(mean(mean_vec_asd_p)/mean(mean_vec_ctl_p))
  cp_mean_table$NumDiff_ASD_FDR[which(cp_mean_table$Comparison==comp)]=mean(mean_vec_asd_fdr)
  cp_mean_table$NumDiff_CTL_FDR[which(cp_mean_table$Comparison==comp)]=mean(mean_vec_ctl_fdr)
  cp_mean_table[which(cp_mean_table$Comparison==comp),7]=(mean(mean_vec_asd_fdr)/mean(mean_vec_ctl_fdr))
  
}

write.csv(cp_mean_table,file="data_user/03_TRI/03_01_A_02_ASD_v_CTL_TRI_Boostrap_SummarizedResults.RData")

### and make an extended table with more statistics

cp_list = cp_mat_list

cp_mean_table <- data.frame(matrix(NA,nrow=length(cp_list[[1]]$Comparison),ncol=8))
colnames(cp_mean_table) <- c("Comparison","Mean_NumDiff_ASD","SD_ASD","95%_ASD",
                             "Mean_NumDiff_CTL","SD_CTL","95%_CTL","%_ASD/CTL")
cp_mean_table$Comparison=cp_list[[1]]$Comparison

for(comp in cp_list[[1]]$Comparison){
  
  mean_vec_ctl <- mean_vec_asd <- rep(NA,10000)
  mean_vec_ctl_fdr <- mean_vec_asd_fdr <- rep(NA,10000)
  
  for(i in c(1:10000)){
    mean_vec_ctl[i] = cp_list[[i]]$Num_DGE_CTL_P[which(cp_list[[i]]$Comparison==comp)]
    mean_vec_asd[i] = cp_list[[i]]$Num_DGE_ASD_P[which(cp_list[[i]]$Comparison==comp)]
    mean_vec_ctl_fdr[i] = cp_list[[i]]$Num_DGE_CTL[which(cp_list[[i]]$Comparison==comp)]
    mean_vec_asd_fdr[i] = cp_list[[i]]$Num_DGE_ASD[which(cp_list[[i]]$Comparison==comp)]
  }
  
  datPlot=data.frame("Control_P"=mean_vec_ctl,"Control_FDR"=mean_vec_ctl_fdr,
                     "ASD_P"=mean_vec_asd,"ASD_FDR"=mean_vec_asd_fdr)
  colnames(datPlot)[1]=paste(comp,": ",colnames(datPlot)[1],sep="")
  for(group in colnames(datPlot)){
    plot(density(datPlot[,group]),main=group,xlim=c(0,20000))
  }
  
  mean_asd=mean(mean_vec_asd)
  mean_ctl=mean(mean_vec_ctl)
  cp_mean_table$Mean_NumDiff_ASD[which(cp_mean_table$Comparison==comp)]=mean(mean_vec_asd)
  cp_mean_table$Mean_NumDiff_CTL[which(cp_mean_table$Comparison==comp)]=mean(mean_vec_ctl)
  cp_mean_table$`%_ASD/CTL`[which(cp_mean_table$Comparison==comp)]=(mean(mean_vec_asd)/mean(mean_vec_ctl))
  
  sd_asd=signif(sd(mean_vec_asd),5)
  sd_ctl=signif(sd(mean_vec_ctl),5)
  n=10000
  
  error_asd=qnorm(0.975)*(sd_asd/sqrt(n))
  error_ctl=qnorm(0.975)*(sd_ctl/sqrt(n))
  
  cp_mean_table$SD_ASD[which(cp_mean_table$Comparison==comp)]=sd_asd
  cp_mean_table$SD_CTL[which(cp_mean_table$Comparison==comp)]=sd_ctl
  cp_mean_table$`95%_ASD`[which(cp_mean_table$Comparison==comp)]=paste("[",signif((mean_asd-error_asd),5),
                                                                       ",",signif((mean_asd+error_asd),5),"]",sep="")
  cp_mean_table$`95%_CTL`[which(cp_mean_table$Comparison==comp)]=paste("[",signif((mean_ctl-error_ctl),5),
                                                                       ",",signif((mean_ctl+error_ctl),5),"]",sep="")
}

dev.off()

write.csv(cp_mean_table,file="data_user/03_TRI/03_01_A_02_ASD_v_CTL_TRI_Boostrap_DetailedResults.RData")

##### (3) Identify ARI Genes #####

rm(list=ls())

### Obtain genes from the Wilcoxon test (with true samples) that are missing in 95% of permutations
### It is recommended to run the iterations below in parallel on a high-performance computing cluster.
  
### First, extract the permuted genes for each comparison
  
idx_list=list(c(1:1000),c(1001:2000),c(2001:3000),c(3001:4000),c(4001:5000),
             c(5001:6000),c(6001:7000),c(7001:8000),c(8001:9000),c(9001:10000))

for(iter_idx in c(1:10)){
  
  idx = idx_list[[as.numeric(iter_idx)]]
  
  perm_diff_genes_ctl=list()
  perm_diff_genes_asd=list()
  
  for(i in c(1:length(idx))){
    if (i%%50 == 0) {print(paste(i,"/1000",sep=""))}
    load(paste("data_user/03_TRI/03_01_A_01_permutation/TRI_list_permutation",as.character(idx[i]),".RData",sep=""))
    names(cp_list) <- gsub(" - ","_v_",names(cp_list))
    for(comp in names(cp_list)){
      if(i==1){
        perm_diff_genes_ctl[[comp]] <- list()
        perm_diff_genes_asd[[comp]] <- list()
      }
      tmp=cp_list[[which(names(cp_list)==comp)]]$ASD$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      perm_diff_genes_asd[[comp]][[i]]=rownames(tmp)[which(tmp$FDR < 0.05)]
      tmp=cp_list[[which(names(cp_list)==comp)]]$CTL$data
      tmp$FDR=p.adjust(tmp$PVal,method="fdr")
      perm_diff_genes_ctl[[comp]][[i]]=rownames(tmp)[which(tmp$FDR < 0.05)]
    }
  }
  
  save(perm_diff_genes_ctl,perm_diff_genes_asd,
       file=paste("03_01_A_03_permutation_genes/TRI_permutation_genes_Perm",iter_idx,".RData",sep=""))

}
  
### now compile based on comparison
  
perm_diff_genes_asd_comp <- list()
perm_diff_genes_ctl_comp <- list()

for(j in c(1:10)){
  print(j)
  load(paste("03_01_A_03_permutation_genes/TRI_permutation_genes_Perm",j,".RData",sep=""))
  for(comp in names(cp_list)){
    print(comp)
    if(j==1){
      perm_diff_genes_asd_comp[[comp]] <- vector(mode="list",length=10000)
      perm_diff_genes_ctl_comp[[comp]] <- vector(mode="list",length=10000)
    }
    perm_diff_genes_asd_comp[[comp]][idx_list[[j]]] = perm_diff_genes_asd[[comp]]
    perm_diff_genes_ctl_comp[[comp]][idx_list[[j]]] = perm_diff_genes_ctl[[comp]]
  }
}

for(comp in names(perm_diff_genes_asd_comp)){
  perm_diff_genes_asd_comp_single = perm_diff_genes_asd_comp[[comp]]
  perm_diff_genes_ctl_comp_single = perm_diff_genes_ctl_comp[[comp]]
  save(perm_diff_genes_asd_comp_single,perm_diff_genes_ctl_comp_single,
       file=paste("03_01_A_03_permutation_genes/Compiled_TRI_permutation_genes_",comp,".RData",sep=""))
}

  
### Now, count the number of true diff genes in the permutations 
  
for(iter_idx in c(1:10)){
  for(dx in c("CTL","ASD")){

    true_diff_genes_ctl=list()
    true_diff_genes_asd=list()
    
    true_diff_genes_count_ctl=list()
    true_diff_genes_count_asd=list()
    
    load("data_user/03_TRI/03_01_A_01_True_TRI_DGE_CompleteList.RData")
    
    names(cp_list) <- gsub(" - ","_v_",names(cp_list))
    
    comp = names(cp_list)[as.numeric(iter_idx)]
    dx = as.character(dx)  
    
   "Getting True CP Genes"
    
    tmp=cp_list[[which(names(cp_list)==comp)]]$ASD$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    true_diff_genes_asd[[comp]]=rownames(tmp)[which(tmp$FDR < 0.05)]
    tmp=cp_list[[which(names(cp_list)==comp)]]$CTL$data
    tmp$FDR=p.adjust(tmp$PVal,method="fdr")
    true_diff_genes_ctl[[comp]]=rownames(tmp)[which(tmp$FDR < 0.05)]
    
    "Counting Number of 'True Gene' Occurences in Permutations"
    
    true_diff_genes_count_asd[[comp]]=rep(NA,length(true_diff_genes_asd[[comp]]))
    names(true_diff_genes_count_asd[[comp]]) <- true_diff_genes_asd[[comp]]
    true_diff_genes_count_ctl[[comp]]=rep(NA,length(true_diff_genes_ctl[[comp]]))
    names(true_diff_genes_count_ctl[[comp]]) <- true_diff_genes_ctl[[comp]]
    
    load(paste("data_user/03_TRI/03_01_A_03_permutation_genes/Compiled_TRI_permutation_genes_",comp,".RData",sep=""))
    
    if(dx == "ASD"){
      
      print("Counting for ASD")
      
      ind=1
      for(gene in true_diff_genes_asd[[comp]]){
        count=0
        if (ind%%50 == 0) {print(paste(ind,"/",length(true_diff_genes_asd[[comp]]),sep=""))}
        for(i in c(1:length(perm_diff_genes_asd_comp_single))){
          iter=length(grep(gene,perm_diff_genes_asd_comp_single[[i]]))
          count=count+iter
        }
        true_diff_genes_count_asd[[comp]][gene]=count
        ind=ind+1
      }
      
      save(true_diff_genes_count_asd,true_diff_genes_asd,
           file=paste("data_user/03_TRI/03_01_A_03_permutation_genes/TrueGenes_in_PermutedGenes_Compiled_",as.character(comp),"_ASD.RData",sep=""))
      
      
    }else if(dx == "CTL"){
      
      print("Counting for CTL")
      
      ind=1
      for(gene in true_diff_genes_ctl[[comp]]){
        count=0
        if (ind%%50 == 0) {print(paste(ind,"/",length(true_diff_genes_ctl[[comp]]),sep=""))}
        for(i in c(1:length(perm_diff_genes_ctl_comp_single))){
          iter=length(grep(gene,perm_diff_genes_ctl_comp_single[[i]]))
          count=count+iter
        }
        true_diff_genes_count_ctl[[comp]][gene]=count
        ind=ind+1
      }
      
      save(true_diff_genes_count_ctl,true_diff_genes_ctl,
           file=paste("data_user/03_TRI/03_01_A_03_permutation_genes/TrueGenes_in_PermutedGenes_Compiled_",as.character(comp),"_CTL.RData",sep=""))
    } 
  }
}
  
### Combine CTL and ASD data
  
  true_diff_genes_ctl_all <- vector(mode="list",length=55)
  true_diff_genes_asd_all <- vector(mode="list",length=55)
  
  true_diff_genes_count_ctl_all <- vector(mode="list",length=55)
  true_diff_genes_count_asd_all <- vector(mode="list",length=55)
  
  load("data_user/03_TRI/03_01_A_01_True_TRI_DGE_CompleteList.RData")
  
  names(cp_list) <- gsub(" - ","_v_",names(cp_list))
  
  names(true_diff_genes_ctl_all) <- names(true_diff_genes_asd_all) <- names(cp_list)
  names(true_diff_genes_count_ctl_all) <- names(true_diff_genes_count_asd_all) <- names(cp_list)
  
  for(comp in names(cp_list)){
    print(comp)
    for(dx in c("ASD","CTL")){
      print(dx)
      if(dx == "ASD"){
        
        load(paste("data_user/03_TRI/03_01_A_03_permutation_genes/TrueGenes_in_PermutedGenes_Compiled_",comp,"_ASD.RData",sep=""))
        
        true_diff_genes_asd_all[[comp]] <- true_diff_genes_asd[[comp]]        
        true_diff_genes_count_asd_all[[comp]] <- true_diff_genes_count_asd[[comp]]
        
      }else if(dx == "CTL"){
        
        load(paste("data_user/03_TRI/03_01_A_03_permutation_genes/TrueGenes_in_PermutedGenes_Compiled_",comp,"_CTL.RData",sep=""))
        
        true_diff_genes_ctl_all[[comp]] <- true_diff_genes_ctl[[comp]]        
        true_diff_genes_count_ctl_all[[comp]] <- true_diff_genes_count_ctl[[comp]]
        
      }
    }
  }
  
  save(true_diff_genes_count_asd_all,true_diff_genes_count_ctl_all,
       true_diff_genes_asd_all,true_diff_genes_ctl_all,
       file="data_user/03_TRI/03_01_A_03_TrueGenes_in_PermutedGenes_All_Compiled.RData")

### Now, make a histogram for every comparison of the number of times a 'true' gene was in a scrambled permutation

rm(list=ls())

load("data_user/03_TRI/03_01_A_03_TrueGenes_in_PermutedGenes_All_Compiled.RData")
  
pdf(file="plots/03_TRI/03_01_A_03_Permutation-Count-Density.pdf",width=10)

par(mfrow=c(1,2))
for(comp in names(true_diff_genes_asd_all)){
  if(length(true_diff_genes_count_asd_all[[comp]]) > 0){
    plot(density(true_diff_genes_count_asd_all[[comp]],bw=100),main=paste(comp,"; ASD",sep=""))
    abline(v=500,col="red")
  }else{
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste(comp,"; \n No True Diff Genes; \n ASD",sep=""), 
         cex = 1, col = "black")
  }
  if(length(true_diff_genes_count_ctl_all[[comp]]) > 0){
    plot(density(true_diff_genes_count_ctl_all[[comp]],bw=100),main=paste(comp,"; CTL",sep=""))
    abline(v=500,col="red")
  }else{
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, paste(comp,"; \n No True Diff Genes; \n CTL",sep=""), 
         cex = 1, col = "black")
  }
}

dev.off()

### Choose filter of 5% (highly conservative)
### Filer is for a gene to be kept as a 'true' cp gene, from perm test; needs to be present in less than 500 out of the 10,000 perms

true_diff_genes_final_asd=list()
true_diff_genes_final_ctl=list()

for(comp in names(true_diff_genes_asd_all)){
  idx = which(true_diff_genes_count_asd_all[[comp]] < 500)
  true_diff_genes_final_asd[[comp]] = true_diff_genes_asd_all[[comp]][idx]
  idx = which(true_diff_genes_count_ctl_all[[comp]] < 500)
  true_diff_genes_final_ctl[[comp]] = true_diff_genes_ctl_all[[comp]][idx]
}

save(true_diff_genes_final_asd,true_diff_genes_final_ctl,file="data_user/03_TRI/03_01_A_03_All_ARI_Genes.RData")

### Next, get gene expression profiles of ARI genes

### First, load relevant data

ttable = read.csv("data_user/03_TRI/03_01_A_01_Permutation_RegComp-TTable.csv") ## permutation results
load("data_user/03_TRI/03_01_A_01_True_TRI_DGE_ASD_v_Control.RData") ## DGE results for TRI analysis with true samples

### Get reference to match regions to cortical lobules

lobe_match=data.frame("Contrast"=names(diff_list),
                      "Lobe_Contrast"=rep(NA,length(names(diff_list))))

frontal = c("BA9","BA4-6","BA24","BA44-45")
parietal = c("BA3-1-2-5","BA7","BA39-40")
temporal = c("BA38","BA20-37","BA41-42-22")
occipital = c("BA17")

for(comp in lobe_match[,1]){
  reg1=strsplit(comp,split=" - ")[[1]][1]
  reg2=strsplit(comp,split=" - ")[[1]][2]
  if( ((reg1 %in% frontal) & (reg2 %in% parietal)) | ((reg2 %in% frontal) & (reg1 %in% parietal)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Frontal - Parietal"
  }else if( ((reg1 %in% frontal) & (reg2 %in% temporal)) | ((reg2 %in% frontal) & (reg1 %in% temporal)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Frontal - Temporal"
  }else if( ((reg1 %in% frontal) & (reg2 %in% occipital)) | ((reg2 %in% frontal) & (reg1 %in% occipital)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Occipital - Frontal"
  }else if( ((reg1 %in% frontal) & (reg2 %in% frontal)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Frontal - Frontal"
  }else if( ((reg1 %in% temporal) & (reg2 %in% occipital)) | ((reg2 %in% temporal) & (reg1 %in% occipital)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Occipital - Temporal"
  }else if( ((reg1 %in% parietal) & (reg2 %in% occipital)) | ((reg2 %in% parietal) & (reg1 %in% occipital)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Occipital - Parietal"
  }else if( ((reg1 %in% parietal) & (reg2 %in% temporal)) | ((reg2 %in% parietal) & (reg1 %in% temporal)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Parietal - Temporal"
  }else if( ((reg1 %in% parietal) & (reg2 %in% parietal)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Parietal - Parietal"
  }else if( ((reg1 %in% temporal) & (reg2 %in% temporal)) ){
    lobe_match$Lobe_Contrast[which(lobe_match$Contrast==comp)]="Temporal - Temporal"
  }
}

### from this point on - just look at the regional comparisons significantly attenuated in ASD

sig_comps = ttable$Comparison[which(ttable$Perm_PVal < 0.05)]
names(sig_comps) = rep("CTL",length(sig_comps))
sig_comps = gsub(" - ","_v_",sig_comps)

sig_comp_list_ctl = true_diff_genes_final_ctl[which(names(true_diff_genes_final_ctl) %in% sig_comps[which(names(sig_comps)=="CTL")])] 

### remove any genes from the controls that are present in asd

sig_comp_list_ctl_unique = list()

for(comp in names(sig_comp_list_ctl)){
  print(comp)
  rm = true_diff_genes_final_asd[[comp]][which(true_diff_genes_final_asd[[comp]] %in% sig_comp_list_ctl[[comp]])]
  print(length(rm))
  if(length(rm)==0){
    sig_comp_list_ctl_unique[[comp]]=sig_comp_list_ctl[[comp]]
  }else{
    sig_comp_list_ctl_unique[[comp]]=sig_comp_list_ctl[[comp]][-match(rm,sig_comp_list_ctl[[comp]])]
  }
}

sig_comp_list_ctl = sig_comp_list_ctl_unique

names(sig_comp_list_ctl) <- paste(names(sig_comp_list_ctl),"_CTL-gt-ASD",sep="")
sig_comp_list = sig_comp_list_ctl

save(sig_comp_list,file="data_user/03_TRI/03_01_A_03_ARI_genes_from_RegComparisons_Attenuated_inASD.RData")

### determine which region each comparison is expressed in

cp_genes = list()

load("data_provided/03_TRI/03_01_A_AllProcessedData_wModelMatrix.RData")
load("data_provided/03_TRI/03_01_A_RegressedExpression_TRI.RData")
load("data_provided/03_TRI/03_01_A_01_Permutation_Index.RData")

for(comp in names(sig_comp_list)){
    
  asd_att_genes = sig_comp_list[[comp]]
  
  idx_name = gsub("_CTL-gt-ASD","",comp)
  idx_name = gsub("_v_"," - ",idx_name)
  
  reg1 = strsplit(comp,split="_")[[1]][1]
  reg2 = strsplit(comp,split="_")[[1]][3]
  
  asd_att_expr = datExpr.reg[which(rownames(datExpr.reg) %in% asd_att_genes),idx_comp_all[[idx_name]]]
  datMeta_asd_att = datMeta[idx_comp_all[[idx_name]],]
  datMeta_asd_att$Dx_Reg = paste(datMeta_asd_att$Diagnosis,datMeta_asd_att$region,sep="_")
  
  ## sort the genes which are more highly expressed in reg1 or reg2
  
  count=rep(NA,length(asd_att_genes))
  names(count) <- rownames(asd_att_expr)
  for(gene in rownames(asd_att_expr)){
    avg_1=median(asd_att_expr[gene,which(datMeta_asd_att$Dx_Reg==paste("CTL",reg1,sep="_"))])
    avg_2=median(asd_att_expr[gene,which(datMeta_asd_att$Dx_Reg==paste("CTL",reg2,sep="_"))])
    if(avg_1 > avg_2){
      count[gene]=reg1
    }else{
      count[gene]=reg2
    }
  }
  
  reg1_genes = names(count)[which(count==reg1)]
  reg2_genes = names(count)[which(count==reg2)]
  cp_genes[[comp]][[reg1]]=reg1_genes
  cp_genes[[comp]][[reg2]]=reg2_genes
    
}
  

save(cp_genes,file="data_user/03_TRI/03_01_A_03_ARI_genes_from_RegComparisons_Attenuated_inASD_RegSorted.RData")

### Sort ARI genes into relevant groups: increasing in control expression from anterior to posterior, and vice versa
### Obtain two 'interesting groups': BA17 and BA39-40 expressed genes v. all other regions

cp_att_groups = list("BA17_BA39-40"=rep(NA,2),"Other_Cortical_Regions"=rep(NA,2))

for(comp in names(cp_genes)){
  regs = names(cp_genes[[comp]])
  for(reg in regs){
    if(reg == "BA17" | reg == "BA39-40"){
      cp_att_groups[["BA17_BA39-40"]] = c(cp_att_groups[["BA17_BA39-40"]],cp_genes[[comp]][[reg]])
    }else{
      cp_att_groups[["Other_Cortical_Regions"]] = c(cp_att_groups[["Other_Cortical_Regions"]],cp_genes[[comp]][[reg]])
    }
  }
}

cp_att_groups[["BA17_BA39-40"]] = cp_att_groups[["BA17_BA39-40"]][!duplicated(cp_att_groups[["BA17_BA39-40"]])]
cp_att_groups[["BA17_BA39-40"]] = cp_att_groups[["BA17_BA39-40"]][-which(is.na(cp_att_groups[["BA17_BA39-40"]])==TRUE)]

cp_att_groups[["Other_Cortical_Regions"]] = cp_att_groups[["Other_Cortical_Regions"]][!duplicated(cp_att_groups[["Other_Cortical_Regions"]])]
cp_att_groups[["Other_Cortical_Regions"]] = cp_att_groups[["Other_Cortical_Regions"]][-which(is.na(cp_att_groups[["Other_Cortical_Regions"]])==TRUE)]

### remove any intersecting genes between the two groups - these are not the ones we are interested in (likely BA39-40 and another region v. highly anterior region)

inter = intersect(cp_att_groups[["BA17_BA39-40"]],cp_att_groups[["Other_Cortical_Regions"]])
cp_att_groups[["BA17_BA39-40"]] = cp_att_groups[["BA17_BA39-40"]][-which(cp_att_groups[["BA17_BA39-40"]] %in% inter)]
cp_att_groups[["Other_Cortical_Regions"]] = cp_att_groups[["Other_Cortical_Regions"]][-which(cp_att_groups[["Other_Cortical_Regions"]] %in% inter)]

all_att_genes = c(cp_att_groups[[1]],cp_att_groups[[2]])

### finally, double check that median of genes isn't higher in the other group - remove any that don't comply

reg1=names(cp_att_groups)[1]
reg2=names(cp_att_groups)[2]

count=rep(NA,length(all_att_genes))
names(count) <- rownames(all_att_genes)
for(gene in all_att_genes){
  avg_1=median(datExpr.reg[gene,which( datMeta$Diagnosis == "CTL" & (datMeta$region == "BA17" | datMeta$region == "BA39-40" ))])
  avg_2=median(datExpr.reg[gene,which( datMeta$Diagnosis == "CTL" & datMeta$region != "BA17" & datMeta$region != "BA39-40" )])
  if(avg_1 > avg_2){
    count[gene]=reg1
  }else{
    count[gene]=reg2
  }
}

plot <- table(count) 

#count
#BA17_BA39-40 Other_Cortical_Regions 
#2094                   1851 

reg1_genes = names(count)[which(count==reg1)] # 2094
reg2_genes = names(count)[which(count==reg2)] # 1851

length(cp_att_groups[[reg1]]) # 2037
length(cp_att_groups[[reg2]]) # 1908

reg1_keep = intersect(reg1_genes,cp_att_groups[[reg1]]) # 1881
reg2_keep = intersect(reg2_genes,cp_att_groups[[reg2]]) # 1695

### finalize list

cp_att_groups[[reg1]] = reg1_keep
cp_att_groups[[reg2]] = reg2_keep
all_att_genes = c(cp_att_groups[[reg1]],cp_att_groups[[reg2]])

### get DE for these genes in each group in ctls

mod=model.matrix(~ 0 + DxReg+SeqBatch+Sex+Ancestry,data=datMeta_model)
design=data.frame(mod,datMeta_model[,c(6:(dim(datMeta_model)[2]))])

load("data_user/02_DEGenesIsoforms/02_01_A_01_lmFit.RData")

this_contrast_cp_att = makeContrasts(contrast=(DxRegCTL_BA17 + DxRegCTL_BA39_40)/2 -         
                                       (DxRegCTL_BA20_37 + DxRegCTL_BA24 + DxRegCTL_BA3_1_2_5 + DxRegCTL_BA38 +          
                                          DxRegCTL_BA4_6 + DxRegCTL_BA41_42_22 + DxRegCTL_BA44_45 +      
                                          DxRegCTL_BA7 + DxRegCTL_BA9)/9, levels=design)

this_contrast_cp_att <- data.frame(this_contrast_cp_att[match(colnames(fit$coefficients),rownames(this_contrast_cp_att)),])

fit2 = contrasts.fit(fit, this_contrast_cp_att)
fit3= eBayes(fit2,trend = T, robust = T)
tt_cp_att= topTable(fit3, coef=1, number=Inf, sort.by = 'none')

tt_cp_att = tt_cp_att[match(all_att_genes,rownames(tt_cp_att)),]
tt_contrast_cp_att = tt_cp_att[,c(1,5)]

### get the external_gene_names as well and make a .csv with these genes and DGE info.

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host = "grch37.ensembl.org")
# get hg19 annotation

identifier <- "ensembl_gene_id"
getinfo <- c("ensembl_gene_id","external_gene_name","hgnc_symbol","gene_biotype",
             "chromosome_name","strand","band","start_position","end_position","description")

reg1 = rep(names(cp_att_groups)[1],length(cp_att_groups[[1]]))
reg2 = rep(names(cp_att_groups)[2],length(cp_att_groups[[2]]))

tmp = data.frame("ensembl_gene_id"=all_att_genes,"group"=c(reg1,reg2))
tmp = data.frame(cbind(tmp,tt_cp_att[match(tmp$ensembl_gene_id,rownames(tt_cp_att)),c(1,3,4,5)]))
colnames(tmp)[3] = "logFC_BA17_BA39_40_minus_Other_Cortical_Regions"  

geneDat <- getBM(attributes = getinfo, filters=identifier, values = substr(tmp$ensembl_gene_id,1,15),mart=ensembl)
geneDat = geneDat[match(substr(tmp$ensembl_gene_id,1,15),geneDat$ensembl_gene_id),]
tmp = data.frame(cbind(tmp,geneDat[,-1]))

write.csv(tmp,file="data_user/03_TRI/03_01_A_03_Annotated_ARI_Gene_InterestingGroups.csv")
