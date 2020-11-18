##### 03_01_A_TRI.R
##### Identify any changes in Transcriptomic Regional Identity (TRI) in ASD (Gene-Level)
##### November 2020, Jillian Haney
##### Note: many variables are named 'cp_xyz'; cp stands for 'cortical patterning', another definition for TRI

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(ggplot2); library(plyr)

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


