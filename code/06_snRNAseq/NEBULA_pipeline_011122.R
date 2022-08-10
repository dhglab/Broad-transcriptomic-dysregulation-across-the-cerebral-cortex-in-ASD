library(nebula)
library(Seurat)

f=list.files("../Seurat_Obj_Subtype/",pattern=".rds")

for(i in 1:length(f)){
  dat.seurat=readRDS(paste0("../Seurat_Obj_Subtype/",f[i]))
  print(f[i])
  print("Working on Parietal...")
  minCell=min(table(dat.seurat@meta.data$BrainRegion))
  cells=rownames(dat.seurat@meta.data[which(dat.seurat@meta.data$BrainRegion=="Parietal"),])
  
  if(length(cells)>minCell){
    cells=cells[sample(length(cells),minCell)]
    print("minCells used")
  }
  dat.s1=subset(dat.seurat,cells=cells)
  
  count=dat.s1@assays$RNA@counts
  meta=dat.s1@meta.data
  meta<-meta[order(meta$Channel),];head(meta)
  meta$nGene=scale(meta$nGene)
  meta$Age=scale(as.integer(meta$Age))
  meta$percent_mito=scale(meta$percent_mito)
  meta$Channel=as.factor(meta$Channel)
  meta$Sex=as.factor(meta$Sex)
  meta$Diagnosis=as.factor(meta$Diagnosis)
  
  count<-count[,rownames(meta)]
  
  df = model.matrix(~Diagnosis+Age+Sex+nGene+percent_mito, data=meta)
  group_cell(count=count,id=meta$Channel,pred=df)
  re = nebula(count,meta$Channel,pred=df,method='LN')
  re$summary$pAdj_DiagnosisCTL=p.adjust(re$summary$p_DiagnosisCTL, method = "BH")
  
  fname=gsub(".rds","",f[i])
  write.csv(re$summary,paste0("DEA_Parietal_", fname,".csv"))
  
  print("Working on Occipital...")
  cells=rownames(dat.seurat@meta.data[which(dat.seurat@meta.data$BrainRegion=="Occipital"),])
  if(length(cells)>minCell){
    cells=cells[sample(length(cells),minCell)]
    print("minCells used")
  }
  dat.s1=subset(dat.seurat,cells=cells)
  
  count=dat.s1@assays$RNA@counts
  meta=dat.s1@meta.data
  meta<-meta[order(meta$Channel),];head(meta)
  meta$nGene=scale(meta$nGene)
  meta$Age=scale(as.integer(meta$Age))
  meta$percent_mito=scale(meta$percent_mito)
  meta$Channel=as.factor(meta$Channel)
  meta$Sex=as.factor(meta$Sex)
  meta$Diagnosis=as.factor(meta$Diagnosis)
  
  count<-count[,rownames(meta)]
  
  df = model.matrix(~Diagnosis+Age+Sex+nGene+percent_mito, data=meta)
  group_cell(count=count,id=meta$Channel,pred=df)
  re = nebula(count,meta$Channel,pred=df,method='LN')
  re$summary$pAdj_DiagnosisCTL=p.adjust(re$summary$p_DiagnosisCTL, method = "BH")

  write.csv(re$summary,paste0("DEA_Occipital_", fname,".csv"))
  
  print("Working on PFC...")
  cells=rownames(dat.seurat@meta.data[which(dat.seurat@meta.data$BrainRegion=="PFC"),])
  if(length(cells)>minCell){
    cells=cells[sample(length(cells),minCell)]
    print("minCells used")
  }
  dat.s1=subset(dat.seurat,cells=cells)
  
  count=dat.s1@assays$RNA@counts
  meta=dat.s1@meta.data
  meta<-meta[order(meta$Channel),];head(meta)
  meta$nGene=scale(meta$nGene)
  meta$Age=scale(as.integer(meta$Age))
  meta$percent_mito=scale(meta$percent_mito)
  meta$Channel=as.factor(meta$Channel)
  meta$Sex=as.factor(meta$Sex)
  meta$Diagnosis=as.factor(meta$Diagnosis)
  
  count<-count[,rownames(meta)]
  
  df = model.matrix(~Diagnosis+Age+Sex+nGene+percent_mito, data=meta)
  group_cell(count=count,id=meta$Channel,pred=df)
  re = nebula(count,meta$Channel,pred=df,method='LN')
  re$summary$pAdj_DiagnosisCTL=p.adjust(re$summary$p_DiagnosisCTL, method = "BH")
  
  write.csv(re$summary,paste0("DEA_PFC_", fname,".csv"))
}
