##### 05_01_A_RegVarAnalysis.R
##### Analyze regional variation in ASD gene expression effect (Gene-Level)
##### November 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="/path/to/my/directory"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

library(gplots); library(reshape2); library(gridExtra); library(limma)
library(grid); library(ggplot2); library(lattice); library(reshape2)
library(igraph); library(scales); library(biomaRt); library(stringi); library(readxl)

##### Pipeline #####

###(1) Macaque NeuN Density ~ Module Eigengene Region-specific ASD Effect
###(2) Estimated Layer 3/4 Cortical Thickness ~ Module Eigengene Region-specific ASD Effect

##### (1) Macaque NeuN Density ~ Module Eigengene Region-specific ASD Effect #####

load("data_provided/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData") ### produced by 04_01_A_WGCNA.R, section 3
load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8

load("data_provided/05_RegVarAnalysis/05_01_A_01_ME_assoc_table.RData")
load("data_provided/05_RegVarAnalysis/05_01_A_01_ME_assoc_table_all_covariates.RData")
### both of theses objects are extracted from 04_01_A_03_ME_Diagnosis_by_Region_Effects.RData, produced by 04_01_A_WGCNA.R, section 3

### Get Neuronal Density Associations with MEs

datMeta_model$Region = datMeta$region
datMeta_model$Diagnosis = datMeta$Diagnosis

ME_dysreg = as.matrix(assoc_table_list$ASDReg_B)
ME_dysreg_orig = ME_dysreg

ME_dysreg = ME_dysreg_orig

macaque_dens = read.csv("data_provided/05_RegVarAnalysis/collins_PNAS2010_macaque_cortical_density.csv")
macaque_dens_wt = macaque_dens[,c(1,5)]

macaque_dens_wt[,2] = gsub(",","",macaque_dens_wt[,2])

s1_dens = mean(as.numeric(macaque_dens_wt[c(2,3),2]))
s1_dens = c("BA3_1_2_5",s1_dens)
macaque_dens_wt = rbind(macaque_dens_wt,s1_dens)

macaque_dens_wt = macaque_dens_wt[-c(2,3),]

colnames(macaque_dens_wt) <-  c("Region","Density")

keep_regs = gsub("-","_",macaque_dens_wt$Region)
idx = which(colnames(ME_dysreg) %in% keep_regs)
ME_dysreg = ME_dysreg[,idx]

reg_neu_density = macaque_dens_wt$Density[match(colnames(ME_dysreg),gsub("-","_",macaque_dens_wt$Region))]
reg_neu_density = gsub(",","",reg_neu_density)
reg_neu_density = as.numeric(reg_neu_density)

mac_dens = data.frame(matrix(NA,nrow(ME_dysreg),ncol=2))
colnames(mac_dens) = c("P-Value","B")
rownames(mac_dens) = rownames(ME_dysreg)

for(i in c(1:nrow(ME_dysreg))){
  p=coef(summary(lm(ME_dysreg[i,] ~ reg_neu_density)))[2,4]
  mac_dens$`P-Value`[which(rownames(mac_dens)==rownames(ME_dysreg)[i])]=p
  b=coef(summary(lm(ME_dysreg[i,] ~ reg_neu_density)))[2,1]
  mac_dens$`B`[which(rownames(mac_dens)==rownames(ME_dysreg)[i])]=b
}

mac_dens$FDR = p.adjust(mac_dens$`P-Value`,method="fdr")

pdf(file="plots/04_WGCNA/05_01_A_01_Macaque_Neuronal_Density_v_ME_PHist.pdf")
hist(mac_dens$`P-Value`,main="Histogram of P-value\nME logFC ~ Mac_Neuron_Density")
dev.off()

rownames(mac_dens)[which(mac_dens$FDR < 0.1)]

#[1] "ME2_blue"           "ME3_brown"          "ME4_yellow"         "ME6_red"           
#[5] "ME12_tan"           "ME13_salmon"        "ME14_cyan"          "ME16_lightcyan"    
#[9] "ME17_grey60"        "ME21_darkred"       "ME22_darkgreen"     "ME28_saddlebrown"  
#[13] "ME30_paleturquoise" "ME34_plum1"         "ME35_mediumpurple3"

# mostly neurons, but also some glia (particularly early-up mods/developmental)

save(mac_dens,file="data_user/05_RegVarAnalysis/05_01_A_01_Macaque_Neu_Density_v_MEs.RData")

### Plot

MEs_keep = MEs
PC_plot_gg_reg1 = data.frame(cbind(MEs_keep,as.factor(datMeta$Diagnosis),as.factor(gsub("-","_",datMeta$region))))
colnames(PC_plot_gg_reg1)[c((ncol(PC_plot_gg_reg1)-1):ncol(PC_plot_gg_reg1))] = c("Diagnosis","Region")
PC_plot_gg_reg1 = melt(PC_plot_gg_reg1,id = c("Diagnosis","Region"))
colnames(PC_plot_gg_reg1)[3] = "PC1"

PC_plot_gg_reg1 = PC_plot_gg_reg1[-which(PC_plot_gg_reg1$Diagnosis=="Dup15q"),]

reg_order = c("BA9","BA44_45","BA4_6","BA24","BA38","BA20_37","BA41_42_22","BA3_1_2_5","BA7","BA39_40","BA17")

PC_plot_gg_reg1$Region = factor(PC_plot_gg_reg1$Region,levels=reg_order)
PC_plot_gg_reg1$Diagnosis = gsub("CTL","Control",PC_plot_gg_reg1$Diagnosis)
PC_plot_gg_reg1$Diagnosis = factor(PC_plot_gg_reg1$Diagnosis,levels=c("Control","ASD"))

ME_dysreg = as.matrix(assoc_table_list$ASDReg_B)

macaque_dens = read.csv("data_provided/05_RegVarAnalysis/collins_PNAS2010_macaque_cortical_density.csv")
macaque_dens_wt = macaque_dens[,c(1,5)]

macaque_dens_wt[,2] = gsub(",","",macaque_dens_wt[,2])

s1_dens = mean(as.numeric(macaque_dens_wt[c(2,3),2]))
s1_dens = c("BA3_1_2_5",s1_dens)
macaque_dens_wt = rbind(macaque_dens_wt,s1_dens)

macaque_dens_wt = macaque_dens_wt[-c(2,3),]

colnames(macaque_dens_wt) <-  c("Region","Density")

keep_regs = gsub("-","_",macaque_dens_wt$Region)
idx = which(colnames(ME_dysreg) %in% keep_regs)
ME_dysreg = ME_dysreg[,idx]

reg_neu_density = macaque_dens_wt$Density[match(colnames(ME_dysreg),gsub("-","_",macaque_dens_wt$Region))]
reg_neu_density = gsub(",","",reg_neu_density)
reg_neu_density = as.numeric(reg_neu_density)

load("data_user/05_RegVarAnalysis/05_01_A_01_Macaque_Neu_Density_v_MEs.RData")

plotDat = data.frame(matrix(NA,nrow=1,ncol=4))
colnames(plotDat) = c("logFC","Density","Region","Module")

p_vec <- fdr_vec <- rep(NA,ncol(MEs))
names(p_vec) <- names(fdr_vec) <- ncol(MEs)

for(mod in rownames(ME_dysreg)){
  
  logFC = ME_dysreg[which(rownames(ME_dysreg)==mod),]
  names(logFC) = colnames(ME_dysreg)
  
  fdr_vec[mod]=mac_dens$FDR[which(rownames(mac_dens)==mod)]
  p_vec[mod]=mac_dens$`P-Value`[which(rownames(mac_dens)==mod)]
  
  tmp_add = data.frame("logFC"=logFC,
                       "Density"=reg_neu_density,
                       "Region"=colnames(ME_dysreg),
                       "Module"=rep(mod,length(logFC)))
  
  plotDat = data.frame(rbind(plotDat,tmp_add))
  
}

plotDat = plotDat[-1,]
plotDat$Region = gsub("_","-",plotDat$Region)

show_col(hue_pal()(9))

discrete_cols = c("BA9"="#AEA200","BA24"="#AEA200","BA4-6"="#AEA200","BA3-1-2-5"="#DB8E00","BA41-42-22"="#00BD5C","BA17"="#B385FF")
discrete_shapes = c("BA9"=16,"BA24"=15,"BA4-6"=17,"BA3-1-2-5"=16,"BA41-42-22"=16,"BA17"=16)
plotDat$Region = factor(plotDat$Region,levels=names(discrete_shapes))

int_mods = rownames(ME_dysreg)
plot_names = paste(int_mods,"\nfdr=",signif(fdr_vec[match(int_mods,names(fdr_vec))],2),
                   ", p=",signif(p_vec[match(int_mods,names(p_vec))],2),sep="")
names(plot_names) = int_mods

plotDat$Module = factor(plotDat$Module,levels=int_mods)

asd_col="#FC8D62"
ctl_col="#66C2A5"

text_size=6
tick_size=6
lab_size=6
title_size=6
leg_size=6

pdf(file="plots/05_RegVarAnalysis/05_01_A_01_MacaqueNeuN_v_allModules_ASD_Effect.pdf",width=8,height=11.5)

nd_plot <- ggplot(plotDat, aes(x=Density, y=logFC)) + 
  facet_wrap(. ~ Module,nrow=8,scales="free_y",labeller=labeller(Module=plot_names)) +
  theme_bw() +
  geom_point(size=1,aes(color=Region,shape=Region)) +
  xlab("Macaque NeuN Density (cells/gram)") +
  ylab("ASD Effect") +
  stat_smooth(linetype=1,method="lm",col="black",size=0.25,se = FALSE,alpha=0.5) +
  scale_color_manual(values=discrete_cols,labels=names(discrete_cols)) +
  scale_shape_manual(values=discrete_shapes,labels=names(discrete_shapes)) +
  scale_y_continuous(expand=expansion(0.25)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(plot.title = element_text(hjust=0.5,size=tick_size),
        axis.text.y = element_text(size=tick_size),
        axis.text.x = element_text(size=tick_size,angle=45,hjust=1),
        axis.title = element_text(size=tick_size),
        legend.text = element_text(size=leg_size),
        legend.title = element_text(size=leg_size),
        legend.margin = margin(r=0,l=0,t=0,b=0),
        legend.position="right",
        legend.direction="vertical",
        strip.text.x = element_text(size = tick_size),
        strip.background = element_rect(fill="white"),
        panel.spacing = unit(0.25, "lines"))

print(nd_plot)

dev.off()

### Leave-one-out Cross-validation

datMeta_model$Region = datMeta$region
datMeta_model$Diagnosis = datMeta$Diagnosis

ME_dysreg = as.matrix(assoc_table_list$ASDReg_B)
ME_dysreg_orig = ME_dysreg

### macaque NeuN

ME_dysreg = ME_dysreg_orig

macaque_dens = read.csv("data_provided/05_RegVarAnalysis/collins_PNAS2010_macaque_cortical_density.csv")
macaque_dens_wt = macaque_dens[,c(1,5)]

macaque_dens_wt[,2] = gsub(",","",macaque_dens_wt[,2])

s1_dens = mean(as.numeric(macaque_dens_wt[c(2,3),2]))
s1_dens = c("BA3_1_2_5",s1_dens)
macaque_dens_wt = rbind(macaque_dens_wt,s1_dens)

macaque_dens_wt = macaque_dens_wt[-c(2,3),]

colnames(macaque_dens_wt) <-  c("Region","Density")

keep_regs = gsub("-","_",macaque_dens_wt$Region)
idx = which(colnames(ME_dysreg) %in% keep_regs)
ME_dysreg = ME_dysreg[,idx]

reg_neu_density = macaque_dens_wt$Density[match(colnames(ME_dysreg),gsub("-","_",macaque_dens_wt$Region))]
reg_neu_density = gsub(",","",reg_neu_density)
reg_neu_density = as.numeric(reg_neu_density)
names(reg_neu_density) = colnames(ME_dysreg)

mac_dens = data.frame(matrix(NA,nrow(ME_dysreg),ncol=12))
colnames(mac_dens) = c(paste("LeaveOut",colnames(ME_dysreg),"MacNeuN_Beta",sep="_"),
                       paste("LeaveOut",colnames(ME_dysreg),"MacNeuN_Pvalue",sep="_"))
rownames(mac_dens) = rownames(ME_dysreg)

for(i in c(1:nrow(ME_dysreg))){
  for(j in c(1:ncol(ME_dysreg))){
    
    reg_out = colnames(ME_dysreg)[j]
    tmp_ME = ME_dysreg[,-j]
    tmp_reg_neu = reg_neu_density[-which(names(reg_neu_density)==reg_out)]
    
    p=coef(summary(lm(tmp_ME[i,] ~ tmp_reg_neu)))[2,4]
    mac_dens[which(rownames(mac_dens)==rownames(ME_dysreg)[i]),
             which(colnames(mac_dens)==paste("LeaveOut",reg_out,"MacNeuN_Pvalue",sep="_"))]=p
    b=coef(summary(lm(tmp_ME[i,] ~ tmp_reg_neu)))[2,1]
    mac_dens[which(rownames(mac_dens)==rownames(ME_dysreg)[i]),
             which(colnames(mac_dens)==paste("LeaveOut",reg_out,"MacNeuN_Beta",sep="_"))]=b
  }
}

for(j in c(1:ncol(ME_dysreg))){
  reg_out = colnames(ME_dysreg)[j]
  idx = which(colnames(mac_dens)==paste("LeaveOut",reg_out,"MacNeuN_Pvalue",sep="_"))
  fdr = p.adjust(mac_dens[,idx],method="fdr")
  mac_dens = data.frame(cbind(mac_dens,fdr))
  colnames(mac_dens)[ncol(mac_dens)] = paste("LeaveOut",reg_out,"MacNeuN_FDR",sep="_")
}

## reorder

for(j in c(1:ncol(ME_dysreg))){
  reg_out = colnames(ME_dysreg)[j]
  idx = grep(reg_out,colnames(mac_dens))
  if(j==1){
    mac_dens_lo = mac_dens[,idx]
  }else{
    mac_dens_lo = data.frame(cbind(mac_dens_lo,mac_dens[,idx]))
  }
}

neu_comp = mac_dens
colnames(neu_comp) = c("MacNeuN_Pvalue","MacNeuN_Beta","MacNeuN_FDR")
neu_comp = neu_comp[,c(2,1,3)]

neu_comp = data.frame(cbind(neu_comp,mac_dens_lo[match(rownames(neu_comp),rownames(mac_dens_lo)),]))

write.csv(neu_comp,file="data_user/05_RegVarAnalysis/NeuN_Density_v_MEs.csv")

##### (2) Estimated Layer 3/4 Cortical Thickness ~ Module Eigengene Region-specific ASD Effect #####

load("data_provided/04_WGCNA/04_01_A_03_rWGCNA_consensusModules.RData") ### produced by 04_01_A_WGCNA.R, section 3
load("data_provided/04_WGCNA/04_01_A_AllProcessedData_wModelMatrix.RData") ### produced by 01_02_A_CountsProcessing.R, section 6
load("data_provided/04_WGCNA/04_01_A_RegressedExpression.RData") ### produced by 01_02_A_CountsProcessing.R, section 8

load("data_provided/05_RegVarAnalysis/05_01_A_01_ME_assoc_table.RData")
load("data_provided/05_RegVarAnalysis/05_01_A_01_ME_assoc_table_all_covariates.RData")
### both of theses objects are extracted from 04_01_A_03_ME_Diagnosis_by_Region_Effects.RData, produced by 04_01_A_WGCNA.R, section 3

L4 = read.csv("data_provided/05_RegVarAnalysis/wagstyl_PlosBiol2020_cortical_L4_thickness.csv")
idx = which(L4$BA=="")
L4 = L4[-idx,]

L4_mat = data.frame(matrix(NA,nrow=11,ncol=3))
colnames(L4_mat) = c("Region","VonEconomo_L4","BigBrain_L4")
L4_mat$Region = unique(L4$BA)

for(reg in L4_mat$Region){
  idx = which(L4$BA==reg)
  if(length(idx) > 1){
    L4_mat$VonEconomo_L4[which(L4_mat$Region==reg)] = mean(L4$ve_4[idx])
    L4_mat$BigBrain_L4[which(L4_mat$Region==reg)] = mean(L4$bigbrain_layer_4[idx])
  }else{
    L4_mat$VonEconomo_L4[which(L4_mat$Region==reg)] = L4$ve_4[idx]
    L4_mat$BigBrain_L4[which(L4_mat$Region==reg)] = L4$bigbrain_layer_4[idx]
  }
}

datMeta_model$Region = datMeta$region
datMeta_model$Diagnosis = datMeta$Diagnosis

ME_dysreg = as.matrix(assoc_table_list$ASDReg_B)
ME_dysreg_orig = ME_dysreg
ME_dysreg = ME_dysreg_orig

ME_dysreg = ME_dysreg[,match(L4_mat$Region,colnames(ME_dysreg))]

ME_L4_assoc = data.frame(matrix(NA,nrow=nrow(ME_dysreg),ncol=4))
colnames(ME_L4_assoc) = c("VonEconomo_L4_Beta","VonEconomo_L4_Pvalue",
                          "BigBrain_L4_Beta","BigBrain_L4_Pvalue")
rownames(ME_L4_assoc) = rownames(ME_dysreg)

pdf(file="plots/05_RegVarAnalysis/05_01_A_02_L4_Thickness_v_allModules_ASD_Effect.pdf",width=12,height=6)

for(i in c(1:nrow(ME_dysreg))){
  
  ME_L4_assoc$VonEconomo_L4_Pvalue[i]=coef(summary(lm(ME_dysreg[i,] ~ L4_mat$VonEconomo_L4)))[2,4]
  ME_L4_assoc$VonEconomo_L4_Beta[i]=coef(summary(lm(ME_dysreg[i,] ~ L4_mat$VonEconomo_L4)))[2,1]
  ME_L4_assoc$BigBrain_L4_Pvalue[i]=coef(summary(lm(ME_dysreg[i,] ~ L4_mat$BigBrain_L4)))[2,4]
  ME_L4_assoc$BigBrain_L4_Beta[i]=coef(summary(lm(ME_dysreg[i,] ~ L4_mat$BigBrain_L4)))[2,1]
  
  ME_plot = data.frame("ME"=as.numeric(ME_dysreg[i,]),"vonEconomo_L4"=L4_mat$VonEconomo_L4,
                       "BigBrain_L4"=L4_mat$BigBrain_L4,
                       "Region"=factor(colnames(ME_dysreg),
                                       levels=c("BA9","BA44_45","BA4_6","BA24","BA38","BA20_37",
                                                "BA41_42_22","BA3_1_2_5","BA7","BA39_40","BA17")))
  plot1 <- ggplot(ME_plot, aes(x=vonEconomo_L4, y=ME,color=Region)) +
    geom_point(size=2, shape=16) +
    ggtitle(paste(rownames(ME_dysreg)[i]," ~ L4 Thickness (VE)\np=",
                  signif(ME_L4_assoc$VonEconomo_L4_Pvalue[i],2),"; b=",
                  signif(ME_L4_assoc$VonEconomo_L4_Beta[i],2),sep="")) +
    stat_smooth(method="lm",col="black",size=1,alpha=0.5,se = FALSE) +
    theme(plot.title = element_text(hjust=0.5,size=16),
          axis.text.y = element_text(size=14),
          axis.text.x = element_text(size=14,angle=45,hjust=1),
          axis.title = element_text(size=16),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16))
  
  plot2 <- ggplot(ME_plot, aes(x=BigBrain_L4, y=ME,color=Region)) +
    geom_point(size=2, shape=16) +
    ggtitle(paste(rownames(ME_dysreg)[i]," ~ L4 Thickness (BB)\np=",
                  signif(ME_L4_assoc$BigBrain_L4_Pvalue[i],2),"; b=",
                  signif(ME_L4_assoc$BigBrain_L4_Beta[i],2),sep="")) +
    stat_smooth(method="lm",col="black",size=1,alpha=0.5,se = FALSE) +
    theme(plot.title = element_text(hjust=0.5,size=16),
          axis.text.y = element_text(size=14),
          axis.text.x = element_text(size=14,angle=45,hjust=1),
          axis.title = element_text(size=16),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16))
  
  grid.arrange(plot1,plot2,nrow=1)
  
}

dev.off()

ME_L4_assoc$VonEconomo_L4_FDR = p.adjust(ME_L4_assoc$VonEconomo_L4_Pvalue,method="fdr")
ME_L4_assoc$BigBrain_L4_FDR = p.adjust(ME_L4_assoc$BigBrain_L4_Pvalue,method="fdr")

write.csv(ME_L4_assoc,file="data_user/05_RegVarAnalysis/05_01_A_02_L4_Thickness_v_allModules_ASD_Effect.csv")
