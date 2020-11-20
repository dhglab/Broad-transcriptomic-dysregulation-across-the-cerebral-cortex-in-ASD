##### 04_01_A_RegVarAnalysis.R
##### Analyze regional variation in ASD gene expression effect (Gene-Level)
##### November 2020, Jillian Haney

options(stringsAsFactors = FALSE)

wkdir="C:/Users/jillh/Dropbox/GitHub/"
setwd(paste(wkdir,"Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD/",sep=""))

##### Pipeline #####

###(1) Macaque NeuN Density ~ Module Eigengene Region-specific ASD Effect
###(2) Estimated Layer 3/4 Cortical Thickness ~ Module Eigengene Region-specific ASD Effect
