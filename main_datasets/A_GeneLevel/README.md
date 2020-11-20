# README

## File Descriptions

### AllGenesByModule_wDGE_wAnnotation.csv

All genes, with DE gene statisticts and annotation


### Gene_NormalizedExpression_Metadata_wModelMatrix.RData

Gene expression data used for DE gene analysis. log2(CPM), with adjustment for size factors (TMM method), gene length, and GC content.
Also included metadata, sequencing statistics, and the metadata used for the DE gene linear mixed model.


### Gene_NormalizedExpression_TechnicalCovariatesRemoved.RData

Regressed gene expression dataset used for WGCNA. Effects of technical covariates (eg. batch, sequencing statistics) are removed.


### GeneModule_Annotation.csv

Annotation and extra information for every gene module, including covariate effects.


### GeneModule_SummaryPlots.pdf

Summary plots with covariate associations and functional enrichments for every gene module.


### RawData_GeneCounts.RData

Gene-level counts, estimated gene length, metadata, and sequencing statistics (including outliers that were ultimately removed before DE gene analysis and WGCNA).
