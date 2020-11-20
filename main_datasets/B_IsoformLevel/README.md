# README

## File Descriptions

### AllIsoformsByModule_wDGE_wAnnotation.csv

All isoforms, with DE isoform statisticts and annotation


### Isoform_NormalizedExpression_Metadata_wModelMatrix.RData

Estimated isoform expression data used for DE gene analysis. log2(CPM), with adjustment for size factors (TMM method), isoform length, and GC content.
Also included metadata, sequencing statistics, and the metadata used for the DE isoform linear mixed model.


### Isoform_NormalizedExpression_TechnicalCovariatesRemoved.RData

Regressed estimated isoform expression dataset used for WGCNA. Effects of technical covariates (eg. batch, sequencing statistics) are removed.


### IsoformModule_Annotation.csv

Annotation and extra information for every isoform module, including covariate effects.
More information is provided for unique isoform modules (distinct from gene modules).


### Unique_IsoformModule_SummaryPlots.pdf

Summary plots with covariate associations and functional enrichments for every unique isoform module (distinct from gene modules).


### RawData_IsoformCounts.RData

Estimated isoform-level counts, estimated isoform length, metadata, and sequencing statistics (including outliers that were ultimately removed before DE isoform analysis and WGCNA).
