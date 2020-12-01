## Broad transcriptomic dysregulation across the cerebral cortex in ASD

Code, data, and plots accompanying the manuscript.

To visualize gene and isoform expression data across diagnoses, regions, and sex, follow this link.

The 'main_datasets' directory includes main datasets and functional analysis summary plots of all gene and isoform modules.

### File Notation:

PipelineStage_ScriptNumber_[OPTIONAL]Dataset_ScriptSection.extension

### Pipeline Stages:

01) RNAseq Processing
02) DE Genes and Isoforms
03) Transcriptomic Regional Identity (TRI) Analysis
04) Weighted Gene Co-expression Network Analysis (WGCNA)
05) Analysis of Regional Variation

### Script Number:

Refers to the script within the pipeline stage. For example, for 'RNA Processing', Script #1 contains FASTQ file processing steps, whereas Script #2 contain quantified counts processing.

### Dataset:
<ol>
A. Gene-level Expression<br />
B. Isoform-level Expression
</ol>
  
If no dataset is indicated, the object applies to both datasets.

### Script Section:

Section of the script where an object came from is indicated.


