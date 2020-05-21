#!/bin/bash

##### 01_01_GetCounts.sh
##### Get RNA-seq Counts (Gene-Level and Isoform-Level).
##### This script contains the main pipeline run for all samples.
##### All fastq file variables below are from the same sample.
##### May 2020, Jillian Haney.

batch=
### ASD3Reg or ASDPan

RSEM=/path/to/RSEM

fastq_dir=/path/to/fastq_files

output_dir=/path/to/output

L1fastq1=Sample1_L001_R1.fastq.gz
### Seq Lane 1, Paired Read 1

L1fastq2=Sample1_L001_R2.fastq.gz
### Seq Lane 1, Paired Read 2

L2fastq1=Sample1_L002_R1.fastq.gz
### Seq Lane 2, Paired Read 1

L2fastq2=Sample1_L002_R2.fastq.gz
### Seq Lane 2, Paired Read 2

##### Pipeline #####

### (1) FastQC
### (2) STAR Alignment
### (3) Picard Tools QC
### (4) RSEM Quantification

####################

##### (1) FastQC #####

fastqc --noextract --threads 1 --outdir ${output_dir}/Sample1L1R1_FastQC ${fastq_dir}/${L1fastq1}
fastqc --noextract --threads 1 --outdir ${output_dir}/Sample1L1R2_FastQC ${fastq_dir}/${L1fastq2}
fastqc --noextract --threads 1 --outdir ${output_dir}/Sample1L2R1_FastQC ${fastq_dir}/${L2fastq1}
fastqc --noextract --threads 1 --outdir ${output_dir}/Sample1L2R2_FastQC ${fastq_dir}/${L2fastq2}

##### (2) STAR Alignment #####

STAR      --runMode                 alignReads \
          --genomeDir               GRCh37_Gencode25_STAR_index \
          --outFileNamePrefix       ${output_dir}/Sample1L1 \
          --readFilesCommand        zcat \
          --runThreadN              4 \
          --outMultimapperOrder     "Random" \
          --outSAMmultNmax          1 \
          --outSAMtype              "BAM SortedByCoordinate" \
          --outSAMunmapped          "Within KeepPairs" \
	  --outSAMattrRGline        "ID:Sample1.L1 SM:Sample1 PL:ILLUMINA" \
          --outBAMcompression       6 \
          --bamRemoveDuplicatesType "UniqueIdentical" \
          --outSJfilterReads        "Unique" \
          --readFilesIn             $L1fastq1 $L1fastq2 \
          --quantMode               TranscriptomeSAM \

STAR      --runMode                 alignReads \
          --genomeDir               GRCh37_Gencode25_STAR_index \
          --outFileNamePrefix       ${output_dir}/Sample1L2 \
          --readFilesCommand        zcat \
          --runThreadN              4 \
          --outMultimapperOrder     "Random" \
          --outSAMmultNmax          1 \
          --outSAMtype              "BAM SortedByCoordinate" \
          --outSAMunmapped          "Within KeepPairs" \
          --outSAMattrRGline        "ID:Sample1.L2 SM:Sample1 PL:ILLUMINA" \
          --outBAMcompression       6 \
          --bamRemoveDuplicatesType "UniqueIdentical" \
          --outSJfilterReads        "Unique" \
          --readFilesIn             $L2fastq1 $L2fastq2 \
          --quantMode               TranscriptomeSAM \

### Merge BAMS for Picard Tools QC

cd $output_dir

samtools merge -l 6 -@ 1 -rfp Sample1.bam Sample1L1.bam Sample1L2.bam

##### (3) Picard Tools QC #####

cd $output_dir

java -Xmx4g -Xms64m -jar picard.jar CollectAlignmentSummaryMetrics \
  R=GRCh37.fa \
  O=Sample1.alignment_metrics.txt \
  EXPECTED_PAIR_ORIENTATIONS=FR \
  EXPECTED_PAIR_ORIENTATIONS=RF \
  I=Sample1.bam

if [ "$batch" = ASD3Reg ]; then

java -Xmx4g -Xms64m -jar picard.jar CollectRnaSeqMetrics \
  I=Sample1.bam \
  O=Sample1.rna_seq_metrics.txt \
  R=GRCh37.fa \
  STRAND=NONE \
  RIBOSOMAL_INTERVALS=GRCh37_ribosomal_exons.interval_list \
  REF_FLAT=GRCh37.annotation.flat.txt

fi

if [ "$batch" = ASDPan ]; then

java -Xmx4g -Xms64m -jar picard.jar CollectRnaSeqMetrics \
  I=Sample1.bam \
  O=Sample1.rna_seq_metrics.txt \
  R=GRCh37.fa \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=GRCh37_ribosomal_exons.interval_list \
  REF_FLAT=GRCh37.annotation.flat.txt

fi

java -Xmx4g -Xms64m -jar picard.jar CollectGcBiasMetrics \
  R=GRCh37.fa \
  I=Sample1.bam \
  O=Sample1.gc_bias_metrics.txt \
  S=Sample1.gc_bias_metrics.summary.txt \
  CHART=Sample1.gc_bias.chart.pdf

java -Xmx4g -Xms64m -jar picard.jar MarkDuplicates \
  I=Sample1.bam \
  O=Sample1.temp.dup.bam \
  M=Sample1.duplication_metrics.txt

java -Xmx4g -Xms64m -jar picard.jar CollectInsertSizeMetrics \
  I=Sample1.bam \
  O=Sample1.insert.size.metrics.txt \
  H=Sample1.insert.size.histogram.pdf

### Remove excess BAM files to save disk space

rm Sample1.temp.dup.bam
rm Sample.bam

##### (4) RSEM Quantification #####

cd output_dir

samtools merge -l 6 -@ 1 -rfp Sample1.tx.bam Sample1L1.bam Sample1L2.bam

${RSEM}/convert-sam-for-rsem -p 1 Sampl1.tx.bam Sample1

if [ "$batch" = ASD3Reg ]; then

${RSEM}/rsem-calculate-expression --paired-end \
    --strandedness none \
    --num-threads 1 \
    --alignments \
    --estimate-rspd \
    --no-bam-output \
    --seed 20161010 \
    Sample1.bam \
    GRCh37.RSEM.index \
    Sample1

fi

if [ "$batch" = ASDPan ]; then

${RSEM}/rsem-calculate-expression --paired-end \
    --strandedness reverse \
    --num-threads 1 \
    --alignments \
    --estimate-rspd \
    --no-bam-output \
    --seed 20161010 \
    Sample1.bam \
    GRCh37.RSEM.index \
    Sample1

fi

### Remove excess BAM files to save disk space

rm Sample1.bam
rm Sample1.tx.bam
