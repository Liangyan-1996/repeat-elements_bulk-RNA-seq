#!/bin/sh

fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
STAR="/public/home/liunangroup/liangyan/software/miniconda3/envs/rnaseq/bin/STAR"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/bin/samtools"
TEcount="/public/home/liunangroup/liangyan/software/miniconda3/bin/TEcount"

STAR_index="/public/home/liunangroup/liangyan/Genome/ensembl/hg38/STAR_index"
ensembl_gtf="/public/home/liunangroup/liangyan/project/bcl11a/20250721.znf_RNA-seq/hg38_ensembl.gtf"
rmsk_gtf="/public/home/liunangroup/liangyan/project/bcl11a/20250721.znf_RNA-seq/hg38_rmsk_TE.gtf"

sample="$1"
thread="$2"

mkdir QC
mkdir mapping
mkdir log
mkdir TEcount

# qc
$fastp \
 -i raw/${sample}_1.fastq.gz \
 -I raw/${sample}_2.fastq.gz \
 -o QC/${sample}_1.fq.gz \
 -O QC/${sample}_2.fq.gz \
 -w $thread 2> log/$sample.log

# alignment
$STAR \
 --runThreadN $thread \
 --genomeDir $STAR_index \
 --readFilesCommand zcat \
 --readFilesIn QC/${sample}_1.fq.gz QC/${sample}_2.fq.gz \
 --winAnchorMultimapNmax 100 \
 --outFilterMultimapNmax 100 \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.1 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \
 --outFilterType BySJout \
 --outFilterScoreMinOverLread 0.33 \
 --outFilterMatchNminOverLread 0.33 \
 --limitSjdbInsertNsj 1200000 \
 --outSAMstrandField intronMotif \
 --outSAMtype BAM Unsorted \
 --quantMode GeneCounts \
 --outFileNamePrefix mapping/${sample}.

# TEcount
$TEcount \
 -b mapping/${sample}.Aligned.out.bam \
 --GTF ${ensembl_gtf} \
 --TE ${rmsk_gtf} \
 --stranded reverse \
 --project $sample \
 --outdir TEcount