# Sequence Alignment Theory

## Learning Objectives

-
-
-

## Content

Convert slides 

## `bwa`

This `bwa` step will be done prior to the course.

### Normal

```
module load gcc/6.2.0 bwa/0.7.17 samtools/1.15.1 

bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-normal\tPL:illumina\tPU:syn3-normal\tSM:syn3-normal' \
/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
fastq_files/synthetic_challenge_set3_normal_NGv3_1.fq.gz \
fastq_files/synthetic_challenge_set3_normal_NGv3_2.fq.gz \
-o alignments/normal_GRCh38.p7.sam

samtools sort \
-@ 8 \
-O bam \
-o alignments/normal_sorted_GRCh38.p7.bam \
alignments/normal_GRCh38.p7.sam

samtools index \
alignments/normal_sorted_GRCh38.p7.bam
```

### Tumor

```
module load gcc/6.2.0 bwa/0.7.17 samtools/1.15.1 

bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-tumor\tPL:illumina\tPU:syn3-tumor\tSM:syn3-tumor' \
/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
fastq_files/synthetic_challenge_set3_tumor_NGv3_1.fq.gz \
fastq_files/synthetic_challenge_set3_tumor_NGv3_2.fq.gz \
-o alignments/tumor_GRCh38.p7.sam

samtools sort \
-@ 8 \
-O bam \
-o alignments/tumor_sorted_GRCh38.p7.bam \
alignments/tumor_GRCh38.p7.sam

samtools index \
alignments/tumor_sorted_GRCh38.p7.bam
```
