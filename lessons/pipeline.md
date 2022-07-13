# Commands used for Variant Calling materials

## `bwa`

This `bwa` step will be done prior to the course.

### Normal

```
module load gcc/6.2.0 bwa/0.7.17 samtools/1.15.1 

bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-normal\tPL:illumina\tPU:syn3-normal\tSM:syn3-normal' \
fastq_files/synthetic_challenge_set3_normal_NGv3_1.fq.gz \
fastq_files/synthetic_challenge_set3_normal_NGv3_2.fq.gz \
-o alignments/normal_hg19.sam

samtools sort \
-@ 8 \
-O bam \
-o alignments/normal_sorted_hg19.bam \
alignments/normal_hg19.sam

samtools index \
alignments/normal_sorted_hg19.bam
```

### Tumor

```
module load gcc/6.2.0 bwa/0.7.17 samtools/1.15.1 

bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-tumor\tPL:illumina\tPU:syn3-tumor\tSM:syn3-tumor' \
/n/groups/shared_databases/bwa_indexes/hg19 \
fastq_files/synthetic_challenge_set3_tumor_NGv3_1.fq.gz \
fastq_files/synthetic_challenge_set3_tumor_NGv3_2.fq.gz \
-o alignments/tumor_hg19.sam

samtools sort \
-@ 8 \
-O bam \
-o alignments/tumor_sorted_hg19.bam \
alignments/tumor_hg19.sam

samtools index \
alignments/tumor_sorted_hg19.bam
```

## `Mutect2`

Workshop starts here:

Call variants using Mutect2 then filter:

```
module load gatk/4.1.9.0

gatk Mutect2 \
--sequence-dictionary /n/groups/shared_databases/picard_dictionaries/hg19.dict \
-R /n/groups/shared_databases/genomes/hg19.fa \
-I alignments/tumor_sorted_hg19.bam \
--tumor-sample syn3-tumor \
-I alignments/normal_sorted_hg19.bam \
--normal-sample syn3-normal \
--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
-O vcf_files/syn3_hg19-raw.vcf.gz


gatk FilterMutectCalls \
--reference /n/groups/shared_databases/genomes/hg19.fa \
--variant vcf_files/syn3_hg19-raw.vcf.gz \
--output vcf_files/syn3_hg19-raw-filt.vcf.gz
```

## `SnpEff`


Annotate variants with SnpEff:

```
module load snpEff/4.3g

java -jar $SNPEFF/snpEff.jar  eff \
-dataDir /n/groups/shared_databases/snpEff.data/ \
-cancer \
-noLog \
-csvStats vcf_files/syn3_hg19-effects-stats.csv \
-s vcf_files/syn3_hg19-effects-stats.html \
hg19 \
vcf_files/syn3_hg19-raw-filt.vcf.gz > vcf_files/syn3_hg19-raw-filt.snpeff.vcf
```

## `bedtools`

Filter out variants for those in the LCR:

Where does this LCR file come from?

```
module load gcc/6.2.0 bedtools/2.27.1

bedtools intersect \
-v \
-a vcf_files/syn3_hg19-raw-filt.snpeff.vcf \
-b LCR_hg19.bed > vcf_files/syn3_hg19-LCR-filt.snpeff.vcf
```

## `SnpSift`

Example of extracting missense mutations for Oncoprint creation:

```
module load snpEff/4.3g

cat vcf_files/syn3_hg19-LCR-filt.snpeff.vcf | \
java -jar $SNPEFF/SnpSift.jar \
filter \
"( ANN[*].EFFECT = 'missense_variant' )" | \
$SNPEFF/scripts/vcfEffOnePerLine.pl | \
java -jar $SNPEFF/SnpSift.jar \
extractFields \
-  \
"ANN[*].GENE" "EFF[*].AA" "ANN[*].EFFECT"  | \
sort | \
uniq | \
awk '$3 == "missense_variant"' | \
sed 's/missense_variant/MISSENSE/g' | \
awk '{print "sample_1","\t",$0}' > vcf_files/syn3_hg19.missense_variants.txt
```

## TCGA filtering

Filter out TCGA file from cBioPortal data_mutations.txt

```
awk '{print $18,"\t",$1,"\t",$39,"\t"$9}' brca_tcga_pan_can_atlas_2018/data_mutations.txt | \
sed '1d' > brca_filtered.txt
```

