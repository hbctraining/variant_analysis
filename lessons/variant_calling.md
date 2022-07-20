# Variant Calling

## Learning Objectives

-
-
-

## `Mutect2`

Call variants using Mutect2 then filter:

```
module load gatk/4.1.9.0

gatk Mutect2 \
--sequence-dictionary /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.dict \
-R /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
-I alignments/tumor_sorted_GRCh38.p7.bam \
--tumor-sample syn3-tumor \
-I alignments/normal_sorted_GRCh38.p7.bam \
--normal-sample syn3-normal \
--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
-O vcf_files/syn3_GRCh38.p7-raw.vcf.gz


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
