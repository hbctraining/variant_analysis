# Variant Filtering

## Learning Objectives

-
-
-

```
gatk FilterMutectCalls \
--reference /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
--variant vcf_files/syn3_GRCh38.p7-raw.vcf.gz \
--output vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz
```

```
module load gcc/6.2.0 bedtools/2.27.1

bedtools intersect \
-v \
-a vcf_files/syn3_GRCh38.p7-raw-filt.snpeff.vcf \
-b LCR_GRCh38.p7.bed > vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf
```


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
