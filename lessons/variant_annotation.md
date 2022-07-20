# Variant Annotation

## Learning Objectives

```
module load snpEff/4.3g

java -jar $SNPEFF/snpEff.jar  eff \
-dataDir /n/groups/shared_databases/snpEff.data/ \
-cancer \
-noLog \
-csvStats vcf_files/syn3_hg19-effects-stats.csv \
-s vcf_files/syn3_hg19-effects-stats.html \
GRCh38.p7 \
vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz > vcf_files/syn3_GRCh38.p7-raw-filt.snpeff.vcf
```

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
