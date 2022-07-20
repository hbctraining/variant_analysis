# Variant Priorization

## Learning Objectives

-
-
-

```
module load snpEff/4.3g

cat vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf | \
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
awk '{print "sample_1","\t",$0}' > vcf_files/syn3_GRCh38.p7.missense_variants.txt
```

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
