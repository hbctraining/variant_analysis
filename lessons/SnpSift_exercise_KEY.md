---
title: "SnpSift Exercises"
author: "Will Gammerdinger, Meeta Mistry"
date: "May 27, 2024"
---

In this exercise we will use the VCF file that we have generated  in order to practice parsing our apart our VCF. The VCF file is availible for download The VCF file is also availible for download [here](https://hbctraining.github.io/variant_analysis/data/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf).

# Exercises

**1)** Extract all of the stop codon losses

```
java -jar $SNPEFF/SnpSift.jar filter   "( ANN[*].EFFECT has 'stop_lost' ) " mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf 
```

**2)** Extract all of the SNPs on Chromosome 2 between 8Mb and 10.2Mb.

```
java -jar $SNPEFF/SnpSift.jar filter   "( CHROM = '2' ) & ( POS > '8000000' ) & ( POS < '10200000' ) " mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf
```

**3)** Extract all of the mutations that alter a splice site on Chromosome 11.

```
java -jar $SNPEFF/SnpSift.jar filter   "( CHROM = 'X' ) & ( ( ANN[*].EFFECT has 'splice_donor_variant' ) | ( ANN[*].EFFECT has 'splice_acceptor_variant' ) )" mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf
```

**4)** Extract all of the high impact mutations on Chromosome 12 or moderate impact mutations on Chromosome 3.

```
java -jar $SNPEFF/SnpSift.jar filter   "( ( CHROM = '12')  &  ( ANN[*].IMPACT has 'HIGH' ) ) | ( ( CHROM = '13') & ( ANN[*].IMPACT has 'MODERATE' ) )" mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf
```

[Next Lesson >>](11_IGV.md)

[Back to Schedule](../schedule/README.md)


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
