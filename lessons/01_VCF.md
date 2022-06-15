---
title: "Variant Call Format"
author: "Radhika Khetani, Meeta Mistry, Jihe Liu, Will Gammerdinger"
date: "May 20, 2022"
---

## Learning Objectives
- Identify different fields within a VCF file
- Discern a data line entry in a VCF file

## Variant Call Format

The [Variant Call Format (VCF)](https://samtools.github.io/hts-specs/VCFv4.2.pdf) is a standardized, text-file format for describing variants identifed from a sequencing experiment. This allows for downstream processes to be streamlined and also allows for researchers to easily collaborate and manipulate a shared set of variant calls. A VCF file is composed of three main parts:
- Meta-information Lines
- Header Line
- Data Lines

### Meta-information Lines

The first lines in a VCF file are called the ***meta-information lines***. These lines contain information regarding how the file was made and what the file contains. A sample of the meta-information lines can be found below:

```
##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta 
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data"> 
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth"> 
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency"> 
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele"> 
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129"> 
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership"> 
##FILTER=<ID=q10,Description="Quality below 10"> 
##FILTER=<ID=s50,Description="Less than 50% of samples have data"> 
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> 
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality"> 
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth"> 
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
```

Meta-information lines always begin with `##`. Let's analyze a few of the meta-information lines to understand what they are describing. 

The first line `##fileformat=VCFv4.2` is telling us the version of the VCF file. 

The third line `##source=myImputationProgramV3.1` is stating the software that was used to create the VCF file. 

The fourth line `##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta` is stating the full path to the reference genome that was used in the creation of this VCF file. 

The lines starting with `##INFO`, `##FILTER` and `##FORMAT` will define abbreviations for us that we will see in the INFO, FILTER and FORMAT fields of the data lines, respectively.

### Header Line

The header line is a single line between the meta-information lines and the data lines that provies a brief desciption of each field in the following data lines. An example header line could look like:

```
#CHROM  POS ID  REF ALT QUAL  FILTER  INFO  FORMAT  NA00001 NA00002 NA00003
```

The header line always starts with just a single `#` followed by eight mandatory fields:

* **CHROM** - Chromosome where the variant was found
* **POS** - A 1-based index for the position on the chromosome where the variant was found. For multibase variants, this corresponds to the first base's position.
* **ID** - If a SNP has an identifier (i.e. such as an rs number(s) from dbSNP), then it is put here. Otherwise, it will be a `.`.
* **REF** - The reference base(s) for the given position
* **ALT** - The variant base(s) for the given position. In the case of multiple variants present, they will be comma separated.
* **QUAL** - A PHRED-scaled quality score for the variant
* **FILTER** - A status of 'PASS' is given for any variant passing all filters. If a variant fails, then a semi-colon separated list will enumerate the filter(s) that the variant failed.
* **INFO** -  Additional information about the variant. Common catergories can be found in the table below:

| Abbreviation | Data Type |
|--------------|-----------|
| AF | Allele Frequency for the ALT allele |
| DP | Combined Depth across all samples |
| NS | Number of samples with data |

* **FORMAT** - A colon-separated list of abbreviated catergories corresponding to the colon-separated genotype fields. For example, this colon-separated list generally begins with 'GT', which stands for genotype. Thus, the first element in the subsequent genotype field(s) for each sample will be the inferred genotype. If the second element in the colon-sperated list is `GQ`, then the second element in each of the colon-separated genotype fieds will be a genotype quality score. Common catergories include:

| Abbreviation | Data Type |
|--------------|-----------|
| GT* | Genotype |
| DP | Depth of the sample |
| GQ | Genotype Quality |
| HQ | Commma-separated list of Haplotype Qualities |

\* The genotype field will have two integers for a diploid sample separated by either a `/` or a `|`. The integers correspond to the alleles with 0 being the reference allele, 1 being the first allele listed in the ALT field, 2 being the second allele listed in the ALT field, etc. A `/` denotes an unphased genotype, while `|` denotes a phased genotype.
* **Genotype Fields** - A colon-separated list of genotype information about a given sample that corresponds to the categories outlined in the FORMAT field. The fields will be equal to the number of samples investigated. In the case of the sample VCF, these are the fields titled NA00001, NA00002 and NA00003.

The abbreviations for the INFO and FORMAT fields given in the aforementioned tables is not exhaustive. Generally, these abbreviations with be outlined in the meta-information lines at the top of the file or can be found on pages 5 and 6 of the [VCF manual](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

### Data Lines

The data lines are where the variant calls will be found with each field corresponding to its column in the header line. 

```
20  14370 rs6054257 G A 29  PASS  NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
20  17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50  0|1:3:5:65,3  0/0:41:3
20  1110696 rs6040355 A G,T 67  PASS  NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27  2|1:2:0:18,2  2/2:35:4
20  1230237  .  T . 47  PASS  NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60  0|0:48:4:51,51  0/0:61:2
20  1234567 microsat1 GTC G,GTCT  50  PASS  NS=3;DP=9;AA=G  GT:GQ:DP  0/1:35:4  0/2:17:2  1/1:40:3   
```

Let's analyze the first line more closely:

```
20  14370 rs6054257 G A 29  PASS  NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51  1|0:48:8:51,51  1/1:43:5:.,.
```

This variant is on Chromosome 20 at position 14370. It is annotated as rs6054257 in dbSNP. The reference allele in the position was a G and a variant allele, A, was found. This variant passed all filters. There were three samples investigated, with a total depth across all three of 14 reads. The allele frequency of the ALT allele in these samples is 0.5. This SNP is included in dbSNP and HapMap2. For each sample we will have information on the genotype, genotype quality, depth of the sample and haplotype quality. 

The first sample has phased information for an individual that is homozygous for the reference allele and a genotype quality of 48. One read is supporting this with haplotype quality scores of 51 supporting it. 

The second sample has phased information for a heterzygous individual with a genotype score of 48. This sample is supported by eight reads and the haplotype quality scores for both haplotypes is 51. 

The last sample has unphased information for an individual who is homozygous for the ALT allele with a genotype quality score of 43. There are five reads supporting this genotype call and since this sample is unphased data, there are no haplotype quality scores as indicated by the `.,.`

## Exercises

Let's analzye our VCF file in `[insert path here]`.

1. Using `grep`, extract only the meta-information lines from the VCF file. 

grep '^##' smaple.vcf

2. What reference genome was used in the creation of this VCF file?

3. For the sample at position A on chromosome B, what is the read depth at this position for sample C?

4. What is the observed allele frequency for the allele at position D on chromosome E?

***
*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*


