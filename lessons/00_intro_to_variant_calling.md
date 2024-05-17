---
title: "Introduction into Variant Analysis"
author: "Will Gammerdinger, Meeta Mistry"
date: "December 9, 2023"
---

Approximate time: 45 minutes

## Learning Objectives 

- Describe different types of variants
- Differentiate between somatic and germline variants
- Describe the key steps in the variant calling pipeline

## What is genomic variation?

Genomic variants are the basis for many diseases and are the raw material for evolutionary processes. Analyzing genomic variants can help inform clinicians and also help discover novel variants responsible for disease. Evolutionary biologists interpret genomic variants to quantify the evolutionary relationships between species. Conservation biologists use genetic variants to measure diversity amongst endangered populations. As such, there is broad interest in analyzing genomic variants to clinicians and biologists alike. 

## Germline versus Somatic Variant Calling

Variant calling can be broadly **broken up into two groups, germline and somatic**. 

**Germline variant calling** refers to the process of calling variants that are ubiquitous across the organism (i.e. almost all cells carry these variants) and these are the types of variants that can be passed through the germline. Studies that evaluate population genetics are often concerned with germline variant calling. 

**Somatic variant calling** refers to the process of calling variants that differ between cells within a single organism and these variants are not passed through the germline. Somatic variant calling is often used when studying the progression of various cancers. 

<p align="center">
<img src="../img/Germline_Somatic_Variants.png" width="600">
</p>

_Example of a germline variant on the left, and somatic variant on the right._ 

These two types of variant calling methods have different assumptions regarding in the input data and thus are handled differently. For example, germline variant calling for the most part expects at most two alleles in relatively equal frequencies, while a single tumor sample could have various cancer lineages with various allele frequencies. This makes somatic variant calling more difficult than germline variant calling because low frequency variants and sequencing artifacts are difficult to distinguish from sequencing errors. Additionally, oftentimes within somatic variant calling, you are also trying to avoid calling the germline variants. 


## Types of variants 

There are several different types of variants that require their own consideration. These include:

- **Single Nucleotide Polymorphisms (SNPs)** - These are positions in the genome where a single base has been mutated. For example, perhaps most individuals in a population have a Thymine in a given position, while an individual of interest has a Adenine in this position.
- **Small Insertions/Deletions (Indels)** - Small indels are loci where a few bases have been added or removed relative to the larger population. For example, if an individual has an extra `GA` at a location relative to the rest of the population, then this would be considered a small insertion.
- **Structural Variants (SVs)** - This class refers to a broad collection of variants, including inversions, translocations and large insertions or deletions. 
- **Copy Number Variants (CNVs)** - These types of variants often occur in repetitive regions of the genome and involve having more of few copies of a given sequence. For example, the *AMY1* gene which encodes for an enzyme that is important in breaking down starches has been shown to have variable numbers of copies across human populations. Further work has shown that the number of copy numbers correlates with with the levels of starch in various cultures (Perry et al., 2007). Depending on the size, copy number variants are sometimes considered a subcategory of structural variants.

Similarly to the tools in a workshop, variant calling for each of these types of variants requires tools created for it. [`GATK`](https://gatk.broadinstitute.org/hc/en-us) is a popular tool for variant calling that was developed and is maintained by the Broad Institute. `GATK` has packages that can address the needs of several types of variant calling:

- [`HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) can be used for germline SNPs and Indels
- [`MuTect2`](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2) can be used for somatic SNPs and Indels
- [`GermlineCNVCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/5358874158235-GermlineCNVCaller) can be used for germline Copy Number Variants
- [Software for Structural Variants is currently in beta testing.](https://gatk.broadinstitute.org/hc/en-us/articles/5358824293659--Tool-Documentation-Index#StructuralVariantDiscovery)

**This course is going to focus on analyzing somatic SNPs, so we are going to use `MuTect2`.**

## The Importance of Coverage

One of the **most important considerations of experimental design** when carrying out a study to identify variants is to sequence your samples to an adequate level of coverage. High coverage is helpful for two reasons:

**1)** It helps distinguish sequencing errors and artifacts from true low frequency alleles in the tumor samples.

**2)** It helps distinugish germline variants from somatic variants.

Coverage simply means for a given position, what is the average number of sequencing reads that span (or "cover") that position and it is abbreviated as the integer value followed by "X". For example, if the average position in the genome was covered by 22 reads, this sample would be considered to have 22X. Generally speaking, it is often encouraged for researchers to **reach a minimum coverage of 30X for variant calling**. However, higher coverage levels can be useful for detecting rare variants, particularly in somatic variant calling.

Given the costs associated with whole genome sequencing (WGS) of each individual to 30X or greater coverage, some researchers opt to simply carry our whole exome sequencing (WES) rather than whole genome sequencing. The Human exome is about ~1% of the human genome's size and many researchers are interested in focusing on transcripts to begin with. Therefore, it can greatly reduce the sequencing costs a researcher might incur if they are willing to forgo the regions not captured in the exome. It should be noted that **due to the uneveness of WES, a greater depth is often encouraged, 70-100X**.

***

**Exercise**

<p align="center">
<img src="../img/Difficulty_of_assignment.png" width="400">
</p>

Use the figure above to try to make inferences answer the following questions:

**1.** If we assume there are no sequencing errors, are you more inclined to speculate that Locus 1 is a germline or somatic variant? Why?

**2.** Given the existence of sequencing errors, how confident are you that Locus 1 represents a heterozygous locus in the germline?

**3.** Given the existence of sequencing errors, how confident are you that Locus 1 represents a polymorphic locus in a somatic tissue?

**4.** How confident are you that Locus 2 is homozygous?

**5.** What additional information might you want in order to better assess these loci?

***

## Variant Calling Workflow

Analyzing variants takes the form of three main steps: 

1. **Data preparation**
2. **Variant calling** 
3. **Variant evaluation** 

The detailed workflow for variant analysis that many researchers follow looks like:

<p align="center">
<img src="../img/Full_variant calling_pipeline.png" width="400">
</p>

**Sequence Samples** - This is the work of your sequencing facility and they should provide you with FASTQ files. 

**Read Quality Control** - It is important to assess the quality of your reads before you further analysis them. In this step we are looking to ensure that there aren't major errors that occurred in the sequencing and that the quality of the reads is high. 

**Align Reads** - We need to match our raw FASTQ reads to the place in the genome where they most likely originiated from. 

**Processing Alignment Files** - While our reads have now been asssigned to a place in the genome, the files holding this information are not in the format that we need in order to call variants. We will be reformatting the alignment files in order to be useable for variant calling.

**Alignment Quality Control** - This is a second opportunity to evaluate the quality of our data. Here we are ensuring that the data aligns well with our genome.

**Call Variants** - Now we are able to call the raw variants from our data. 

**Filter Variants** - Our raw variants may have some artifacts in them and need to be processed to remove them.

**Annotate Variants** - We would be interested to know where our variants fall within the context of current gene models. Here we will be using existing gene models to predict the impact of a variant. Here we will be asking questions such as, is this a variant that creates a premature stop codon or is it in a intergenic region?

**Prioritize Variants** - Assess the annotated variants based upon their predicted functional impacts and narrow the search for important genes that may be causing a disease or trait.

[Next Lesson >>]()

[Back to Schedule](../schedule/README.md)
  
***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*



