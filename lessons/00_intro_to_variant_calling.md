---
title: "Set up and overview for gene-level differential expression analysis"
author: "Will Gammerdinger, Meeta Mistry"
date: "December 9, 2023"
---

Approximate time: 45 minutes

## Learning Objectives 

- Describe different types of variants
- Differentiate between somatic and germline variants
- Describe the key steps in the variant calling pipeline

Genomic variants are the basis for many diseases and are the raw material for evolutionary processes. As such, inteerest in analyzing genomic variants is of great interest to clinical physicians, evolutionary biologiests and many more. This workshop will aims to guide participants through the process of:

1. Calling variants
2. Annotating variants
3. Prioritizing variants

However, before we can get started, we first must talk about the different types of variants that exist and the challenges that researchers face when analyzing them.

## Germline versus Somatic Variant Calling

Variant calling can be broadly broken up into two groups, germline and somatic. Germline variant calling refers to the process of calling variants that are ubiquitous across the organism (i.e. almost all cells carry these variants) and these are the types of variants that can be passed through the germline. Studies that evaluate population genetics are often concerned with germline variant calling. Somatic variant calling refers to the process of calling variants that differ between cells within a single organism and these variants are not passed through the germline. Somatic variant calling is often used when studying the progression of various cancers. These two types of variant calling methods have different assumptions regarding in the input data and thus are handled differently. 

For example, germline variant calling for the most part expects at most two alleles in relatively equal frequencies, while a single tumor sample could have various cancer lineages with various allele frequencies. This makes somatic variant calling more difficult than germline variant calling because low frequency variants and sequencing artifacts are difficult to distinguish from sequencing errors. Additionally, oftentimes within somatic variant calling, you are also trying to avoid calling the germline variants. The image below should help distinguish between germline and somatic variants:

<p align="center">
<img src="../img/Germline_Somatic_Variants.png" width="600">
</p>

In the above image we can see an example of a germline variant on the left. Approximately half of the reads support each allele in both the tumor and normal sample reads. When compared to the somatic variant on the right, where we observe no variants in the normal sample reads but there is a variant present in the tumor sample reads. 

As you can also presume, high coverage is helpful for two reasons:

**1)** It helps distinguish sequencing errors and artifacts from true low frequency alleles in the tumor samples

**2)** It helps distinugish germline variants from somatic variants

## Types of variants that can be called

There are several different types variants that require their own consideration. These include:

- Single Nucleotide Polymorphisms (SNPs)
- Small Insertions/Deletions (Indels)
- Copy Number Variants (CNVs)
- Structural Variants (SVs)

Similarly to the tools in a workshop, variant calling for each of these types of variants requires tools created for it. [`GATK`](https://gatk.broadinstitute.org/hc/en-us) has packages that can address the needs of several of these:

- [`HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) can be used for germline SNPs and Indels
- [`MuTect2`](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2) can be used for somatic SNPs and Indels
- [`GermlineCNVCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/5358874158235-GermlineCNVCaller) can be used for germline Copy Number Variants
- [Software for Structural Variants is currently in beta testing.](https://gatk.broadinstitute.org/hc/en-us/articles/5358824293659--Tool-Documentation-Index#StructuralVariantDiscovery)

This course is going to focus on analyzing somatic SNPs, so we are going to use `MuTect2`.


