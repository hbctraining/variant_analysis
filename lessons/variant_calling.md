# Variant Calling

## Learning Objectives

-
-
- Call somatics variants from `bam` files using Mutect2

## Germline versus Somatic Variant Calling

Variant calling can be broken up into two groups, germline and somatic. Germline variant calling refers to the process of calling variants that are ubiquitous across the organism (i.e. almost all cells carry these variants) and these are the types of variants that can be passed through the germline. Studies that evaluate population genetics are often concerned with germline variant calling. Somatic variant calling refers to the process of calling variants that differ between cells within a single organism and these variants are not passed through the germline. Somatic variant calling is often used when studying the progression of various cancers. These two types of variant calling methods have different assumptions regarding in the input data and thus are handled differently. 

For example, germline variant calling for the most part expects at most two alleles in relatively equal frequencies, while a single tumor sample could have various cancer lineages with various allele frequencies. This makes somatic variant calling more difficult than germline variant calling since low frequency variants and sequencing artifacts are difficult to distinguish from sequencing error. Additionally, oftentimes within somatic variant calling, you are also trying to avoid calling the germline variants. The image below should help distinguish between germline and somatic variants:

<p align="center">
<img src="../img/Germline_Somatic_Variants.png" width="600">
</p>

In the above image we can see an example of a germline variant on the left. Approximately half of the reads support each allele in both the tumor and normal sample reads. When compared to the somatic variant on the right, where we observe no variants in the normal sample reads but there is a variant present in the tumor sample reads. 

As you can also presume, high coverage is helpful for two reasons:

1) It helps distinguish sequencing errors and artifacts from true low frequency alleles in the tumor samples
2) It helps distinugish germline variants from somatic variants

## Types of variants that can be called

There are several different types variants that require their own consideration. These include:

- Single Nucleotide Polymorphisms (SNPs)
- Small Insertions/Deletions (Indels)
- Copy Number Variants (CNVs)
- Structural Variants (SVs)

Similarly to the tools in a workshop, variant calling for each of these types of variants requires tools created for it. `GATK` has packages that can address the needs of several of these:

- [`HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) can be used for germline SNPs and Indels
- [`MuTect2`](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2) can be used for somatic SNPs and Indels
- [`GermlineCNVCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/5358874158235-GermlineCNVCaller) can be used for germline Copy Number Variants
- [Software for Structural Variants is currently in beta testing.](https://gatk.broadinstitute.org/hc/en-us/articles/5358824293659--Tool-Documentation-Index#StructuralVariantDiscovery)

This course is going to focus on analyzing somatic SNPs, so we are going to use `MuTect2`.

## `MuTect2`

### Basic workflow

When using `MuTect2` will first re-evaluate the alignments of the normal and tumor samples and create "active regions" that appear to need a local re-assembly. During this process of local re-assembly the tumor sample's reads are interrogated to a higher degree for their quality than the normal samples and a *de Bruijn* graph of the region is created with an assembler. From here, most likely haplotypes are assembled and variants are called from these haplotypes. In order to call variants using Mutect2:

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
```

## Tumor-only Mode

MuTect2 is capable of running without a matched normal sample, otherwise called "Tumor-only mode". However, it's ability to reliably call somatic variants is greatly diminished as it has diificulty distinguishing between high frequency variants and germline variants. 

## Discuss PoNs?

Panel of Normals are is a VCF of normal samples run through "tumor-only mode". If the variant is seen in multiple samples then this variant is included in the Panel of Normals. If these variants are then found in the tumor-sample then the variant is ignored. If panels of normal are used, then they should be gathered using a similiar sequeuncing design as the tumor samples.

## Discuss Common General Population Allele Frequencies

gnomAD


## Additional Resources

[BroadE: GATK - Intro to Somatic Variant Discovery](https://www.youtube.com/watch?v=0q5_e2Nfph4)

[BroadE: GATK - Somatic SNVs and Indels](https://www.youtube.com/watch?v=T5IqadGoxow)

[Next Lesson >>>](variant_annotation.md)

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
