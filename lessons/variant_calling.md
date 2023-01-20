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

When using `MuTect2` will first re-evaluate the alignments of the normal and tumor samples and create "active regions" that appear to need a local re-assembly. During this process of local re-assembly the tumor sample's reads are interrogated to a higher degree for their quality than the normal samples and a *de Bruijn* graph of the region is created with an assembler. From here, most likely haplotypes are assembled and variants are called from these haplotypes. Let's start writing out a new `sbatch` submission script for MuTect2`:

```
cd ~/variant_calling/scripts/
vim mutect2_normal_tumor.sbatch
```

First, let's add our shebang line, script description, `sbatch` directives and GATK module.

```
#!/bin/bash
# This sbatch script is for variant calling with GATK's MuTect2

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o mutect2_variant_calling_normal_tumor_%j.out
#SBATCH -e mutect2_variant_calling_normal_tumor_%j.err

# Load the GATK module
module load gatk/4.1.9.0
```

Next, we need to add our variables:

```
# Assign variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
ALIGNMENT_DIRECTORY=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/
VCF_DIRECTORY=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/
NORMAL_SAMPLE_NAME=syn3-normal
TUMOR_SAMPLE_NAME=syn3-tumor

REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
REFERENCE_SEQUENCE_NAME=`basename $REFERENCE_SEQUENCE _genomic.fa`
NORMAL_BAM_FILE=${ALIGNMENT_DIRECTORY}${NORMAL_SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.coordinate_sorted.bam
TUMOR_BAM_FILE=${ALIGNMENT_DIRECTORY}${TUMOR_SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.coordinate_sorted.bam
VCF_OUTPUT_FILE=${VCF_DIRECTORY}mutect2_${NORMAL_SAMPLE_NAME}_${TUMOR_SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-raw.vcf
```

> NOTE: Sometimes when there are many input variables that re-use many of the same textual elements (i.e. paths, sample names and reference genome names), like we have above, it is sometimes cleaner, less typo-prone and more reproducible to assign those repeated items to variables and then use text manipulation tools and variable subsitution in `bash` to create the rest of the variables. In the above example, the first five lines of variable assignment (`REFERENCE_SEQUENCE` to `TUMOR_SAMPLE_NAME`) are likely lines you might edit from run-to-run, but the final five lines (`REFERENCE_DICTIONARY` to `VCF_OUTPUT_FILE`) will likely stay the same. Standardizing your paths and nomenclature will help you keep track of your files much easier.

Lastly, we need to add the `MuTect2` command:

```
# Run MuTect2
gatk Mutect2 \
--sequence-dictionary $REFERENCE_DICTIONARY \
-R $REFERENCE_SEQUENCE \
-I $NORMAL_BAM_FILE \
--normal-sample $NORMAL_SAMPLE_NAME \
-I $TUMOR_BAM_FILE \
--tumor-sample $TUMOR_SAMPLE_NAME \
--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
-O $VCF_OUTPUT_FILE
```

Let's breakdown this command:

- `gatk Mutect2` Calls the `Mutect2` package from `GATK`

- `--sequence-dictionary $REFERENCE_DICTIONARY` `GATK` requires a sequence directory (`.dict`) file of the reference sequence. We have gone ahead and already created this for you

<details>
  <summary></b>Click here for the commands to create a sequence directory</b></summary>
  We can create the required sequence dictionary in <code>Picard</code>. But first, let's double check we have the <code>Picard</code> module loaded:
  <pre>
  module load picard/2.8.0</pre>
  
  The command to do create the sequence dictionary is:<br>
  <pre>
  java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary \
  REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
  OUTPUT=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.dict</pre>
  
  The components of this command are:
  <ul><li><code>java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary</code> This calls the <code>CreateSequenceDictionary</code> command within <code>Picard</code></li>
  <li><code>REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa</code> This is the reference sequence to create the sequence dictionary from.</li>
  <li><code>OUTPUT=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.dict</code> This is the output sequence dictionary from.</li></ul>
  
  Like indexing, once you have created the sequence dictionary for a reference genome once, you won't need to do it again.
</details>

- `-R $REFERENCE_SEQUENCE` This is the genome reference sequence

- `-I $NORMAL_BAM_FILE` This is the first `bam` file that we are providing GATK and it happens to be the normal sample

- `--normal-sample $NORMAL_SAMPLE_NAME` This is the name of that normal sample and will be used as a column header in the VCF file

- `-I $TUMOR_BAM_FILE` This is the second `bam` file that we are providing GATK and it happens to be the tumor sample

- `--tumor-sample $TUMOR_SAMPLE_NAME` This is the name of that tumor sample and will be used as a column header in the VCF file

> NOTE: It is **VERY IMPORTANT** that the sample names (`--normal-sample $NORMAL_SAMPLE_NAME` and `--tumor-sample $TUMOR_SAMPLE_NAME`) are provided in the same order as the `-I` input BAM files!

- ```--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage``` These are a variety of additional annotations that we are going to add to the output VCF file. These are not required for `MuTect2` to run, but they provide additional details about our variants.

- `-O $VCF_OUTPUT_FILE` This is our output VCF file

> NOTE: In order to run `MuTect2` we also need to have a FASTA index file of our reference sequence in addition to our sequence dictionary. Similarly to the sequence dictionary and `bwa` indicies, we have already created this index for you. However, the dropdown below will walk you through how to do it, should you ever need to do it on your own:
>
><details>
>  <summary><b>Click here for details for creating a FASTA index file in <code>samtools</code></b></summary>
>    <br>FASTA index files for reference sequences are fairly common requirements for a variety of NGS software packages. <code>Picard</code> currently does not feature an ability to create a FASTA index file. However, <code>samtools</code> is a very popular tool that is used for a variety of processes for processing BAM/SAM files, but it also includes functionality for the creation of FASTA index files. First, we will need to load the `gcc` and `samtools` modules:
>  
>  <pre>
>  module load gcc/6.2.0
>  module load samtools/1.15.1</pre>
>  
>  The command for indexing a FASTA file is straightforward and should run pretty quickly:
>  <pre>
>    # YOU DON'T NEED TO RUN THIS
>    samtools faidx \
>    reference_sequence.fa</pre>
>  
>    We can breakdown this code:
>    <ul><li><code>samtools faidx</code> This calls the <code>faidx</code> software from <code>samtools</code></li>
>    <li><code>reference_sequence.fa</code> This is the reference sequence FASTA file that you would like to index</li></ul>
>  
>  Once the indexing is complete, then you shoudl have a index file (<code>reference_sequence.fa.fai</code>) in same directory as your reference sequence <code>reference_sequence.fa</code>.
></details>


The final `sbatch` submission script for `MuTect2` should look like:

```
#!/bin/bash
# This sbatch script is for variant calling with GATK's MuTect2

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o mutect2_variant_calling_normal_tumor_%j.out
#SBATCH -e mutect2_variant_calling_normal_tumor_%j.err

# Load the GATK module
module load gatk/4.1.9.0

# Assign variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
ALIGNMENT_DIRECTORY=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/
VCF_DIRECTORY=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/
NORMAL_SAMPLE_NAME=syn3-normal
TUMOR_SAMPLE_NAME=syn3-tumor

REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
REFERENCE_SEQUENCE_NAME=`basename $REFERENCE_SEQUENCE _genomic.fa`
NORMAL_BAM_FILE=${ALIGNMENT_DIRECTORY}${NORMAL_SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.coordinate_sorted.bam
TUMOR_BAM_FILE=${ALIGNMENT_DIRECTORY}${TUMOR_SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.coordinate_sorted.bam
VCF_OUTPUT_FILE=${VCF_DIRECTORY}mutect2_${NORMAL_SAMPLE_NAME}_${TUMOR_SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-raw.vcf

# Run MuTect2
gatk Mutect2 \
--sequence-dictionary $REFERENCE_DICTIONARY \
-R $REFERENCE_SEQUENCE \
-I $NORMAL_BAM_FILE \
--normal-sample $NORMAL_SAMPLE_NAME \
-I $TUMOR_BAM_FILE \
--tumor-sample $TUMOR_SAMPLE_NAME \
--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
-O $VCF_OUTPUT_FILE
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
