# Variant Filtering

## Learning Objectives

- Filter raw variant calls using `FilterMutectCells` to reduce errors
- Remove Low-Complexity Regions from the called variants using `SnpSift` to further reduce errors

## FilterMutectCalls

The output from `Mutect2` is a raw variant calling output and the calls need to be filtered to ensure against errors such as:

- Technical artifacts
- Non-somatic mutations
- Sequencing Errors

`FilterMutectCalls` evaluates the raw variant calls for each of these types of errors using a probabilitic model for errors. It then uses this model to determine the probability of an error and applies this filter across all of the variants. There are also "hard filters" that immediately flag a variant call for filtering. These include:

- Too many alternate alleles
- Low median base quality scores 
- Low median alignment quality scores

> NOTE: While we are not concerned with cross-sample contamination for this dataset, if you were concerned about cross-sample contamination, then you would need to run `CalculateContamination` program within `GATK` to obtain a contamination table which you can input into `FilterMutectCells` with the `--contamination-table` option.



```
gatk FilterMutectCalls \
--reference /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
--variant vcf_files/syn3_GRCh38.p7-raw.vcf.gz \
--output vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz
```

More information on `FilterMutectCalls` can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls) and a more technical guide to the filtering can be found [here](https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf) in Section II.

## Low-Complexity Regions

Low-complexity regions of the genome represent regions that have simple sequence repeats and variant callers are prone to make errors within these regions (see [Li, 2014](https://academic.oup.com/bioinformatics/article/30/20/2843/2422145)). They identified many insertions and deletions (Indels) were erroreously called within these low-complexity regions by various variant callers. They found:

> "low-complexity regions (LCRs), 2% of the human genome, harbor 80–90% of heterozygous INDEL calls and up to 60% of heterozygous SNPs" with false positive rates ranging from "10% to as high as 40%".

As a result of these high error rates in low-complexity regions, they recommend removal of these regions with high error rates until better methods for variant calling in low-complexity regions can become established.

### BED file format

**B**rowser **E**xtensible **D**ata (BED) is a tab-delimited file format that contains information on genomic features. A BED file's first three columns (Chromosome, Starting Position and Ending Position) are required fields. Some BED files have additional columns but these are not required.

<p align="center">
<img src="../img/bed.png" width="600">
</p>

It is important to note that BED files positioning have ***zero-based indexing***. What does this mean? There are two ways to express genomic coordinates:

- **Zero-based** is shown at the top of the image
- **One-based** is shown at the bottom of the image

<p align="center">
<img src="../img/Interbase.png" width="300">
</p>

The benefits to having a **zero-based** system is the ease of calculating distance or length of sequences. We can easily determine the length of the `ATG` sequence using the zero-based coordinates by subtracting the start from the end, whereas for one-based coordinates we would need to add one after the subtraction. Therefore, many file formats used in computation, including the BED file format, use zero-based coordinates.

### Using `SnpSift` to remove LCRs

The BED file (explained below) containing the LCRs for GRCh38 can be obtained from the supplementary files or directly downloaded from this link:

```
curl -o LCR-hs38.bed.gz -L https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
```

We can then unpack the gzipped files using the following command:

```
gunzip -c LCR-hs38.bed.gz > LCR-hs38.bed
```

When we can inspect our BED file we can see that it simply has the required 3 columns denoting the positioning of low-complexitiy regions in GRCh38.

```
less LCR-hs38.bed
```

In order to remove the LCRs from the VCF file, we will be using `SnpSift`, which is part of the [`SnpEff and SnpSift suite`](http://pcingola.github.io/SnpEff/) of tools. We will be later be using SnpEff to annotate our variants and SnpSift to priotize our variants, but for now were are going to focus on using the `intervals` command build into `SnpSift`.

The syntax for running the `filter` command in `SnpSift` on a VCF file and remove all sites that overlap with the BED file is:

```
java -jar $SNPEFF/SnpSift.jar intervals -x -i input_file.vcf[.gz] blacklisted_sites.bed > output_file.vcf
```

- `-x` This option tells `SnpSift` to *exclude* sites found in the BED file. The default behavior of `SnpSift filter` is to only *include* sites found in the BED file.

- `-i input_file.vcf` This is the VCF file that we would like to be filtered. I can either be `.gz` compressed or not. 

- `blacklisted_sites.bed` This represents the BED file you want to use to filter your VCF file with. While in this example we only have one BED file, you can use multiple BED files if you have several filters that you wanted to apply. 

- `> output_file.vcf` Lastly, this is just redirecting the output into a new, filtered VCF to a file.

So, to analyze our data, we are going to use the following command:

```
module load snpEff/4.3g

java -jar $SNPEFF/SnpSift.jar intervals -x -i vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz LCR-hs38.bed > vcf_files/syn3_GRCh38.p7-LCR-filt.vcf
```

We can look at this VCF file and note a few items that have been added to the meta-information lines:

```
less vcf_files/syn3_GRCh38.p7-LCR-filt.vcf
```

Scroll down the VCF file past all of the contigs and you should see a line starting with:

```
##filtering_status=
```

Let's inspect these lines a little:

```
##filtering_status=These calls have been filtered by FilterMutectCalls to label false positives with a list of failed filters and true positives with PASS.
##normal_sample=syn3-normal
##source=FilterMutectCalls
##source=Mutect2
##tumor_sample=syn3-tumor
##SnpSiftVersion="SnpSift 4.3g (build 2016-11-28 08:32), by Pablo Cingolani"
##SnpSiftCmd="SnpSift int -x -i syn3_GRCh38.p7-raw-filt.vcf.gz ../LCR-hs38.bed"
```

The first five lines have been added to our VCF file by GATK. They give information on the programs that have been run on the data, which is listed on the `##source=` lines. These lines also define the column header in the VCF file that corresponds to the normal (`##normal_sample=`) and tumor sample (`##tumor_sample=`). 

`##SnpSiftVersion=` states the version of SnpSift that was used to produce this file.

`##SnpSiftCmd=` provides the command that was used in SnpSift to carry out the LCR filtering.

---

Likely to be removed.

### Using `bedtools` to remove LCRs

The BED file (explained below) containing the LCRs for GRCh38 can be obtained from the supplementary files or directly downloaded from this link:

```
curl -o LCR-hs38.bed.gz -L https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
```

We can then unpack the gzipped files using the following command:

```
gunzip -c LCR-hs38.bed.gz > LCR-hs38.bed
```

When we can inspect our BED file we can see that it simply has the required 3 columns denoting the positioning of low-complexitiy regions in GRCh38.

[`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html) is an useful suite of tools to use when handling BED files. It also has functionality for handling VCF files. We will be using the `intersect` command with the `-v` option. 

```
module load gcc/6.2.0 bedtools/2.27.1

bedtools intersect \
-v \
-a vcf_files/syn3_GRCh38.p7-raw-filt.snpeff.vcf \
-b LCR_GRCh38.p7.bed > vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf
```

`bedtools intersect` calls `bedtools` to run the `intersect` command

`-v` Traditionally, `bedtools intersect` will report the intersection of the file following `-a`and the file following `-b`. However, the `-v` option alters this behavior to find positions in the `-a` file not in `-b` file.

`-a vcf_files/syn3_GRCh38.p7-raw-filt.snpeff.vcf` VCF file that we want filtered

`-b LCR_GRCh38.p7.bed > vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf` BED file containing genomic coordinates for sites in the VCF file to exclude and redirected into an output file.

---

Now, we have successfuly filtered our raw VCF file to only include high-qulaity variant calls.

[Next lesson >>](variant_prioritization.md)

[Back to Schedule](../schedule/README.md)


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
