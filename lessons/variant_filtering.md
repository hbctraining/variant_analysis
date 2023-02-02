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

Let's begin by navigating to our scripts directory and creating the `sbatch` submission script that we will be using for filtering our VCF file:

```
cd ~/variant_calling/scripts/
vim variant_filtering_normal_tumor.sbatch
```

The first step is to add our shebang line, description and `sbatch` directives: 

```
#!/bin/bash
# This sbatch script is for variant filtering 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_filtering_normal_tumor_%j.out
#SBATCH -e variant_filtering_normal_tumor_%j.err
```

Next, we need to add the modules that we will be loading:

```
# Load modules
module load gatk/4.1.9.0
```

Next, we will add our variables:

```
# Assign variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
RAW_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf
MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}filt.vcf
```

Next, we can add the `FilterMutectCells` command:

```
# Filter Mutect Calls
gatk FilterMutectCalls \
--reference $REFERENCE_SEQUENCE \
--variant $RAW_VCF_FILE \
--output $MUTECT_FILTERED_VCF
```

Let's breakdown this command:

- `gatk FilterMutectCalls` Calls the `FilterMutectCalls` package within `GATK`

- `--reference $REFERENCE_SEQUENCE` This is our reference genome

- `--variant $RAW_VCF_FILE` This is our raw VCF file from `MuTect2`

- `--output $MUTECT_FILTERED_VCF` This is our filtered output file from `FilterMutectCalls`

More information on `FilterMutectCalls` can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls) and a more technical guide to the filtering can be found [here](https://github.com/broadinstitute/gatk/blob/master/docs/mutect/mutect.pdf) in Section II.

Save and exit `vim`.

## Low-Complexity Regions

Low-complexity regions of the genome represent regions that have simple sequence repeats and variant callers are prone to make errors within these regions (see [Li, 2014](https://academic.oup.com/bioinformatics/article/30/20/2843/2422145)). It has been identified some insertions and deletions (Indels) are erroreously called within these low-complexity regions by various variant callers. They found:

> "low-complexity regions (LCRs), 2% of the human genome, harbor 80–90% of heterozygous INDEL calls and up to 60% of heterozygous SNPs" with false positive rates ranging from "10% to as high as 40%".

As a result of the high error rates in these low-complexity regions, it is recommended to remove of these regions with high error rates until better methods for variant calling in low-complexity regions can become established.

### Using `SnpSift` to remove LCRs

#### Downloading and Unpacking BED with LCRs
The BED file (explained below) containing the LCRs for GRCh38 can be directly downloaded from this link:

```
cd ~/variant_calling/reference/
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

#### Editing variant filtering script to include LCR filtering

In order to remove the LCRs from the VCF file, we will be using `SnpSift`, which is part of the [`SnpEff and SnpSift suite`](http://pcingola.github.io/SnpEff/) of tools. We will be later be using `SnpEff` to annotate our variants and `SnpSift` to priotize our variants, but for now were are going to focus on using the `intervals` command build into `SnpSift`. Let's go back to our scripts directory and edit our variant filtering script.

```
cd ~/variant_calling/scripts/
vim variant_filtering_normal_tumor.sbatch
```

Let's go ahead and make some edits to our variant filtering scripts. First, we need to add `SnpEff` to the modules we will be loading. So the modules loaded will now look like:

```
# Load modules
module load gatk/4.1.9.0
module load snpEff/4.3g
```

Next, we need to add some additional `bash` variables:

```
# Assign variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
RAW_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf
LCR_FILE=/n/groups/hbctraining/variant_calling/reference/LCR-hs38.bed
MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}filt.vcf
LCR_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}LCR-filt.vcf
```

Add the `filter` command from `SnpSift` in order remove all sites that overlap with the BED file:

```
# Filter LCR
java -jar $SNPEFF/SnpSift.jar intervals \
-x \
-i $MUTECT_FILTERED_VCF \
$LCR_FILE > $LCR_FILTERED_VCF
```

- `-x` This option tells `SnpSift` to *exclude* sites found in the BED file. The default behavior of `SnpSift filter` is to only *include* sites found in the BED file.

- `-i $MUTECT_FILTERED_VCF` This is the VCF file that we would like to be filtered. I can either be `.gz` compressed or not. 

- `$LCR_FILE` This represents the BED file you want to use to filter your VCF file with. While in this case we only have one BED file, you can use multiple BED files if you have several filters that you wanted to apply. 

- `> $LCR_FILTERED_VCF` Lastly, this is just redirecting the output into a new, filtered VCF to a file.


[`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html) is an useful suite of tools to use when handling BED files. It also has functionality for handling VCF files. A similar approach for filtering out low-complexity regions can also be done within the tool `bedtools` and is shown in a dropdown box below.

<details>
  <summary>Click here to see how to use <code>bedtools</code> to exclude sites</summary>

First, we will need to load the <code>bedtools</code> module. The <code>bedtools</code> module 

<pre>
module load gcc/9.2.0
module load bedtools/2.30.0</pre>
  
The command for running <code>bedtools</code> to filter out low-complexity regions is:

<b><i>TEST THIS STILL</i></b>
<pre>
bedtools intersect \
-header \
-v \
-a $MUTECT_FILTERED_VCF \
-b $LCR_FILE > $LCR_FILTERED_VCF
</pre>
  
We can breakdown this command:
  
<ul><li><code>-header</code>This will maintain the headers from the VCF-file</li>
  
<li><code>-v</code>Traditionally, <code>bedtools intersect</code> will report the intersection of the file following <code>-a</code> and the file following <code>-b</code>. However, the <code>-v</code> option alters this behavior to find positions in the <code>-a</code> file not in <code>-b</code> file.</li>
  
<li><code>-a $MUTECT_FILTERED_VCF</code> VCF file that we want filtered</li>
  
<li><code>-b $LCR_FILE</code>BED file containing genomic coordinates for sites in the VCF file to exclude</li>
  
<li><code>> $LCR_FILTERED_VCF</code></li> Redirecting the output of this filtering command to a new file.</li></ul>
</details>

Our final `sbatch` script should look like:

```
#!/bin/bash
# This sbatch script is for variant filtering 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_filtering_normal_tumor_%j.out
#SBATCH -e variant_filtering_normal_tumor_%j.err

# Load modules
module load gatk/4.1.9.0
module load snpEff/4.3g

# Assign variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
RAW_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf
LCR_FILE=/n/groups/hbctraining/variant_calling/reference/LCR-hs38.bed
MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}filt.vcf
LCR_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}LCR-filt.vcf

# Filter Mutect Calls
gatk FilterMutectCalls \
--reference $REFERENCE_SEQUENCE \
--variant $RAW_VCF_FILE \
--output $MUTECT_FILTERED_VCF

# Filter LCR
java -jar $SNPEFF/SnpSift.jar intervals \
-x \
-i $MUTECT_FILTERED_VCF \
$LCR_FILE > $LCR_FILTERED_VCF
```

We can now save and exit `vim`, then submit our `sbatch` script.

```
sbatch variant_filtering_normal_tumor.sbatch
```

## Inspecting Filtered VCF File

Once it has completed, which should be quickly. We can look at the output VCF file and note a few items that have been added to the meta-information lines:

```
less /n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_normal_syn3_tumor_GRCh38.p7-LCR-filt.vcf
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
##SnpSiftCmd="SnpSift int -x -i syn3_GRCh38.p7-raw-filt.vcf ../LCR-hs38.bed"
```

The first five lines have been added to our VCF file by GATK. They give information on the programs that have been run on the data, which is listed on the `##source=` lines. These lines also define the column header in the VCF file that corresponds to the normal (`##normal_sample=`) and tumor sample (`##tumor_sample=`). 

`##SnpSiftVersion=` states the version of SnpSift that was used to produce this file.

`##SnpSiftCmd=` provides the command that was used in SnpSift to carry out the LCR filtering.

---

Now, we have successfuly filtered our raw VCF file to only include high-qulaity variant calls.

[Next lesson >>](variant_annotation.md)

[Back to Schedule](../schedule/README.md)


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
