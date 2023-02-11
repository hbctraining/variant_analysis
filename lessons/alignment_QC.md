# Alignment Quality Control

## Learning Objectives

- Verify alignment rates using `Picard`
- Merge `Picard` QC metrics with `FastQC` metrics using `MultiQC`

## Factors Impacting Alignment

One of the most important metrics for your alignment file is the alignment rate. Alignment rates can vary based upon many factors, including:

- **Quality of reference assembly** - A high-quality assembly like GRCh38.p7 will provide a more than adequate reference geneome for alignment. However, if you were studying a organism with a poorly assembled genome, parts of the reference genome could be missing from the assembly. Therefore, high-quality reads might not align because they there is missing reference sequence to align to that corresponds to their sequence.
- **Quality of libraries** - If the library generation was poor and there wasn't enough input DNA, then your sequencing could be filled with low-quality reads
- **Quality of the reads** - If the reads are poor quality, then it can make alignment more uncertain. If your `FASTQC` report shows any anomalous signs, contact your sequencing center for support.
- **Contamination** - If your samples are contaminated, then it can also skew your alignment. For example, if your samples were heavily contaminated with some bacteria, then much of what you will sequence will be bacteria DNA and not your sample DNA. As a result, most of the sequence reads will not align to your target sequence. If you suspect contamination might be the source of a poor alignment, you could consider running [Kraken](https://ccb.jhu.edu/software/kraken/) to evaluate the levels of contamination in your samples.
- **Evolutionary distance between the sampled organism and reference genome** If a reference genome doesn't exist for your species of interest, you are able to align reads to a closely-related organism. However, it does come at the cost of lowering the alignment rate. 
- **Aligner and alignment parameters** Different aligners work differently and are specialized for different types of data. Additionally, many aligners have a variety of parameters that are able to be adjusted. As a result, different aligners or different parameters for the same aligner will give different alignment rates, but they usually should be within the same approximate alignment rate. Generally speaking, the default parameters for most alignment tools are usually fine and they shouldn't need to be manually adjusted/optimized unless there is a specific reason to do so.

When aligning high-quality reads to a high quality reference genome, one should expect to see alignment rates at 90% or better. If alignment rates dipped below 80-85%, then there could be reason for furtheer inspection. 

## Collecting Alignment Statistics

We are going to use `Picard` once again in order to collect our alignment statistics. `Picard` has many packages for collecting different types of data, but the one we will be using is `CollectAlignmentSummaryMetrics`. Let's start creating an `sbatch` script that can utilize this package:

```
cd ~/variant_calling/scripts/
vim picard_CollectAlignmentMetrics_normal.sbatch
```

First, we need to add our shebang line, description and `sbatch` directives to the script:

```
#!/bin/bash
# This sbatch script is for collecting alignment metrics using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o picard_CollectAlignmentMetrics_normal_%j.out
#SBATCH -e picard_CollectAlignmentMetrics_normal_%j.err
```

Next we need to load `Picard`:

```
# Load picard
module load picard/2.8.0
```

Next, let's assign our files to variables:

```
# Assign variables
INPUT_BAM=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa
OUTPUT_METRICS_FILE=/home/${USER}/variant_calling/reports/picard/syn3_normal/syn3_normal_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
```

Next, we can add the `Picard` command to gather the alignment metrics:

```
# Run Picard CollectAlignmentSummaryMetrics
picard CollectAlignmentSummaryMetrics \
INPUT=$INPUT_BAM \
REFERENCE_SEQUENCE=$REFERENCE \
OUTPUT=$OUTPUT_METRICS_FILE 
```

We can breakdown this command into each of it's components:

- `picard CollectAlignmentSummaryMetrics` Calls the `CollectAlignmentSummaryMetrics` from within `Picard`

- `INPUT=$INPUT_BAM` This is the output from our previous `Picard` alignment processing steps.

- `REFERENCE_SEQUENCE=$REFERENCE` This isn't a required parameter, but `picard` can do a subset of mismatch-related metrics if this is provided.

- `OUTPUT=$OUTPUT_METRICS_FILE` This is the file to write the output metrics to.

The `sbatch` submission script for collecting the alignment metrics should look like:

```
#!/bin/bash
# This sbatch script is for collecting alignment metrics using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o picard_CollectAlignmentMetrics_normal_%j.out
#SBATCH -e picard_CollectAlignmentMetrics_normal_%j.err

# Load picard
module load picard/2.8.0

# Assign variables
INPUT_BAM=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa
OUTPUT_METRICS_FILE=/home/${USER}/variant_calling/reports/picard/syn3_normal/syn3_normal_GRCh38.p7.CollectAlignmentSummaryMetrics.txt

# Run Picard CollectAlignmentSummaryMetrics
picard CollectAlignmentSummaryMetrics \
INPUT=$INPUT_BAM \
REFERENCE_SEQUENCE=$REFERENCE \
OUTPUT=$OUTPUT_METRICS_FILE 
```

Create the tumor version of this submission script using `sed`:

```
sed 's/normal/tumor/g' picard_CollectAlignmentMetrics_normal.sbatch > picard_CollectAlignmentMetrics_tumor.sbatch
```

Now that we have created these files to submit, let's check the status of our previous `Picard` alignment processing steps:

```
squeue -u $USER
```

**If your `Picard` alignment processing steps are not completed yet**, wait until they have finished before submitting these jobs to collect alignment metrics.

**If your `Picard` alignment processing steps are completed**, then submit these jobs to collect alignment metrics:

```
sbatch picard_CollectAlignmentMetrics_normal.sbatch
sbatch picard_CollectAlignmentMetrics_tumor.sbatch
```

> **NOTE:** The syntax that `Picard` uses is quite particular and you may note in your error file the warning:
> ```
>********* NOTE: Picard's command line syntax is changing.
>*********
>********* For more information, please see:
>********* https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)
>*********
>********* The command line looks like this in the new syntax:
>*********
>*********    CollectAlignmentSummaryMetrics -INPUT /n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam -OUTPUT /home/${USER}/variant_calling/reports/picard/syn3_normal/syn3_normal.CollectAlignmentSummaryMetrics.txt -REFERENCE_SEQUENCE /n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa
>*********
> ```
> The version of `Picard` on the cluster is a bit old, but most of the functionality is still the same. At the time of version 2.8.0, not all packages were converted over to using this syntax. For example, in O2's version of `SortSam`, it doesn't accept this syntax. In order for all of our `Picard` syntax to be consistent in this workshop, we are showing the older syntax. Importantly, the more recent versions `Picard` *likely* all accept the new syntax, but [it should still accept the older syntax](https://github.com/broadinstitute/picard/wiki/Barclay-Transition-Notes). However, the `Picard` manuals themselves are not always the most consistent when differentiating between these two syntaxes and sometimes show both, so we tend to recommend using the old syntax as we find it works more consistently.


## Collecting Coverage Metrics

Coverage is the average level of alignment for any random locus in the genome.  `Picard` also has a package called `CollectWgsMetrics` which is also very nice for collecting data about coverage for our alignments. However, since our data set is whole exome sequencing rather than whole genome sequencing and thus only compromises about 1-2% of the human genome, average coverage across the whole genome is not a very useful metric. However, if one did have whole genome data, then running `CollectWgsMetrics` would be useful and even could be incorporated easily into the downstream <code>MultiQC</code> HTML report. In the dropdown box below be provide the code that you could use to collect this information.

<details>
<summary>Click here to find out more on collecting coverage metrics for WGS datasets in <code>Picard</code></summary>
<br>The tool in <code>Picard</code> used for collecting coverage metrics for WGS datasets is called <code>CollectWgsMetrics</code>.<br><br>
  <pre>
  # Assign paths to bash variables
  $COORDINATE_SORTED_BAM_FILE=/path/to/sample.coordinate_sorted.bam
  $OUTPUT=/home/$USER/variant_calling/reports/picard/sample.CollectWgsMetrics.txt
  $REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa<br>
  # Run Picard CollectWgsMetrics \
  picard CollectWgsMetrics \
  INPUT=$COORDINATE_SORTED_BAM_FILE \
  OUTPUT=$METRICS_OUTPUT_FILE \
  REFERENCE_SEQUENCE=$REFERENCE
  </pre>
        
  <ul><li><code>picard CollectWgsMetrics</code> This calls the <code>CollectWgsMetrics</code> package within <code>Picard</code></li>
  <li><code>INPUT=$COORDINATE_SORTED_BAM_FILE</code> Assign the input as the coordinate sorted BAM file</li>
  <li><code>OUTPUT=$METRICS_OUTPUT_FILE</code> Assign the report output file </li>
  <li><code>REFERENCE_SEQUENCE=$REFERENCE</code> This is the path to the reference genome that was used for the alignment.</li></ul>
<hr />
</details>

## Inspecting procressed `BAM` files

We discussed BAM/SAM file formatting in the [file format lesson](file_formats.md), but didn't discuss the header section and we could now inspect our processed BAM files to see what that looks like. In order to do this, we are going to use `Picard`, but once again this can be done in `samtools` as well (and more oftentimes is done in `samtools`) and that will be shown in a dropdown at the end of this section.

First, we will need to make sure that the `Picard` module is loaded:

```
module load picard/2.8.0
```

Next, we will run the `ViewSam` package in `Picard`:

```
picard ViewSam INPUT=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam | less
```

This will open up the fully processed normal sample's BAM file in a human-readable format.

We can see the top line is the header and it tells you information about the alignment file:

```
@HD     VN:1.5  SO:coordinate
```

`VN` is telling you the version number for the BAM/SAM formatting and `SO` is telling you the sort order. If we were to open up, our query-sorted BAM file, we would see that this says "queryname" instead of "coordinate".

Then, you have all of your sequence lines (@SQ) outlining all of the reference sequence chromosomes and contigs.

```
@SQ     SN:1    LN:248956422
@SQ     SN:HSCHR1_CTG1_UNLOCALIZED      LN:175055
.
.
.
@SQ     SN:MT   LN:16569
```

After that you can get to the Read Group (@RG) and it has the read group information that we provided `bwa`:

```
@RG     ID:syn3_normal  PL:illumina     PU:syn3_normal  SM:syn3_normal
```

Next, we get some very useful lines describing how this alignment file has been processed and the software used to do it in the program lines (@PG):

```
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem -M -t 8 -R @RG\tID:syn3_normal\tPL:illumina\tPU:syn3_normal\tSM:syn3_normal /n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa /home/${USER}/variant_calling/raw_data/syn3_normal_1.fq.gz /home/${USER}/variant_calling/raw_data/syn3_normal_2.fq.gz -o /n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.sam
@PG     ID:MarkDuplicates       VN:2.8.0-SNAPSHOT       CL:picard.sam.markduplicates.MarkDuplicates INPUT=[/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.query_sorted.bam] OUTPUT=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.remove_duplicates.bam METRICS_FILE=/home/${USER}/variant_calling/reports/picard/syn3_normal/syn3_normal.remove_duplicates_metrics.txt REMOVE_DUPLICATES=true    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag ASSUME_SORTED=false DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json        PN:MarkDuplicates
```

These program lines are added in the order that they are implemented. The @PG lines do not include sorting commands when using <code>Picard</code>, but does include the sorting commands when using <code>samtools</code>. It's not critical and all of the downstream analyses will work either way, but it is a difference to note. We can see that `bwa` was run first because it is in the `ID` and `PN` (Program Name) fields. We are provided the version number (`VN` field) and the exact command-line command that was used (`CL` field). This can all be really useful if you have to try to reproduce someone's work from an alignment file, but the scripts that were originally used to produce it aren't availible.

Lastly, we can see the aligned reads after all of this metadata concerning the file.

<details>
<summary><b>Click here to see how to inspect BAM files in <code>samtools</code></b></summary>
If we wanted to use <code>samtools</code> to view the BAM/SAM file, we would first need to make sure the <code>samtools</code> module is loaded (Note: that <code>samtools</code> does require <code>gcc</code> to be loaded as well:

<pre>
module load gcc/6.2.0
module load samtools/1.15.1
</pre>

Then we can view our normal sample alignment file with:
  
<pre>
samtools view /n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam  | less
</pre>
  
We can break down this command:
  
<ul><li><code>samtools view</code> Calls the <code>view</code> package from within <code>samtools</code></li>

<li><code>/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam</code> This is the alignment file we want to view</li>

<li><code>| less</code> Finally, we are going to pipe this output into a <code>less</code> command so that it doesn't fill up our screen</li></ul>

However, <code>samtools</code> is a bit strange in that it skips over the header lines when you do this, so in order to see the header lines we need to use the <code>-H</code> option.

<pre>
samtools view -H /n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam  | less
</pre>

The <code>-H</code> option modifies the output to only print the header lines. You won't get the read lines as well.

The information here is the same as with <code>Picard</code>, so we won't rehash it. There are only differences that <code>samtools</code> adds the <code>view</code> command into the @PG lines.
<hr />
</details>

## Options for Inspecting `Picard` Alignment Metrics

Once the job has finished we would inspect the output files. This could be done in one of a few ways:

1) View each metrics file in a `less` buffer
  
    **Pros:**
      - Simple to do
  
    **Cons:**
      - Hard to compare across samples
      - Tedious to parse the columns
    
2) Download each metrics from the O2 cluster and import them into Excel/Excel-like program that puts tab-delmited files into a grid
  
    **Pros:**
      - Easier to interpret the data than the `less` buffer approach
  
    **Cons:**
      - Hard to compare across samples
      - Have to download the metrics files from the O2 cluster
    
3) Collate metrics files using [`MultiQC`](https://multiqc.info) and download the `MultiQC` HTML report from the O2 cluster
  
    **Pros:**
      - Alignment metrics are easy to compare across samples
      - Easy to interpret results
      - We can combine the `picard CollectAlignmentSummaryMetrics` output with our `FASTQC` reports
  
    **Cons:**
      - Have to download a file from the O2 cluster
      - Have to run samples through an extra `MultiQC` step

None of the above methods are wrong, but some are more elegant than others. One might use **Method 1)** if they only had a handful of samples (~<5) to analyze and only wanted a single statistic, like alignment rate, from each. Then, it might be fastest just to open them up in a `less` buffer. However, if one has lots of samples (>5) then the advantages of `MultiQC` collating the results starts to become really helpful and one might choose **Method 3)**. Unfortunately, `MultiQC` doesn't display *ALL* of the data contained in the metrics file, so one may be inclined to do **Method 2** and downloard the directory full of metrics files in order to view the metrics not included in the `MultiQC` report in a program like Microsoft Excel.

However, for this workshop, we are going to collate our results in `MultiQC` and download the HTML report to our local computers.

## Inspecting `Picard` Alignment Metrics

One nice feature of `MultiQC` is that it accepts many different file formats. It figures out which format was submitted and tailors the report to that type of analysis. Collating our `MultiQC` results would be relatively quick to just run from the command-line, but it's best practice to write our steps to scripts so that we always have a record of what we did and how we created our reports. We will start by writing a `sbatch` script in `vim` for submission:

```
vim multiqc_alignment_metrics.sbatch
```

First, we will add our sheband line, description and `sbatch` directives

```
#!/bin/bash
# This sbatch script is for collating alignment metrics from Picard using MultiQC 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:10:00
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -o multiqc_alignment_metrics_%j.out
#SBATCH -e multiqc_alignment_metrics_%j.err
```

Next, we will load out modules:

```
# Load modules
module load gcc/9.2.0
module load multiqc/1.12
```
> NOTE: `MultiQC` version 1.12 requires `gcc/9.2.0` on the O2 cluster.

Next, we will assign our variables:

```
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/
NORMAL_SAMPLE_NAME=syn3_normal
TUMOR_SAMPLE_NAME=syn3_tumor

NORMAL_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
TUMOR_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
NORMAL_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_1_fastqc.zip
NORMAL_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_2_fastqc.zip
TUMOR_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_1_fastqc.zip
TUMOR_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_2_fastqc.zip
OUTPUT_DIRECTORY=${REPORTS_DIRECTORY}/multiqc/
```

Next, we need to add the output directory:

```
mkdir -p $OUTPUT_DIRECTORY
```

Then, we will add the command to run `MultiQC`:

```
multiqc \
$NORMAL_PICARD_METRICS \
$TUMOR_PICARD_METRICS \
$NORMAL_FASTQC_1 \
$NORMAL_FASTQC_2 \
$TUMOR_FASTQC_1 \
$TUMOR_FASTQC_2 \
-o $OUTPUT_DIRECTORY
```

So our final `sbatch` script should look like:

```
#!/bin/bash
# This sbatch script is for collating alignment metrics from Picard using MultiQC 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:10:00
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -o multiqc_alignment_metrics_%j.out
#SBATCH -e multiqc_alignment_metrics_%j.err

# Load modules
module load gcc/9.2.0
module load multiqc/1.12

REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/
NORMAL_SAMPLE_NAME=syn3_normal
TUMOR_SAMPLE_NAME=syn3_tumor

NORMAL_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
TUMOR_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
NORMAL_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_1_fastqc.zip
NORMAL_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_2_fastqc.zip
TUMOR_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_1_fastqc.zip
TUMOR_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_2_fastqc.zip
OUTPUT_DIRECTORY=${REPORTS_DIRECTORY}/multiqc/

mkdir -p $OUTPUT_DIRECTORY

multiqc \
$NORMAL_PICARD_METRICS \
$TUMOR_PICARD_METRICS \
$NORMAL_FASTQC_1 \
$NORMAL_FASTQC_2 \
$TUMOR_FASTQC_1 \
$TUMOR_FASTQC_2 \
-o $OUTPUT_DIRECTORY
```

Like the previous step, we will need to check to ensure that the previous `Picard` step for collecting metrics for each sample is down before we can submit this script. To do this, we will check out `squeue`:

```
squeue -u $USER
```

**If your `Picard` collect alignment metric steps are not completed yet**, wait until they have finished before submitting these jobs to `MultiQC`.

**If your `Picard` collect alignment metric steps are completed**, then submit this `MultiQC` job to collate the alignment metrics:

```
sbatch multiqc_alignment_metrics.sbatch
```

This job should finish fairly quickly and then we can proceed to downloading it with `FileZilla` in the next lesson.

[Next Lesson >>>](evaluate_QC.md)

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
