# Automation of Variant Calling

## Learning objectives
- Construct a flexible pipeline for automating variant calling
- Integrate the `--dependency` option for `sbatch` into workflows

## Automating the Variant Calling Workflow

Now that we have completed most of the variant calling workflow, we would like to explore ways to automate it to make future runs identical and also minimize downtime in the workflow. Automating this workflow with consist of three parts:

**1)** Editing our `sbatch` submission scripts to accept positional parameters for variables

**2)** Creating a wrapper script to execute each script with the appropriate parameters

**3)** Using the `--dependency` option in `sbatch` to construct a dependent pipeline for scripts that have diffent computational requirements

The workflow and methodolgies we used up to this point are completely fine and could work for you. However, you may have many scripts that you need to create if you have many samples then that can be cumbersome, so we are going to introduce this way of automating the process.

### Editing our sbatch Submission Scripts

Currently, within our `sbatch` submission scripts, we have `bash` variables set to hold various files. Instead of explicitly defining those variables now as a given file, we are going to use positional parameters to allow ourselves to plug-in various files. First let's navigate to our scripts directory:

```
cd ~/variant_calling/scripts/
```

#### FastQC

Because many of some variables already weren't explictly defined (`$RIGHT_READS` and `$SAM_FILE`), but were rather dependent on our explicitly defined variables (`$REFERENCE_SEQUENCE` and `LEFT_READS`) we won't need to edit our script as much. **We are able to do this because `sbatch` accepts positional parameters.** So let's create a new script and make this edit.

```
cp fastqc_normal.sbatch fastqc_automated.sbatch
vim fastqc_automated.sbatch
```
In order for us to specify "tumor" or "normal" samples in the standard error and standard output files we will need to pass `sbatch` these arguments directly when submit the jobs in the wrapper. However, SLURM sometimes has issues when some arguments come from the command line and others come from the script. As a result, we are going remove all `sbatch` directives from scripts and rather provide them all in the command-line as part of the wrapper. Therefore, we need to remove the following lines of `sbatch` directives:

```
#REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -o fastqc_normal_%j.out
#SBATCH -e fastqc_normal_%j.err
```

Next, we need to change our variables to accept positional parameters. So we need to change:

```
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/reports/fastqc/syn3_normal/
THREADS=4
```

To:

```
LEFT_READS=$1
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=$2
THREADS=$3
```

Our `FastQC` submission script should now look like:

```
#!/bin/bash
# This sbatch script is for running FastQC to evaluate read qualities

# Load module
module load fastqc/0.11.9

# Assign variables
LEFT_READS=$1
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=$2
THREADS=$3

# Create directory to hold output
mkdir -p $OUTPUT_DIRECTORY

# Run FastQC
fastqc \
$LEFT_READS \
$RIGHT_READS \
--outdir $OUTPUT_DIRECTORY \
--threads $THREADS
```

#### bwa Alignment

Now, we will need to make similar edits to a new `bwa` script:

```
cp bwa_alignment_normal.sbatch bwa_alignment_automated.sbatch
vim bwa_alignment_automated.sbatch
```

Next, we will need to remove all of the `SBATCH` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_normal_%j.out
#SBATCH -e bwa_alignment_normal_%j.err
```

Next, we will also need to change the variables from:

```
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAMPLE=`basename $LEFT_READS _1.fq.gz
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.sam
```

To:

```
REFERENCE_SEQUENCE=$1
LEFT_READS=$2
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAMPLE=$3
SAM_FILE=$4
```

Our final `bwa` submission script should look like:

```
#!/bin/bash
# This script is for aligning sequencing reads against a reference genome using bwa

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17

# Assign files to bash variables
REFERENCE_SEQUENCE=$1
LEFT_READS=$2
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAMPLE=$3
SAM_FILE=$4

# Align reads with bwa
bwa mem \
-M \
-t 8 \
-R "@RG\tID:$SAMPLE\tPL:illumina\tPU:$SAMPLE\tSM:$SAMPLE" \
$REFERENCE_SEQUENCE \
$LEFT_READS \
$RIGHT_READS \
-o $SAM_FILE
```

#### Picard processing

Similarly to the previous scripts, we will just need to remove our `SBATCH` directives and edit our variables to automate the behavior of the script. Let's start by creating a new script:

```
cp picard_alignment_processing_normal.sbatch picard_alignment_processing_automated.sbatch
vim picard_alignment_processing_automated.sbatch
```

Once again, we will remove the `sbatch` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o picard_alignment_processing_normal_%j.out
#SBATCH -e picard_alignment_processing_normal_%j.err
```

Next, we will change the following lines to accept positional parameters by changing:

```
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.sam
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/picard/syn3_normal/
SAMPLE_NAME=syn3_normal
```

To:

```
SAM_FILE=$1
REPORTS_DIRECTORY=$2
SAMPLE_NAME=$3
```

Now your `Picard` alignment processing script should look like:

```
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Picard 

# Load module
module load picard/2.27.5

# Assign file paths to variables
SAM_FILE=$1
REPORTS_DIRECTORY=$2
SAMPLE_NAME=$3
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=${REPORTS_DIRECTORY}/${SAMPLE_NAME}.remove_duplicates_metrics.txt
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`

# Make reports directory
mkdir -p $REPORTS_DIRECTORY

# Query-sort alginment file and convert to BAM
java -jar $PICARD/picard.jar SortSam \
--INPUT $SAM_FILE \
--OUTPUT $QUERY_SORTED_BAM_FILE \
--SORT_ORDER queryname

# Mark and remove duplicates
java -jar $PICARD/picard.jar MarkDuplicates \
--INPUT $QUERY_SORTED_BAM_FILE \
--OUTPUT $REMOVE_DUPLICATES_BAM_FILE \
--METRICS_FILE $METRICS_FILE \
--REMOVE_DUPLICATES true

# Coordinate-sort BAM file and create BAM index file
java -jar $PICARD/picard.jar SortSam \
--INPUT $REMOVE_DUPLICATES_BAM_FILE \
--OUTPUT $COORDINATE_SORTED_BAM_FILE \
--SORT_ORDER coordinate \
--CREATE_INDEX true
```

#### Picard Metrics

Now we will conver the `Picard` metrics step:

```
cp picard_metrics_normal.sbatch picard_metrics_automated.sbatch
vim picard_metrics_automated.sbatch
```

Remove the `SBATCH` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o picard_metrics_normal_%j.out
#SBATCH -e picard_metrics_normal_%j.err
```

Then change the `bash` variables from:

```
# Assign variables
INPUT_BAM=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa
OUTPUT_METRICS_FILE=/home/${USER}/variant_calling/reports/picard/syn3_normal/syn3_normal_GRCh38.p7.CollectAlignmentSummaryMetrics.txt
```

To:

```
INPUT_BAM=$1
REFERENCE=$2
OUTPUT_METRICS_FILE=$3
```

The automated `Picard` metrics submission script should now look like:

```
#!/bin/bash
# This sbatch script is for collecting alignment metrics using Picard 

# Load picard
module load picard/2.27.5

# Assign variables
INPUT_BAM=$1
REFERENCE=$2
OUTPUT_METRICS_FILE=$3

# Run Picard CollectAlignmentSummaryMetrics
java -jar $PICARD/picard.jar CollectAlignmentSummaryMetrics \
--INPUT $INPUT_BAM \
--REFERENCE_SEQUENCE $REFERENCE \
--OUTPUT $OUTPUT_METRICS_FILE
```

#### MultiQC

Let's create our new `MultiQC` submission script:

```
cp multiqc_alignment_metrics_normal_tumor.sbatch multiqc_alignment_metrics_automated.sbatch 
vim multiqc_alignment_metrics_automated.sbatch
```

Remove the `SBATCH` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:10:00
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -o multiqc_alignment_metrics_%j.out
#SBATCH -e multiqc_alignment_metrics_%j.err
```

Change the following variables to accept positional parameters from:

```
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/
NORMAL_SAMPLE_NAME=syn3_normal
TUMOR_SAMPLE_NAME=syn3_tumor
REFERENCE=GRCh38.p7
```

To:

```
REPORTS_DIRECTORY=$1
NORMAL_SAMPLE_NAME=$2
TUMOR_SAMPLE_NAME=$3
REFERENCE=$4
```

The final `MultiQC` automation submission script should look like:

```
#!/bin/bash
# This sbatch script is for collating alignment metrics from Picard using MultiQC

# Load modules
module load gcc/9.2.0
module load multiqc/1.12

# Assign variables
REPORTS_DIRECTORY=$1
NORMAL_SAMPLE_NAME=$2
TUMOR_SAMPLE_NAME=$3
REFERENCE=$4
NORMAL_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_${REFERENCE}.CollectAlignmentSummaryMetrics.txt
TUMOR_PICARD_METRICS=${REPORTS_DIRECTORY}picard/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_${REFERENCE}.CollectAlignmentSummaryMetrics.txt
NORMAL_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_1_fastqc.zip
NORMAL_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${NORMAL_SAMPLE_NAME}/${NORMAL_SAMPLE_NAME}_2_fastqc.zip
TUMOR_FASTQC_1=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_1_fastqc.zip
TUMOR_FASTQC_2=${REPORTS_DIRECTORY}fastqc/${TUMOR_SAMPLE_NAME}/${TUMOR_SAMPLE_NAME}_2_fastqc.zip
OUTPUT_DIRECTORY=${REPORTS_DIRECTORY}/multiqc/

# Create directory for output
mkdir -p $OUTPUT_DIRECTORY

# Run MultiQC
multiqc \
$NORMAL_PICARD_METRICS \
$TUMOR_PICARD_METRICS \
$NORMAL_FASTQC_1 \
$NORMAL_FASTQC_2 \
$TUMOR_FASTQC_1 \
$TUMOR_FASTQC_2 \
--outdir $OUTPUT_DIRECTORY
```

#### Variant calling with MuTect2

Let's make a copy of the `MuTect2` script that we wrote so that we can adapt it for our automation.

```
cp mutect2_normal_tumor.sbatch mutect2_automated.sbatch
vim mutect2_automated.sbatch
```

Once again, we need to remove the `sbatch` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o mutect2_variant_calling_%j.out
#SBATCH -e mutect2_variant_calling_%j.err
```

Next, change the rest of the variable lines from:

```
# Assign variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_SAMPLE_NAME=syn3_normal
NORMAL_BAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${NORMAL_SAMPLE_NAME}_GRCh38.p7.coordinate_sorted.bam
TUMOR_SAMPLE_NAME=syn3_tumor
TUMOR_BAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${TUMOR_SAMPLE_NAME}_GRCh38.p7.coordinate_sorted.bam
VCF_OUTPUT_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/mutect2_${NORMAL_SAMPLE_NAME}_${TUMOR_SAMPLE_NAME}_GRCh38.p7-raw.vcf
```

To:

```
# Assign variables
REFERENCE_SEQUENCE=$1
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_SAMPLE_NAME=$2
NORMAL_BAM_FILE=$3
TUMOR_SAMPLE_NAME=$4
TUMOR_BAM_FILE=$5
VCF_OUTPUT_FILE=$6
```

The final automated `sbatch` script for `MuTect2` should look like:

```
#!/bin/bash
# This sbatch script is for variant calling with GATK's MuTect2

# Load the GATK module
module load gatk/4.1.9.0

# Assign variables
REFERENCE_SEQUENCE=$1
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_SAMPLE_NAME=$2
NORMAL_BAM_FILE=$3
TUMOR_SAMPLE_NAME=$4
TUMOR_BAM_FILE=$5
VCF_OUTPUT_FILE=$6

# Run MuTect2
gatk Mutect2 \
--sequence-dictionary $REFERENCE_DICTIONARY \
--reference $REFERENCE_SEQUENCE \
--input $NORMAL_BAM_FILE \
--normal-sample $NORMAL_SAMPLE_NAME \
--input $TUMOR_BAM_FILE \
--tumor-sample $TUMOR_SAMPLE_NAME \
--annotation ClippingRankSumTest --annotation DepthPerSampleHC --annotation MappingQualityRankSumTest --annotation MappingQualityZero --annotation QualByDepth --annotation ReadPosRankSumTest --annotation RMSMappingQuality --annotation FisherStrand --annotation MappingQuality --annotation DepthPerAlleleBySample --annotation Coverage \
--output $VCF_OUTPUT_FILE
```

#### Variant Filtering

Next, we can automate our variant filtering script. First, make a copy of our variant filtering `sbatch` script and open it in `vim`:

```
cp variant_filtering_normal_tumor.sbatch variant_filtering_automated.sbatch 
vim variant_filtering_automated.sbatch 
```

Once again we will remove the `sbatch` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:15:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_filtering_%j.out
#SBATCH -e variant_filtering_%j.err
```

Next we need to alter some of our variables from:

```
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
RAW_VCF_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf.gz
LCR_FILE=/n/groups/hbctraining/variant_calling/reference/LCR-hs38.bed
```

To:

```
REFERENCE_SEQUENCE=$1
RAW_VCF_FILE=$2
LCR_FILE=$3
```

The final automated `sbatch` should look like:

```
#!/bin/bash
# This sbatch script is for variant filtering 

# Load modules
module load gatk/4.1.9.0
module load snpEff/4.3g

# Assign variables
REFERENCE_SEQUENCE=$1
RAW_VCF_FILE=$2
LCR_FILE=$3
MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf}filt.vcf
PASSING_FILTER_VCF=${RAW_VCF_FILE%raw.vcf}pass-filt.vcf
LCR_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf}pass-filt-LCR.vcf

# Filter Mutect Calls
gatk FilterMutectCalls \
--reference $REFERENCE_SEQUENCE \
--variant $RAW_VCF_FILE \
--output $MUTECT_FILTERED_VCF

# Filter for only SNPs with PASS in the FILTER field
java -jar $SNPEFF/SnpSift.jar filter \
-noLog \
"( FILTER = 'PASS' )" \
$MUTECT_FILTERED_VCF > $PASSING_FILTER_VCF

# Filter LCR
java -jar $SNPEFF/SnpSift.jar intervals \
-noLog \
-x \
-i $PASSING_FILTER_VCF \
$LCR_FILE > $LCR_FILTERED_VCF
```

#### Variant Annotation

Make a copy of our variant annotation `sbatch` submission script and open it up in `vim`:

```
cp variant_annotation_normal_tumor.sbatch variant_annotation_automated.sbatch
vim variant_annotation_automated.sbatch
```

We will need to remove the `sbatch` directives:

```
# REMOVE THESE LINES
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_annotation_%j.out
#SBATCH -e variant_annotation_%j.err
```

Next, change the variables from:

```
# Assign variables
REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/snpeff/
SAMPLE_NAME=mutect2_syn3_normal_syn3_tumor
REFERENCE_SEQUENCE_NAME=GRCh38.p7
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=GRCh38.p7.RefSeq
DATADIR=/n/groups/hbctraining/variant_calling/reference/snpeff/data/
FILTERED_VCF_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-pass-filt-LCR.vcf
PEDIGREE_HEADER_FILE=/home/$USER/variant_calling/scripts/syn3_normal_syn3_tumor_pedigree_header.txt
FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER=${FILTERED_VCF_FILE%.vcf}.pedigree_header.vcf
SNPEFF_ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER%.vcf}.snpeff.vcf
DBSNP_DATABASE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.dbSNP.vcf.gz
DBSNP_ANNOTATED_VCF_FILE=${SNPEFF_ANNOTATED_VCF_FILE%.vcf}.dbSNP.vcf
```

To:

```
# Assign variables
REPORTS_DIRECTORY=$1
SAMPLE_NAME=$2
REFERENCE_SEQUENCE_NAME=$3
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=$4
DATADIR=$5
FILTERED_VCF_FILE=$6
PEDIGREE_HEADER_FILE=$7
FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER=${FILTERED_VCF_FILE%.vcf}.pedigree_header.vcf
SNPEFF_ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER%.vcf}.snpeff.vcf
DBSNP_DATABASE=$8
DBSNP_ANNOTATED_VCF_FILE=${SNPEFF_ANNOTATED_VCF_FILE%.vcf}.dbSNP.vcf
```

Our automated `sbatch` submission for variant annotation should look like:

```
#!/bin/bash
# Using SnpEff to annotate our variants

# Load modules
module load gcc/9.2.0
module load bcftools/1.14
module load snpEff/4.3g

# Assign variables
REPORTS_DIRECTORY=$1
SAMPLE_NAME=$2
REFERENCE_SEQUENCE_NAME=$3
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=$4
DATADIR=$5
FILTERED_VCF_FILE=$6
PEDIGREE_HEADER_FILE=$7
FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER=${FILTERED_VCF_FILE%.vcf}.pedigree_header.vcf
SNPEFF_ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER%.vcf}.snpeff.vcf
DBSNP_DATABASE=$8
DBSNP_ANNOTATED_VCF_FILE=${SNPEFF_ANNOTATED_VCF_FILE%.vcf}.dbSNP.vcf

# Create reports directory
mkdir -p $REPORTS_DIRECTORY

# Append Header
bcftools annotate \
--header-lines $PEDIGREE_HEADER_FILE \
$FILTERED_VCF_FILE \
> $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER

# Run SnpEff
java -jar -Xmx4g $SNPEFF/snpEff.jar  eff \
-dataDir $DATADIR \
-cancer \
-noLog \
-csvStats $CSV_STATS \
-s $HTML_REPORT \
$REFERENCE_DATABASE \
$FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER \
> $SNPEFF_ANNOTATED_VCF_FILE

# Use dbSNP VCF to annotate our VCF
java -jar $SNPEFF/SnpSift.jar annotate \
$DBSNP_DATABASE \
-tabix \
-noLog \
$SNPEFF_ANNOTATED_VCF_FILE \
> $DBSNP_ANNOTATED_VCF_FILE
```

### Developing a Wrapper Script

Now that we have each of our individual `sbatch` submission scripts written, we can tie them all together in a wrapper script. Before we do this though, we are going to briefly discuss the  `--dependency` option in the `sbatch` command. The `--dependency` option is very helpful for a few cases:

**1)** When designing a pipeline to run on a cluster where different parts of the pipeline have different computational requirements.

**2)** When running a pipeline where you expect the output from one job to finish at a time when you won't be availible to submit it as input to the next step in the pipeline.

**3)** When multiple parts of the pipeline all need to be completed before the next step can begin.

#### Making one sbatch submission dependent on another

Using the `--dependency` option in `sbatch` will likely make the execution of your pipelines faster and also help you be a better citizen of the cluster because each part of your pipeline can have customized computational requests. So how does the `--dependency` option work? Let's imagine that we have two jobs `Job_A.sbatch` and `Job_B.sbatch`. Additionally, we want to use the output from `Job_A.sbatch` as input for `Job_B.sbatch`. First, we would submit `Job_A.sbatch` to the cluster like normal:

```
# Don't run this, it is an example
sbatch Job_A.sbatch
```

And once it has been submitted it should return the job ID with text like:

`Submitted batch job 10928`

Where `10928` represent your unique job ID for the job.

You can also find this job ID by looking in your `squeue` like:

```
squeue -u $USER
```

This will show the actuve jobs running and it will show you the job IDs for the jobs that you currently have running.

Now with out job ID for `Job_A.sbatch` in hand, we can submit `Job_B.sbatch`

```
# Don't run this, it is an example
sbatch --dependency=afterok:10928 Job_B.sbatch
```

There are different types of dependencies besides `afterok`, but their usage is rare. However, `afterok` tells SLURM that as longs as job ID `10928` ended without an error, then go ahead and put `Job_B.sbatch` in the queue.

#### Making one sbatch submission dependent multiple others

Now let's imagine we have a case that we will experience shortly. We have a script (`Job_C.sbatch`) that is dependent on two job IDs (`Job_A.sbatch` and `Job_B.sbatch`), as is the case with `GATK` waiting for `Picard` metrics for the normal and tumor samples to finish. In this case, each job ID just needs to be separated by a `:` like:

```
$ sbatch Job_A.sbatch
Submitted batch job 10928
$ sbatch Job_B.sbatch
Submitted batch job 10931
$ sbatch --dependency=afterok:10928:10931 Job_C.sbatch
```

In this case, `Job_C.sbatch` won't queue until `Job_A.sbatch` and `Job_B.sbatch` exit without an error. Furthermore, you can queue a dependent job on a job that is dependent, thus creating an entire pipeline of dependent jobs.

> NOTE: The behavior on O2 is that if any part of your pipeline gets an error, then the pipeline will stop and the downstream dependent jobs will be cancelled. Some high-performance computing clusters will not cancel the downstream dependent jobs though and you will need to manually remove them from the queue with `scancel`.

#### Writing the Wrapper script

Let's start by opening a file in `vim` called `variant_calling_wrapper.sh`:

```
vim variant_calling_wrapper.sh
```

Inside of this file copy and paste the following wrapper: 

```
#!/bin/bash
# This is the wrapper script for the variant calling pipeline. 
# It is designed for a single individual, so if you wanted to run this across multiple individuals you would need to loop this script across all of them.

# Assign variables
FASTQ_DIRECTORY=/home/$USER/variant_calling/raw_data/
NORMAL_SAMPLE=syn3_normal
TUMOR_SAMPLE=syn3_tumor
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.fa
REFERENCE_SEQUENCE_NAME=`basename $REFERENCE_SEQUENCE .fa`
LCR_FILE=/n/groups/hbctraining/variant_calling/reference/LCR-hs38.bed
REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/
PEDIGREE_HEADER_FILE=/home/$USER/variant_calling/scripts/${NORMAL_SAMPLE}_${TUMOR_SAMPLE}_pedigree_header.txt
SNPEFF_DIRECTORY=/n/groups/hbctraining/variant_calling/reference/snpeff/data/
SNPEFF_DATABASE=GRCh38.p7.RefSeq
DBSNP_DATABASE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.dbSNP.vcf.gz

# Intiate bash arrays to hold sample names, Picard Job IDs and coordinate-sorted BAM files
SAMPLE_NAME_ARRAY=()
PICARD_METRICS_JOB_ID_ARRAY=()
COORDINATE_SORTED_BAM_ARRAY=()

# Submit samples for Alignment and Picard Processing
# For each Read 1 file do:
for SAMPLE in $FASTQ_DIRECTORY*_1.fq.gz; do
  # Parse out the sample name
  SAMPLE_NAME=`basename $SAMPLE _1.fq.gz`
  # Add the sample name to our array holding sample names
  SAMPLE_NAME_ARRAY+=($SAMPLE_NAME)
  # Assign FastQC report directory
  FASTQC_REPORT_DIRECTORY=`echo -e "${REPORTS_DIRECTORY}fastqc/${SAMPLE_NAME}/"`
  # Submit the FastQC script
  FASTQC_JOB_SUBMISSION=$(sbatch -p priority -t 0-00:30:00 -c 4 --mem 8G -o fastqc_${SAMPLE_NAME}_%j.out -e fastqc_${SAMPLE_NAME}_%j.err fastqc_automated.sbatch $SAMPLE $FASTQC_REPORT_DIRECTORY 4)
  # Parse out the job ID from outout from the FastQC submission
  FASTQC_JOB_ID=`echo $FASTQC_JOB_SUBMISSION | cut -d ' ' -f 4`
  # Print to standard output the job that has been submitted
  echo -e "FastQC job for sample $SAMPLE_NAME submitted as job ID $FASTQC_JOB_ID"
  # Assign a path and name for the alignments
  SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.sam
  # Submit the bwa sbatch script and save the output to a variable named $BWA_JOB_SUBMISSION
  BWA_JOB_SUBMISSION=$(sbatch -p priority -t 0-04:00:00 -c 8 --mem 16G -o bwa_alignment_${SAMPLE_NAME}_%j.out -e bwa_alignment_${SAMPLE_NAME}_%j.err --dependency=afterok:$FASTQC_JOB_ID bwa_alignment_automated.sbatch $REFERENCE_SEQUENCE $SAMPLE $SAMPLE_NAME $SAM_FILE)
  # Parse out the job ID from outout from the bwa submission
  BWA_JOB_ID=`echo $BWA_JOB_SUBMISSION | cut -d ' ' -f 4`
  # Print to standard output the job that has been submitted
  echo -e "bwa job for sample $SAMPLE_NAME submitted as job ID $BWA_JOB_ID"
  # Assign Picard report directory
  PICARD_METRICS_REPORTS_DIRECTORY=`echo -e "${REPORTS_DIRECTORY}picard/${SAMPLE_NAME}/"`
  # Submit the picard sbatch script and save the output to a variable named $PICARD_JOB_SUBMISSION
  PICARD_PROCESSING_JOB_SUBMISSION=$(sbatch -p priority -t 0-04:00:00 -c 1 --mem 32G -o picard_alignment_processing_${SAMPLE_NAME}_%j.out -e picard_alignment_processing_${SAMPLE_NAME}_%j.err --dependency=afterok:$BWA_JOB_ID picard_alignment_processing_automated.sbatch $SAM_FILE $PICARD_METRICS_REPORTS_DIRECTORY $SAMPLE_NAME)
  # Parse out the job ID from output from the Picard processing submission
  PICARD_PROCESSING_JOB_ID=`echo $PICARD_PROCESSING_JOB_SUBMISSION | cut -d ' ' -f 4`
  # Print to standard output the job that has been submitted
  echo -e "Picard processing job for sample $SAMPLE_NAME submitted as job ID $PICARD_PROCESSING_JOB_ID"
  # Create a variable to hold the the path and output BAM file from the Picard processing sbatch script
  COORDINATE_SORTED_BAM_FILE=`echo -e "${SAM_FILE%sam}coordinate_sorted.bam"`
  # Add this BAM file path to an array
  COORDINATE_SORTED_BAM_ARRAY+=($COORDINATE_SORTED_BAM_FILE)
  # Assign Picard report file path
  PICARD_METRICS_REPORT_FILE=`echo -e "${PICARD_METRICS_REPORTS_DIRECTORY}${SAMPLE_NAME}.CollectAlignmentSummaryMetrics.txt"`
  #Submit the picard metrics sbatch script and save the output to a variable named $PICARD_METRICS_JOB_SUBMISSION
  PICARD_METRICS_JOB_SUBMISSION=$(sbatch -p priority -t 0-00:30:00 -c 1 --mem 16G -o picard_metrics_${SAMPLE_NAME}_%j.out -e picard_metrics_${SAMPLE_NAME}_%j.err --dependency=afterok:$PICARD_PROCESSING_JOB_ID picard_metrics_automated.sbatch $COORDINATE_SORTED_BAM_FILE $REFERENCE_SEQUENCE $PICARD_METRICS_REPORT_FILE)
  # Parse out the job ID from the Picard metrics submission
  PICARD_METRICS_JOB_ID=`echo $PICARD_METRICS_JOB_SUBMISSION | cut -d ' ' -f 4`
  # Add the Picard metrics job ID to an array holding all of the Picard metrics job IDs
  PICARD_METRICS_JOB_ID_ARRAY+=($PICARD_METRICS_JOB_ID)
  # Print to standard output the job that has been submitted
  echo -e "Picard metrics job for sample $SAMPLE_NAME submitted as job ID $PICARD_METRICS_JOB_ID"
done

# Create a string that contains all of the samples preceded by and separated by an "_"
SAMPLE_NAME_STRING=`for i in ${SAMPLE_NAME_ARRAY[@]}; do  echo -ne "_$i"; done`
# Create a string that contains all of the Picard procession job IDs preceded by and separated by a ":"
DEPENDENT_PICARD_METRICS_JOB_IDS=`for i in ${PICARD_METRICS_JOB_ID_ARRAY[@]}; do  echo -ne ":$i"; done`

# for loop to discover which position in the $SAMPLE_NAME_ARRAY contains the normal sample
NORMAL_ARRAY_POSITION=`for SAMPLE_POSITION in {0..1}; do if [[ "${SAMPLE_NAME_ARRAY[$SAMPLE_POSITION]}" == "$NORMAL_SAMPLE" ]];then echo $SAMPLE_POSITION; fi; done`
# for loop to discover which position in the $SAMPLE_NAME_ARRAY contains the tumor sample
TUMOR_ARRAY_POSITION=`for SAMPLE_POSITION in {0..1}; do if [[ "${SAMPLE_NAME_ARRAY[$SAMPLE_POSITION]}" == "$TUMOR_SAMPLE" ]];then echo $SAMPLE_POSITION; fi; done`

# Assign variables that will be used by Mutect2 
NORMAL_BAM_MUTECT_INPUT=${COORDINATE_SORTED_BAM_ARRAY[$NORMAL_ARRAY_POSITION]}
TUMOR_BAM_MUTECT_INPUT=${COORDINATE_SORTED_BAM_ARRAY[$TUMOR_ARRAY_POSITION]}
MUTECT2_VCF_OUTPUT=`echo -e "/n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/mutect2${SAMPLE_NAME_STRING}_${REFERENCE_SEQUENCE_NAME}-raw.vcf"`

# Submit the Mutect2 sbatch script and save the output to a variable named $MUTECT2_JOB_SUBMISSION 
MUTECT2_JOB_SUBMISSION=$(sbatch -p priority -t 1-00:00:00 -c 1 --mem 16G -o mutect2_variant_calling${SAMPLE_NAME_STRING}_%j.out -e mutect2_variant_calling${SAMPLE_NAME_STRING}_%j.err --dependency=afterok${DEPENDENT_PICARD_METRICS_JOB_IDS} mutect2_automated.sbatch  $REFERENCE_SEQUENCE $NORMAL_SAMPLE $NORMAL_BAM_MUTECT_INPUT $TUMOR_SAMPLE $TUMOR_BAM_MUTECT_INPUT $MUTECT2_VCF_OUTPUT)

# Parse out the job ID from output from the Mutect2 submission
MUTECT2_JOB_ID=`echo $MUTECT2_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "Mutect2 job submitted as job ID $MUTECT2_JOB_ID"

# Submit MultiQC sbatch script and save the output to a variable named $MULTIQC_JOB_SUBMISSION
MULTIQC_JOB_SUBMISSION=$(sbatch -p priority -t 0-00:10:00 -c 1 --mem 1G -o multiqc${SAMPLE_NAME_STRING}_%j.out -e multiqc${SAMPLE_NAME_STRING}_%j.err --dependency=afterok${DEPENDENT_PICARD_METRICS_JOB_IDS} multiqc_alignment_metrics_automated.sbatch $REPORTS_DIRECTORY $NORMAL_SAMPLE $TUMOR_SAMPLE $REFERENCE_SEQUENCE_NAME)

# Parse out the job ID from output from the Mutect2 submission
MULTIQC_JOB_ID=`echo $MULTIQC_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "MultiQC job submitted as job ID $MULTIQC_JOB_ID"

# Submit the variant filtering sbatch script and save the output to a variable named $VARIANT_FILTERING_JOB_SUBMISSION
VARIANT_FILTERING_JOB_SUBMISSION=$(sbatch -p priority -t 0-00:10:00 -c 1 --mem 8G -o variant_filtering${SAMPLE_NAME_STRING}_%j.out -e variant_filtering${SAMPLE_NAME_STRING}_%j.err --dependency=afterok:$MUTECT2_JOB_ID variant_filtering_automated.sbatch $REFERENCE_SEQUENCE $MUTECT2_VCF_OUTPUT $LCR_FILE)

# Parse out the job ID from output from the variant filtering submission
VARIANT_FILTERING_JOB_ID=`echo $VARIANT_FILTERING_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "Variant filtering job submitted as job ID $VARIANT_FILTERING_JOB_ID"

# Assign a variable that was used in the Mutect2 filtering step
MUTECT2_VCF_OUTPUT_FILTERED=`echo -e "${MUTECT2_VCF_OUTPUT%raw.vcf}pass-filt-LCR.vcf"`

# Create PEDIGREE header file 
echo -e "##PEDIGREE=<Derived=${TUMOR_SAMPLE},Original=${NORMAL_SAMPLE}>" > $PEDIGREE_HEADER_FILE

# Create SnpEff reports directory
SNPEFF_REPORTS_DIRECTORY=`echo -e "${REPORTS_DIRECTORY}snpeff/"`

# Submit the variant annotation sbatch script
VARIANT_ANNOTATION_JOB_SUBMISSION=$(sbatch -p priority -t 0-02:00:00 -c 1 --mem 8G -o variant_annotation${SAMPLE_NAME_STRING}_%j.out -e variant_annotation${SAMPLE_NAME_STRING}_%j.err --dependency=afterok:$VARIANT_FILTERING_JOB_ID variant_annotation_automated.sbatch $SNPEFF_REPORTS_DIRECTORY $SAMPLE_NAME_STRING $REFERENCE_SEQUENCE_NAME $SNPEFF_DATABASE $SNPEFF_DIRECTORY $MUTECT2_VCF_OUTPUT_FILTERED $PEDIGREE_HEADER_FILE $DBSNP_DATABASE)

# Parse out the job ID from output from the variant annotation submission
VARIANT_ANNOTATION_JOB_ID=`echo $VARIANT_ANNOTATION_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "Variant annotation job submitted as job ID $VARIANT_ANNOTATION_JOB_ID"
```

If we were to type:

```
### DON'T DO THIS
sh variant_calling_wrapper.sh
```

This it would execute our automated pipeline. But since we already have the files we neeed, we don't want to do this.

## Exrecises

**1.** If we submit `Job_A.sh` and `SLURM` returns:

`Submitted batch job 213489`

And we want to start `Job_B.sh` after `Job_A.sh` finishes without error, what command would we use to do this?

**2.** If we submit `Job_X.sh` and `SLURM` returns:

`Submitted batch job 213489`

Then we submit `Job_Y.sh` and `SLURM` returns:

`Submitted batch job 213496`

And we want to start `Job_Z.sh` after `Job_X.sh` and `Job_Y.sh` finish without error, what command would we use to do this?

**3.** If you want to submit `Job_M.sh`, and provide `reference.fa` as the first positional parameter and `input.bam` as the second positional parameter, how could you do this?

***

[Next Lesson >>](13_variant_prioritization.md)

[Back to Schedule](../schedule/README.md)
