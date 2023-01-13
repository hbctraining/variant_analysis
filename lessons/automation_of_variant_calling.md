# Automation of Variant Calling

## Learning objectives
- Construct a flexible pipeline for automating variant calling
- Integrate the `--dependency` option for `sbatch` into workflows

## Automating the Variant Calling Workflow

Now that we have completed most of the variant calling workflow, we would like to explore ways to automate it to make future runs identical and also minimize downtime in the workflow. Automating this workflow with consist of three parts:

1) Editing our `sbatch` submission scripts to accept positional parameters for variable
2) Creating a wrapper script to execute each script with the appropriate parameters
3) Using the `--dependency` option in `sbatch` to construct dependent pipeline for scripts that have diffent computational requirements

The workflow and methodolgies we used up to this point are completely fine and could work for you. However, you may have many scripts that you need to create if you have many samples and that can be cumbersome, so we are going to introduce this way of automating the process.

### Editing our `sbatch` Submission Scripts

Currently, within our `sbatch` submission scripts, we have `bash` variables set to hold various files. Instead of explicitly defining those variables now as a given file, we are going to use positional parameters to allow ourselves to plug-in various files. First let's navigate to our scripts directory:

```
cd ~/variant_calling/scripts/
```

#### `bwa` Alignment

Because many of some variables already weren't explictly defined (`$RIGHT_READS` and `$SAM_FILE`), but were rather dependent on our explicitly defined variables (`$REFERENCE_SEQUENCE` and `LEFT_READS`) we won't need to edit our script as much. We are able to do this because `sbatch` accepts positional parameters. So let's create a new script and make this edit.

```
cp bwa_alignment_normal.sbatch bwa_alignment_automated.sbatch
vim bwa_alignment_automated.sbatch
```

In order for us to specify "tumor" or "normal" samples in the standard error and standard output files we will need to pass `sbatch` these arguments directly when submit the jobs in the wrapper. However, SLURM sometimes has issues when some arguments come from the command line and others come from the script. As a result, we are going remove all `sbatch` directives from scripts and rather provide them all in the command-line as part of the wrapper. Therefore, we need to remove the following lines of `sbatch` directives:

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

And we will also need to change the variables:

```
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
REFERENCE_SEQUENCE_NAME=`basename $REFERENCE_SEQUENCE _genomic.fa`
SAMPLE_NAME=`basename $LEFT_READS _1.fq.gz`
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.sam
```

And we will set these to take positional parameters:

```
REFERENCE_SEQUENCE=$1
LEFT_READS=$2
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=$3
```

Now our `bwa` submission script should look like:

```
#!/bin/bash
# This script is the automated version for aligning sequencing reads against a reference genome using bwa

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17

# Assign files to bash variables
REFERENCE_SEQUENCE=$1
LEFT_READS=$2
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=$3

# Align reads with bwa
bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-normal\tPL:illumina\tPU:syn3-normal\tSM:syn3-normal' \
$REFERENCE_SEQUENCE \
$LEFT_READS \
$RIGHT_READS \
-o $SAM_FILE
```

#### `Picard` processing

Because `Picard` and `bwa` have different computational requirements (i.e. different numbers of CPUs, memory), we are going to need to automate the `sbatch` submission scripts separately and then use a wrapper script to tie them all together. Similarly to the `bwa` script, we will just need to edit our variables to automate the behavior of the script. Let's start by creating a new script:

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

Because all of the variables we set except `$SAM_FILE` are dependent on the value of `$SAM_FILE`, we only need to change this to be a positional parameter.

Change the following line:

```
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.sam
```

To:

```
SAM_FILE=$1
```

Now your `Picard` alignment processing script should look like:

```
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Picard 

# Assign file paths to variables 
SAM_FILE=$1
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=`echo ${SAM_FILE%sam}remove_duplicates_metrics.txt`
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`

module load picard/2.8.0

java -jar $PICARD/picard-2.8.0.jar SortSam \
INPUT=$SAM_FILE \
OUTPUT=$QUERY_SORTED_BAM_FILE \
SORT_ORDER=queryname

java -jar $PICARD/picard-2.8.0.jar MarkDuplicates \
INPUT=$QUERY_SORTED_BAM_FILE \
OUTPUT=$REMOVE_DUPLICATES_BAM_FILE \
METRICS_FILE=$METRICS_FILE \
REMOVE_DUPLICATES=true

java -jar $PICARD/picard-2.8.0.jar SortSam \
INPUT=$REMOVE_DUPLICATES_BAM_FILE \
OUTPUT=COORDINATE_SORTED_BAM_FILE \
SORT_ORDER=coordinate \
CREATE_INDEX=true
```

#### Variant calling with GATK

```
cp mutect2_normal_tumor.sbatch mutect2_automated.sbatch
```

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

module load gatk/4.1.9.0

REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_BAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_tumor_GRCh38.p7.coordinate_sorted.bam
NORMAL_SAMPLE_NAME=syn3-normal
TUMOR_BAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.coordinate_sorted.bam
TUMOR_SAMPLE_NAME=syn3-tumor
VCF_OUTPUT_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf.gz

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


```
#!/bin/bash
# This sbatch script is for variant calling with GATK's MuTect2

module load gatk/4.1.9.0

REFERENCE_SEQUENCE=$1
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_BAM_FILE=$2
NORMAL_SAMPLE_NAME=$3
TUMOR_BAM_FILE=$4
TUMOR_SAMPLE_NAME=$5
VCF_OUTPUT_FILE=$6

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

#### Variant Filtering

```
#!/bin/bash
# This sbatch script is for variant filtering 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:15:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_filtering_%j.out
#SBATCH -e variant_filtering_%j.err

module load gatk/4.1.9.0
module load snpEff/4.3g

REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
RAW_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_normal_syn3_tumor_GRCh38.p7-raw.vcf.gz
LCR_FILE=/n/groups/hbctraining/variant_calling/reference/LCR-hs38.bed

MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}filt.vcf.gz
LCR_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}LCR-filt.vcf

gatk FilterMutectCalls \
--reference  $REFERENCE_SEQUENCE \
--variant $RAW_VCF_FILE \
--output $MUTECT_FILTERED_VCF

java -jar $SNPEFF/SnpSift.jar intervals -x -i $MUTECT_FILTERED_VCF $LCR_FILE > $LCR_FILTERED_VCF
```

```
cp variant_filtering.sbatch variant_filtering_automated.sbatch 
vim variant_filtering_automated.sbatch 
```

Once again we need to remove the `sbatch` directives:

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


```
#!/bin/bash
# This sbatch script is for variant filtering 

module load gatk/4.1.9.0
module load snpEff/4.3g

REFERENCE_SEQUENCE=$1
RAW_VCF_FILE=$2
LCR_FILE=$3

MUTECT_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}filt.vcf.gz
LCR_FILTERED_VCF=${RAW_VCF_FILE%raw.vcf.gz}LCR-filt.vcf

gatk FilterMutectCalls \
--reference $REFERENCE_SEQUENCE \
--variant $RAW_VCF_FILE \
--output $MUTECT_FILTERED_VCF

java -jar $SNPEFF/SnpSift.jar intervals -x -i $MUTECT_FILTERED_VCF $LCR_FILE > $LCR_FILTERED_VCF
```

#### Variant Annotation

```
#!/bin/bash
# This sbatch script is for variant annotation 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_annotation_%j.out
#SBATCH -e variant_annotation_%j.err

module load snpEff/4.3g

CSV_STATS=/home/$USER/variant_calling/reports/syn3_GRCh38.p7-effects-stats.csv
HTML_REPORT=/home/$USER/variant_calling/reports/syn3_GRCh38.p7-effects-stats.html
REFERENCE_DATABASE=GRCh38.p7.RefSeq
FILTERED_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/syn3_GRCh38.p7-LCR-filt.vcf
ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE%vcf}snpeff.vcf

java -jar $SNPEFF/snpEff.jar  eff \
-dataDir /n/groups/shared_databases/snpEff.data/ \
-cancer \
-noLog \
-csvStats $CSV_STATS\
-s  $HTML_REPORT \
$REFERENCE_DATABASE \
$FILTERED_VCF_FILE > $ANNOTATED_VCF_FILE
```
```
cp variant_annotation.sbatch variant_annotation_automated.sbatch
vim variant_annotation_automated.sbatch
```

For the last time we will need to remove the `sbatch` directives:

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

```
#!/bin/bash
# This sbatch script is for variant annotation 

module load snpEff/4.3g

CSV_STATS=$1
HTML_REPORT=$2
REFERENCE_DATABASE=$3
FILTERED_VCF_FILE=$4
ANNOTATED_VCF_FILE=$5

java -jar $SNPEFF/snpEff.jar  eff \
-dataDir /n/groups/shared_databases/snpEff.data/ \
-cancer \
-noLog \
-csvStats $CSV_STATS\
-s  $HTML_REPORT \
$REFERENCE_DATABASE \
$FILTERED_VCF_FILE > $ANNOTATED_VCF_FILE
```


### Developing a Wrapper Script

Now that we have each of our individual `sbatch` submission scripts written, we can tie them all together in a wrapper script. Before we do this though, we are going to briefly discuss the  `--dependency` option in the `sbatch` command. The `--dependency` option is very helpful for a few cases:

1) When designing a pipeline to run on a cluster where different parts of the pipeline have different computational requirements 
2) When running a pipeline where you expect the output from one job to finish at a time when you won't be availible to submit it as input to the next step in the pipeline
3) When multiple parts of the pipeline all need to be completed before the next step can begin

#### Making one sbatch submission dependent on another

Using the `--dependency` option in `sbatch` will likely make the execution of your pipelines faster and also help you be a better citizen of the cluster because each part of your pipeline can have customized computational requests. So how does the `--dependency` option work? Let's imagine that we have two jobs `Job_A.sbatch` and `Job_B.sbatch`. Additionally, we want to use the output from `Job_A.sbatch` as input for `Job_B.sbatch`. First we would submit `Job_A.sbatch` to the cluster like normal:

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

Now let's imagine we have a case that we will experience shortly. We have a script (Job_C.sbatch) that is dependent on two job IDs (Job_A.sbatch and Job_B.sbatch), as is the case with `GATK` waiting for `Picard` for the normal and tumor samples to finish. In this case, each job ID just needs to be separated by a `:` like:

```
$ sbatch Job_A.sbatch
Submitted batch job 10928
$ sbatch Job_B.sbatch
Submitted batch job 10931
$ sbatch --dependency=afterok:10928:10931 Job_C.sbatch
```

In this case, `Job_C.sbatch` won't queue until `Job_A.sbatch` and `Job_B.sbatch` exit without an error. Furthermore, you can queue a dependent job on a job that is dependent, thus creating an entire pipeline of dependent jobs.

> NOTE: If any part of your pipeline gets an error, then the pipeline will stall and you will need to use `scancel` to cancel your dependent jobs and requeue them.

#### Writing the Wrapper script

```
#!/bin/bash
# This is the wrapper script for the variant calling pipeline. 
# It is designed for a single individual, so if you wanted to run this across multiple individuals you would need to loop this script across all of them.

# Assign variables
FASTQ_DIRECTORY=/home/$USER/variant_calling/raw_fastq/
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
REFERENCE_SEQUENCE_NAME=`basename $REFERENCE_SEQUENCE _genomic.fa`
NORMAL_SAMPLE=syn3_normal
TUMOR_SAMPLE=syn3_tumor
LCR_FILE=/n/groups/hbctraining/variant_calling/reference/LCR-hs38.bed
REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/
SNPEFF_DATABASE=GRCh38.p7.RefSeq

# Intiate bash arrays to hold sample names, Picard Job IDs and coordinate-sorted BAM files
SAMPLE_NAME_ARRAY=()
PICARD_JOB_ID_ARRAY=()
COORDINATE_SORTED_BAM_ARRAY=()

# Submit samples for Alignment and Picard Processing
# For each Read 1 file do:
for SAMPLE in $FASTQ_DIRECTORY*_1.fq.gz; do
  # Parse out the sample name
  SAMPLE_NAME=`basename $SAMPLE _1.fq.gz`
  # Add the sample name to our array holding sample names
  SAMPLE_NAME_ARRAY+=($SAMPLE_NAME)
  # Assign a path and name for the alignments
  SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.sam
  # Submit the bwa sbatch script and save the output to a variable named $BWA_JOB_SUBMISSION
  BWA_JOB_SUBMISSION=$(sbatch -p priority -t 0-04:00:00 -c 8 --mem 16G -o bwa_alignment_${SAMPLE_NAME}_%j.out -e bwa_alignment_${SAMPLE_NAME}_%j.err bwa_alignment_automated.sbatch $REFERENCE_SEQUENCE $SAMPLE $SAM_FILE)
  # Parse out the job ID from outout from the bwa submission
  BWA_JOB_ID=`echo $BWA_JOB_SUBMISSION | cut -d ' ' -f 4`
  # Print to standard output the job that has been submitted
  echo -e "bwa job for sample $SAMPLE_NAME submitted as job ID $BWA_JOB_ID"
  # Submit the picard sbatch script and save the output to a variable named $PICARD_JOB_SUBMISSION
  PICARD_JOB_SUBMISSION=$(sbatch -p priority -t 0-02:00:00 -c 1 --mem 8G -o picard_alignment_processing_${SAMPLE_NAME}_%j.out -e picard_alignment_processing_${SAMPLE_NAME}_%j.err --dependency=afterok:$BWA_JOB_ID picard_alignment_processing_automated.sbatch $SAM_FILE)
  # Parse out the job ID from output from the Picard submission
  PICARD_JOB_ID=`echo $PICARD_JOB_SUBMISSION | cut -d ' ' -f 4`
  # Print to standard output the job that has been submitted
  echo -e "Picard processing job for sample $SAMPLE_NAME submitted as job ID $PICARD_JOB_ID"
  # Add the Picard job ID to an array holding all of the Picard job IDs
  PICARD_JOB_ID_ARRAY+=($PICARD_JOB_ID)
  # Create a variable to hold the the path and output BAM file from the Picard processing sbatch script
  COORDINATE_SORTED_BAM_FILE=`echo -e "${SAM_FILE%sam}coordinate_sorted.bam"`
  # Add this BAM file path to an array 
  COORDINATE_SORTED_BAM_ARRAY+=($COORDINATE_SORTED_BAM_FILE)
done

# Create a string that contains all of the samples preceded by and separated by an "_"
SAMPLE_NAME_STRING=`for i in ${SAMPLE_NAME_ARRAY[@]}; do  echo -ne "_$i"; done`
# Create a string that contains all of the Picard procession job IDs preceded by and separated by a ":"
DEPENDENT_PICARD_JOB_IDS=`for i in ${PICARD_JOB_ID_ARRAY[@]}; do  echo -ne ":$i"; done`

# for loop to discover which position in the $SAMPLE_NAME_ARRAY contains the normal sample
NORMAL_ARRAY_POSITION=`for SAMPLE_POSITION in {0..1}; do if [[ "${SAMPLE_NAME_ARRAY[$SAMPLE_POSITION]}" == "$NORMAL_SAMPLE" ]];then echo $SAMPLE_POSITION; fi; done`
# for loop to discover which position in the $SAMPLE_NAME_ARRAY contains the tumor sample
TUMOR_ARRAY_POSITION=`for SAMPLE_POSITION in {0..1}; do if [[ "${SAMPLE_NAME_ARRAY[$SAMPLE_POSITION]}" == "$TUMOR_SAMPLE" ]];then echo $SAMPLE_POSITION; fi; done`

# Assign variables that will be used by Mutect2 
NORMAL_BAM_MUTECT_INPUT=${COORDINATE_SORTED_BAM_ARRAY[$NORMAL_ARRAY_POSITION]}
TUMOR_BAM_MUTECT_INPUT=${COORDINATE_SORTED_BAM_ARRAY[$TUMOR_ARRAY_POSITION]}
MUTECT2_VCF_OUTPUT=`echo -e "/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/mutect2${SAMPLE_NAME_STRING}_${REFERENCE_SEQUENCE_NAME}-raw.vcf.gz"`

# Submit the Mutect2 sbatch script and save the output to a variable named $MUTECT2_JOB_SUBMISSION 
MUTECT2_JOB_SUBMISSION=$(sbatch -p priority -t 1-00:00:00 -c 1 --mem 16G -o mutect2_variant_calling${SAMPLE_NAME_STRING}_%j.out -e mutect2_variant_calling${SAMPLE_NAME_STRING}_%j.err --dependency=afterok${DEPENDENT_PICARD_JOB_IDS} mutect2_automated.sbatch  $REFERENCE_SEQUENCE $NORMAL_BAM_MUTECT_INPUT $NORMAL_SAMPLE $TUMOR_BAM_MUTECT_INPUT $TUMOR_SAMPLE $MUTECT2_VCF_OUTPUT)

# Parse out the job ID from output from the Mutect2 submission
MUTECT2_JOB_ID=`echo $MUTECT2_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "Mutect2 job submitted as job ID $MUTECT2_JOB_ID"

# Assign a variable that will be used in the Mutect2 filtering step
MUTECT2_VCF_OUTPUT_FILTERED=`echo -e "${MUTECT2_VCF_OUTPUT%raw.vcf.gz}filt.vcf.gz"`

# Submit the variant filtering sbatch script and save the output to a variable named $VARIANT_FILTERING_JOB_SUBMISSION
VARIANT_FILTERING_JOB_SUBMISSION=$(sbatch -p priority -t 0-02:00:00 -c 1 --mem 8G -o variant_filtering${SAMPLE_NAME_STRING}_%j.out -e variant_filtering${SAMPLE_NAME_STRING}_%j.err --dependency=afterok:$MUTECT2_JOB_ID variant_filtering_automated.sbatch $SAMPLE_NAME_STRING $REFERENCE_SEQUENCE $MUTECT2_VCF_OUTPUT $MUTECT2_VCF_OUTPUT_FILTERED)

# Parse out the job ID from output from the variant filtering submission
VARIANT_FILTERING_JOB_ID=`echo $VARIANT_FILTERING_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "Variant filtering job submitted as job ID $VARIANT_FILTERING_JOB_ID"

# Assign variables that will be used in the Snpeff annotation step
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation${SAMPLE_NAME_STRING}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation${SAMPLE_NAME_STRING}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
ANNOTATED_VCF_FILE=`echo -e "${MUTECT2_VCF_OUTPUT_FILTERED%vcf.gz}snpeff.vcf"`

# Submit the variant annotation sbatch script
VARIANT_ANNOTATION_JOB_SUBMISSION=$(sbatch -p priority -t 0-02:00:00 -c 1 --mem 8G -o variant_annotation${SAMPLE_NAME_STRING}_%j.out -e variant_annotation${SAMPLE_NAME_STRING}_%j.err --dependency=afterok:$VARIANT_FILTERING_JOB_ID variant_annotation_automated.sbatch $CSV_STATS $HTML_REPORT $SNPEFF_DATABASE $MUTECT2_VCF_OUTPUT_FILTERED $ANNOTATED_VCF_FILE)

# Parse out the job ID from output from the variant annotation submission
VARIANT_ANNOTATION_JOB_ID=`echo $VARIANT_ANNOTATION_JOB_SUBMISSION | cut -d ' ' -f 4`
# Print to standard output the job that has been submitted
echo -e "Variant annotation job submitted as job ID $VARIANT_ANNOTATION_JOB_ID"
```


