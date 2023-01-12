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

Now we need to replace the following lines of the `sbatch` directives:

```
#SBATCH -o bwa_alignment_normal_%j.out
#SBATCH -e bwa_alignment_normal_%j.err
```

With these lines:

```
#SBATCH -o bwa_alignment_${1}_%j.out
#SBATCH -e bwa_alignment_${1}_%j.err
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
REFERENCE_SEQUENCE=$2
LEFT_READS=$3
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=$4
```

Now our `bwa` submission script should look like:

```
#!/bin/bash
# This script is the automated versin for aligning sequencing reads against a reference genome using bwa

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_${1}_%j.out
#SBATCH -e bwa_alignment_${1}_%j.err

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17

# Assign files to bash variables
REFERENCE_SEQUENCE=$2
LEFT_READS=$3
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=$4

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

Once again, we will change the `sbatch` directives for standard error and standard output from:

```
#SBATCH -o picard_alignment_processing_normal_%j.out
#SBATCH -e picard_alignment_processing_normal_%j.err
```

To:

```
#SBATCH -o picard_alignment_processing_${1}_%j.out
#SBATCH -e picard_alignment_processing_${1}_%j.err
```

Because all of the variables except `$SAM_FILE` are dependent on the value of `$SAM_FILE`, we only need to change this to be a positional parameter.

Change the following line:

```
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.sam
```

To:

```
SAM_FILE=$2
```

Now your `Picard` alignment processing script should look like:

```
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o picard_alignment_processing_${1}_%j.out
#SBATCH -e picard_alignment_processing_${1}_%j.err

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
#SBATCH -o mutect2_variant_calling_%j.out
#SBATCH -e mutect2_variant_calling_%j.err

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

```
#!/bin/bash
# This sbatch script is for variant calling with GATK's MuTect2

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 1-00:00:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH -o mutect2_variant_calling${1}_%j.out
#SBATCH -e mutect2_variant_calling${1}_%j.err

module load gatk/4.1.9.0

REFERENCE_SEQUENCE=$2
REFERENCE_DICTIONARY=`echo ${REFERENCE_SEQUENCE%fa}dict`
NORMAL_BAM_FILE=$3
NORMAL_SAMPLE_NAME=$4
TUMOR_BAM_FILE=$5
TUMOR_SAMPLE_NAME=$6
VCF_OUTPUT_FILE=$7

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
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_filtering_%j.out
#SBATCH -e variant_filtering_%j.err

module load gatk/4.1.9.0
module load snpEff/4.3g

gatk FilterMutectCalls \
--reference /n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
--variant vcf_files/syn3_GRCh38.p7-raw.vcf.gz \
--output vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz

java -jar $SNPEFF/SnpSift.jar intervals -x -i vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz LCR-hs38.bed > vcf_files/syn3_GRCh38.p7-LCR-filt.vcf
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
NORMAL_SAMPLE=syn3_normal
TUMOR_SAMPLE=syn3_tumor

SAMPLE_NAME_ARRAY=()
PICARD_JOB_ID_ARRAY=()
COORDINATE_SORTED_BAM_ARRAY=()
REFERENCE_SEQUENCE_NAME=`basename $REFERENCE_SEQUENCE _genomic.fa`

# Sumbit samples for Alignment and Picard Processing
for SAMPLE in $FASTQ_DIRECTORY*_1.fq.gz; do
  SAMPLE_NAME=`basename $SAMPLE _1.fq.gz`
  SAMPLE_NAME_ARRAY+=($SAMPLE_NAME)
  SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}.sam
  BWA_JOB_SUBMISSION=$(echo -e "sbatch bwa_alignment_automated.sbatch $SAMPLE_NAME $REFERENCE_SEQUENCE $SAMPLE $SAM_FILE")
  BWA_JOB_ID=`echo $BWA_JOB_SUBMISSION | cut -d ' ' -f 4`
  PICARD_JOB_SUBMISSION=$(echo -e "sbatch --dependency=afterok:$BWA_JOB_ID picard_alignment_processing_automated.sbatch $SAMPLE_NAME $SAM_FILE")
  PICARD_JOB_ID=`echo $PICARD_JOB_SUBMISSION | cut -d ' ' -f 4`
  PICARD_JOB_ID_ARRAY+=($PICARD_JOB_ID)
  COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`
  COORDINATE_SORTED_BAM_ARRAY+=($COORDINATE_SORTED_BAM_FILE)
done

SAMPLE_NAME_STRING=`for i in ${SAMPLE_NAME_ARRAY[@]}; do  echo -ne "_$i"; done`
DEPENDENT_PICARD_JOB_IDS=`for i in ${PICARD_JOB_ID_ARRAY[@]}; do  echo -ne ":$i"; done`

NORMAL_ARRAY_POSITION=`for SAMPLE_POSITION in {0..1}; do if [[ "${SAMPLE_NAME_ARRAY[$SAMPLE_POSITION]}" == "$NORMAL_SAMPLE" ]];then echo $SAMPLE_POSITION; fi; done`
TUMOR_ARRAY_POSITION=`for SAMPLE_POSITION in {0..1}; do if [[ "${SAMPLE_NAME_ARRAY[$SAMPLE_POSITION]}" == "$TUMOR_SAMPLE" ]];then echo $SAMPLE_POSITION; fi; done`

MUTECT2_SUBMISSION=$(echo -e "sbatch --dependency=afterok${DEPENDENT_PICARD_JOB_IDS} mutect2_automated.sbatch $SAMPLE_NAME_STRING $REFERENCE_SEQUENCE ${COORDINATE_SORTED_BAM_ARRAY[$NORMAL_ARRAY_POSITION]} ${SAMPLE_NAME_ARRAY[$NORMAL_ARRAY_POSITION]} ${COORDINATE_SORTED_BAM_ARRAY[$TUMOR_ARRAY_POSITION]} ${SAMPLE_NAME_ARRAY[$TUMOR_ARRAY_POSITION]} mutect2${SAMPLE_NAME_STRING}_${REFERENCE_SEQUENCE_NAME}-raw.vcf.gz")



```



---

## Creating many `sbatch` scripts

***YOU DO NOT NEED TO RUN THIS FOR THE WORKSHOP!***

The example we have given above works fine for our limited example. However, you may encounter issues if you have:

- Multiple read groups per sample
- Many samples

In the former case, you will need to manually go through each sample and assign it a unique read group ID, which will be error-prone and cumbersome. In the latter case, it could be time-consuming, tedious and also error-prone. As a result, we have created a script below (alignment_and_processing_script_generator.sh) that will create an `sbatch` script for each read group that you provide it with in a tab-delimited text file. The aim of this is to provide you with a reproducible set of `sbatch` scripts that will carry out your alignment and BAM file processing in parallel. ***ONCE AGAIN, YOU DO NOT NEED TO RUN THIS FOR THE WORKSHOP!***

```
#!/bin/bash                                                                                                                                    
# This is a script that will write your bwa alignment and samtools processing steps                                                            

#######                                                                                                                                        

# EDIT THESE THREE VARIABLES                                                                                                                   
REFERENCE_SEQUENCE=/PATH/TO/YOUR/reference_sequence.fa
METADATA_FILE=/PATH/TO/YOUR/metadata_file.txt
ALIGNMENT_DIRECTORY=/PATH/TO/YOUR/alignments/

#######                                                                                                                                        

# Retrieve column number for various metadata from the metadata file                                                                           
SAMPLE_COLUMN=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "sample") { print i } }}' $METADATA_FILE`
FASTQ_FILE_COLUMN=`awk  'NR==1{for (i=1; i<=NF; i++) { if ($i == "fastq_1_file") { print i } }}' $METADATA_FILE`
RGID_COLUMN=`awk  'NR==1{for (i=1; i<=NF; i++) { if ($i == "RGID") { print i } }}' $METADATA_FILE`
RGPL_COLUMN=`awk  'NR==1{for (i=1; i<=NF; i++) { if ($i == "RGPL") { print i } }}' $METADATA_FILE`
RGPU_COLUMN=`awk  'NR==1{for (i=1; i<=NF; i++) { if ($i == "RGPU") { print i } }}' $METADATA_FILE`
RGSM_COLUMN=`awk  'NR==1{for (i=1; i<=NF; i++) { if ($i == "RGSM") { print i } }}' $METADATA_FILE`

# Set header boolean                                                                                                                           
HEADER=TRUE

# Read in each line of the metadata file                                                                                                       
while read -r line; do
    # If the line is the header, then we are going set out header boolean to FALSE and skip this line                                          
    if [[ $HEADER = "TRUE" ]]; then
        HEADER=FALSE
        continue
    fi
    # Retrieve data from each column in the metadata file and set it equal to variables                                                        
    sample_name=`echo $line | awk -v sample_column=$SAMPLE_COLUMN '{ print $sample_column }'`
    sample_file_1=`echo $line | awk -v fastq_file_column=$FASTQ_FILE_COLUMN '{ print $fastq_file_column }'`
    RGID=`echo $line | awk -v rgid_column=$RGID_COLUMN '{ print $rgid_column }'`
    RGPL=`echo $line | awk -v rgpl_column=$RGPL_COLUMN '{ print $rgpl_column }'`
    RGPU=`echo $line | awk -v rgpu_column=$RGPU_COLUMN '{ print $rgpu_column }'`
    RGSM=`echo $line | awk -v rgsm_column=$RGSM_COLUMN '{ print $rgsm_column }'`
    # Create a variable to hold the sbatch file                                                                                                
    SBATCH_FILE=${sample_name}.alignment_processing.sbatch

    # Write the shebang line and sbatch directives                                                                                             
    echo -e "#!/bin/bash\n" > $SBATCH_FILE
    echo -e "#SBATCH -p short" >> $SBATCH_FILE
    echo -e "#SBATCH -t 0-04:00:00" >> $SBATCH_FILE
    echo -e "#SBATCH -c 8" >> $SBATCH_FILE
    echo -e "#SBATCH --mem 16G" >> $SBATCH_FILE
    echo -e "#SBATCH -o alignment_and_processing_${sample_name}_%j.out" >> $SBATCH_FILE
    echo -e "#SBATCH -e alignment_and_processing_${sample_name}_%j.err\n" >> $SBATCH_FILE

    # Write the modules that need to be loaded                                                                                                 
    echo -e "module load gcc/6.2.0" >> $SBATCH_FILE
    echo -e "module load bwa/0.7.17" >> $SBATCH_FILE
    echo -e "module load samtools/1.15.1\n" >> $SBATCH_FILE
    
    # Assign more variables that we will use that are derived from variables that we've already assigned                                       
    LEFT_READS=$sample_file_1
    RIGHT_READS=`echo ${sample_file_1%1.fq.gz}2.fq.gz`
    SAM_FILE=${ALIGNMENT_DIRECTORY}${sample_name}.sam
    BAM_FILE=`echo ${ALIGNMENT_DIRECTORY}${sample_name}.bam`
    REMOVED_DUPLICATES_BAM_FILE=`echo ${ALIGNMENT_DIRECTORY}${sample_name}.removed_duplicates.bam`

    # Write the bwa alignment command                                                                                                          
    echo -e "bwa mem \\" >> $SBATCH_FILE
    echo -e "\t-M \\" >> $SBATCH_FILE
    echo -e "\t-t 8 \\" >> $SBATCH_FILE
    echo -e "\t-R '@RG\tID:$RGID\tPL:$RGPL\tPU:$RGPU\tSM:$RGSM' \\" >> $SBATCH_FILE
    echo -e "\t$REFERENCE_SEQUENCE \\" >> $SBATCH_FILE
    echo -e "\t$LEFT_READS \\" >> $SBATCH_FILE
    echo -e "\t$RIGHT_READS \\" >> $SBATCH_FILE
    echo -e "\t-o $SAM_FILE\\n" >> $SBATCH_FILE
    
    # Write the samtools sort command                                                                                                          
    echo -e "samtools sort \\" >> $SBATCH_FILE
    echo -e "\t-@ 8 \\" >> $SBATCH_FILE
    echo -e "\t-O bam \\" >> $SBATCH_FILE
    echo -e "\t-o $BAM_FILE \\" >> $SBATCH_FILE
    echo -e "\t$SAM_FILE\\n" >> $SBATCH_FILE

    # Write the samtools markdup command                                                                                                       
    echo -e "samtools markdup \\" >> $SBATCH_FILE
    echo -e "\t-r \\" >> $SBATCH_FILE
    echo -e "\t--write-index \\" >> $SBATCH_FILE
    echo -e "\t-@ 8 \\" >> $SBATCH_FILE
    echo -e "\t$BAM_FILE \\" >> $SBATCH_FILE
    echo -e "\t$REMOVED_DUPLICATES_BAM_FILE\n" >> $SBATCH_FILE
done < $METADATA_FILE
```

The `metadata_file.txt` for our in-class example would look like:

```
sample	fastq_1_file	RGID	RGPL	RGPU	RGSM
syn3-normal /home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz    syn3-normal illumina    syn3-normal syn3-normal
syn3-tumor  /home/$USER/variant_calling/raw_data/syn3_tumor_1.fq.gz syn3-tumor illumina   syn3-tumor  syn3-tumor
```
> NOTE: The order of the columns does not matter.

To run this script you need enter the following command in the directory with the `alignment_and_processing_script_generator.sh` script:

```
sh alignment_and_processing_script_generator.sh
```

We won't exhaustively go through each line of this script as most of it should be readable to you. However, we will briefly discuss a  few lines:

1. There are several lines like this one: 

```
SAMPLE_COLUMN=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "sample") { print i } }}' $METADATA_FILE`
```

This is an `awk` command that is looping through every element in the header line and finding the one that matches a particular pattern (in this case, "sample"). When it finds that match it saves the column number to a variable (in this case, "SAMPLE_COLUMN"). The reason we do this way is that it allows there be more flexible in the order of the columns. 

2. The `while` loop 
 
```
while read -r line; do
...
done < $METADATA_FILE
```

This is called a `while` loop and it is a great way to read in lines from a file. Essentially, we are reading in data line-by-line from the file attached to the variable `METADATA_FILE` *while* there are still lines to read in. Each line is assigned to the variable `line`. Once there are no more lines to read in, the loop stops.

3. The header boolean and the `if` statement

```
HEADER=TRUE
...
if [[ $HEADER = "TRUE" ]]; then
    HEADER=FALSE
    continue
fi
```

We are setting a variable `HEADER` equal to `TRUE` outside of the loop, so that the first time it goes through the loop this conditional statement will return true and it will change the `HEADER` variable to `FALSE` and ignore the rest of the loop for this line of input.

4. Extracting information from a column

```
sample_name=`echo $line | awk -v sample_column=$SAMPLE_COLUMN '{ print $sample_column }'`
```

Since, we know which column holds which information (see Explanation 1 above), we are going to provide that column number to awk and ask it to save the content of that column to a variable (in this case, `sample_name`).
