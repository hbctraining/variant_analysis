# Sequence Alignment Theory

## Learning Objectives

- Enumerate difficulties with alignment
- Create an `sbatch` script to align reads
- Implement `sed` to search-and-replace text in a text file

<p align="center">
<img src="../img/Variant_calling_workflow.png" width="300">
</p>

## Understanding Alignment

Sequence Read Alignment may appear to be a simple task at first glance. One has set of reads where each read represents a small section of the genome and therefore, it would simply appear that one just needs to evaluate each read and assign it to a location in the genome matching that sequence. However, this task is complicated by a multitude of factors including:

- Errors in the sequencing reads
- Errors in the reference genome
- Repetitive areas of the genome
- Variants in the sample

<p align="center">
<img src="../img/Alignment_errors.png" width="800">
</p>

However, likely the largest and most difficult part of this task is creating an alignment algorithm that can account for these factors and still provide alignment for ***millions*** of reads within a reasonable time frame. There are a multitude of alignment tools availible today, including [bwa](https://bio-bwa.sourceforge.net), [HISAT2](http://daehwankimlab.github.io/hisat2/) and [STAR](https://github.com/alexdobin/STAR). Many of them has strengths and weakness and are more appropriate for a certain type of analysis. In this course, we will be using `bwa`. 

## `bwa`

Many modern alignment tools rely on the Burrows-Wheeler Transform as part of their alignment and `bwa` does as well. `bwa` runs in several parts:

1. Create an index of the reference sequence (done once)
    1. Create a Suffix Array Index
    2. Generate a Burrows-Wheeler Transform from the suffix array
2. Query reads against the index to get their alignments

> The Burrows-Wheeler Transform and Suffix Arrays are outside of the scope of this course, but a great place to get information on these methods are on [Dr. Ben Langmead's YouTube Channel](https://www.youtube.com/user/BenLangmead/videos). Dr. Langmead developed the algorithm for a similar alignment tool, Bowtie, and he makes many useful videos explaining how the components of these alignment tools work. 

### Creating a `bwa` index

While we will be using an index that has already been made for us, if you need to create and index for a reference sequence using `bwa`, the steps for this are laid out below.

**1) Navigate to your reference sequence directory**

```
cd ~/path/to/reference/sequence/directory/
```

**2) Create a `bwa` index**

```
bwa index reference_sequence.fasta
```

This process may take up to 30+ minutes to run depending on the reference sequence size. The output of this command will produce five files:

- `reference_sequence.fasta.sa` This is the binary version of the Suffix Array Index
- `reference_sequence.fasta.bwt` This is the binary version of the Burrows-Wheeler Transform of the reference sequence
- `reference_sequence.fasta.pac` A special binary compression of the reference sequence
- `reference_sequence.fasta.ann` Notations regarding the reference sequence
- `reference_sequence.fasta.amb` Notations regarding base ambiguities (mostly Ns, but also other base ambiguities) in the reference sequence

### Setting up our `sbatch` Submission Script for alignment

Now that we have indexed the reference sequence, we can align sequence reads to our indexed reference sequence. To align reads to the reference sequence, we will need to create and `sbatch` submission script in `vim`.

```
cd ~/variant_calling/scripts/

vim bwa_alignment_normal.sbatch
```

Now that we have opened up `vim` we need to enter `insert-mode` by pressing `i`. Once in `insert-mode`, we can copy and paste the following shebang line and `sbatch` directives:

```
#!/bin/bash
# This script is for aligning sequencing reads against a reference genome using bwa

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_normal_%j.out
#SBATCH -e bwa_alignment_normal_%j.err
```
Next, we need to add the modules that we will be using in this command:

```
# Load modules
module load gcc/6.2.0
module load bwa/0.7.17
```

> Remember that on O2, many of the common tools were compiled using `GCC` version 6.2.0, so to be able to access them, we first need to load the `GCC` module.

> Also, you can see that we also loaded the `samtools` version 1.15.1 module. That is not required for alignment, however, this script will have a few steps after alignment that will require `samtools`. It is common to load all modules at the top of a `sbtach` submission script. The only exception to this practice happen in the rare case that they are dependency conflicts between modules. In this case, you may see people loading/unloading/swapping modules within their `sbatch` scripts.

Now that we have the module load command for `bwa` in our SBATCH script, we are going to declare some bash variables that we are going to use:

```
# Assign files to bash variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam
```

Some of these variable assignment are straightforward and are simply assigning paths to known files to `bash` variables. However, `RIGHT_READS` uses a little `bash` trick in it in order to swap the last parts of their name. Instead, we could have simply written:

```
RIGHT_READS=/home/$USER/variant_calling/syn3_normal_2.fq.gz
```

However, we chose to do it this way for two reasons:
1. Each time we use this script moving forward, we will never need to edit the `RIGHT_READS` variable
2. Reduces the chance for typos
3. Can help keep filenames nonmenclature consistent across files


#### Brief `bash` Text Manipulation 

Let's briefly explore how this works and discuss why we choose this way of doing it. In order to discuss this, let's create a toy example of a path:

```
TOY_EXAMPLE_PATH=~/variant_calling/sam_alignments/filename.sam
```

`bash` has some clever text manipulation tools that are going to help us. Let's introduce the `%` tool for text manipulation. `%` is placed within the `{}` of a variable and tells `bash` to remove the shortest match from the end that contains the text that follows the `%`.

```
echo ${TOY_EXAMPLE_PATH%sam}
```

We can see that we have maintained the path and stripped the `sam` extension and now all we need to do is add the new `.bam` extension. To do this, just add it to the end of the variable outside of the `{}`

```
echo ${TOY_EXAMPLE_PATH%sam}bam
```

A brief overview of some `bash` text manipulation shortcuts are in the table below:

| Shortcut | Effect |
|------|------|
| % | Remove shortest match from the end of the string |
| %% | Remove longest match from the end of the string|
| # | Remove the shortest match from the beginning of the string|
| ## | Remove the longest match from the beginning of the string|

> NOTE: In the above example `%` and `%%` would give the same result, however, their differences become more clear if we used an `*` character. In the dropdown below, we provide a few more examples for those who are curious for examples of each of these text manipulation tools.

> NOTE: The `#` and `##` within `{}` for text manipulation is one of the rare times in `bash` that uses `#` for a purpose other than commenting code!

<details>
  <summary><b>Click here for more text manipulation examples in <code>bash</code></b></summary>
  <br><code><b>%</b></code> Removes the shortest match of "sam" from the end of the string
  <pre>echo ${TOY_EXAMPLE_PATH%sam*}</pre>
  Returning:<br>
  <pre>/home/wig051/variant_calling/sam_alignments/filename.</pre>
  <code><b>%%</b></code> Removes the longest match of "sam" from the end of the string<br>
  <pre>echo ${TOY_EXAMPLE_PATH%%sam*}</pre>
  Returning:<br>
  <pre>/home/wig051/variant_calling/</pre>
  <code><b>#</b></code> Removes the shortest match of "sam" from the start of the string<br>
  <pre>echo ${TOY_EXAMPLE_PATH#*sam}</pre>
  Returning:<br>
  <pre>_alignments/filename.sam</pre>
  <code><b>##</b></code> Removes the longest match of "sam" from the start of the string<br>
  <pre>echo ${TOY_EXAMPLE_PATH##*sam}</pre>
  Which deletes the entire string!
</details>
  
# Short Read Alignment

Now that we have added the bash variables to our `sbatch` submission script, we can now add our command for running `bwa`:

```
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

Let's breakdown this `bwa` command.

- `mem` This is the specific tools we want to use within `bwa` and this is one of the tools that `bwa` provides for alignment.

- `-M` This will mark shorter split hits as secondary

- `-t 8` We are going to take advantage of multithreading. In order to do this, we need to specifiy the number of ***t***hreads. We are going to use 8.

- `-R '@RG\tID:syn3-normal\tPL:illumina\tPU:syn3-normal\tSM:syn3-normal'` This adds what is called Read Group information. Some software packages, such as `GATK`, require Read Groups, while others are agnostic towards them. However, they can provide important metadata about your reads and thus it is considered best practice to include Read Group Information at the alignment step. This Read Group information consists of several fields separated by tab-characters (\t), including:

    - **ID**: This is the identification for a given batch of reads. This ***MUST*** be unique to your experiment. 

    - **PL**: This is the platform that the sequencing was run on. For aligning Illumina reads, you should use "illumina" here.

    - **PU**: This is the platform unit and it is ideally supposed to hold `<FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>`, where `<FLOWCELL_BARCODE>` is the barcode of the flowcell, `<LANE>` is the lane the data was run on and `<SAMPLE_BARCODE>` is supposed to be a library/sample specific identifer. In some software packages **PU** can take precedence over the **ID** field. If you don't happen to have the `<FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>`, just make this field something useful that will help identify the sample. In this case, we didn't have that information so we are re-using the **ID** field here. 

    - **SM**: This is to mark which *sample* your reads are coming from. Note, this does not need to be unique like the **ID** field, since you may have multiple Read Group **ID**s coming from a single sample. For example, you may be merging alignment (BAM/SAM) files that originated from the same individual, but we sequenced on different lanes or machines. In this case, those would have different Read Group **ID**s, but the same **SM** value.

    - More information about Read Groups and some fields we didn't discuss can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).
  
- `$REFERENCE_SEQUENCE` This calls the variable that holds the path to our reference sequence (`/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa`).
  
- `$LEFT_READS` and `$RIGHT_READS` These are variables with paths to the FASTQ sequence read files that we wish to align (`~/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_1.fq.gz` and `~/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_2.fq.gz`). In this case you can see that the reads are paired-end reads because 1) we have provided two sets of reads and 2) they have `_1`/`_2` right before the `.fq`/`.fastq` extension. This is a very common annotation for pair-end reads. Alternatively, sometimes they might also be annotated as `_R1`/`_R2`, but once again, they will almost always be placed right before the `.fq`/`.fastq` extension. We can also, note that these reads are currently compressed with gzip as annotated by the `.gz` extension. This is not a problem as `bwa` can read both gzip compressed and uncompressed FASTQ files.

- `-o $SAM_FILE` This calls a variable that holds the path to your output file (`~/variant_calling/alignments/normal_GRCh38.p7.sam`) where your alignments written. Notice this is a SAM file (an uncompressed alignment file).

Now you have written your command to run `bwa` you are ready to run alignment. However, there are a few steps that we are going to add to the script so that they run immediately after the alignemnt finishes.

### Submitting `sbatch` bwa script

Now your `sbatch` script for `bwa` should look like this:

```
#!/bin/bash
# This script is for aligning sequencing reads against a reference genome using bwa

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_normal_%j.out
#SBATCH -e bwa_alignment_normal_%j.err

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17

# Assign files to bash variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam

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

If so, go ahead and submit the our `sbatch` script to the cluster:

```
sbatch bwa_alignment_normal.sbatch
```

## Alignment file processing with `samtools` and `Picard`

The processing of the alignment files (SAM/BAM files) can be done either with [`samtools`](https://github.com/samtools/samtools) or [`Picard`](https://broadinstitute.github.io/picard/) and they are for the most part interchangable. The arguments for either are below:

`samtools`
- The software we will be using to call variants (GATK) can be picky about how the files are formatted and so you need to be careful about the formatting as to not produce an error in GATK
- More user-friendly
- Multi-threaded

`Picard`
- Maintained by the Broad Institute, so it has excellent integration with GATK (also maintained by the Broad Institute)
- Written in Java and the error messages are not as easily interpretable 
- Single-threaded

If implemented correctly, they largely provide the same output. In this workshop, we will be using `Picard`, but the `samtools` code will be provided in a dropdown for each section if you would like to know how to do the step in `samtools`.

### Pipeline for processing alignment file with `Picard`

Before we start processing our alignment SAM file provided by `bwa`, let's briefly discuss the steps that we will be doing in this pipeline. Several goals need to be accomplished:

1) **Compress SAM file to BAM file.** The output of `bwa` is a SAM file and it is human readbale. However, it is quite large and we need to compress it to a binary version (BAM) which is much smaller.
2) **Query-sort our alignment file.** Alignment file are initally ordered by the order of the reads in the FASTQ file, which is not particularly useful. `Picard` can more exhaustively look for duplicates if the file is sorted by read-name (query-sorted). We will discuss **query**-sorted and **coordinate**-sorted alignment files soon.
3) **Mark and Remove Duplicates.** Duplicates can introduce bias into our analysis so it is considered best practice to remove them prior to variant calling.
4) **Coordinate-sort our alignment file.** Most downstream software packages require that alignment files be **coordinate**-sorted, so we will need to re-sort our alignment file by **coordinates** now that we have remove the duplicates.
5) **Index the alignment file.** Like the index for a book, indicies for alignment files help direct downstream software packages to where to they can find specific reads. Many software packages require the alignment file that you are analyzing to have an index file, usually with the same name as you alignment file, with the additional `.bai` (BAM-index) extension. Both `Picard` and `samtools` have a way of integrating this indexing as part of their sorting protocols and that is what we will be using. However, both packages have commands for indexing a BAM file independent of their sorting protocol ([BuildBamIndex](https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-) for `Picard` and [index](http://www.htslib.org/doc/samtools-index.html) for `samtools`).

Below is a flow chart of the `Picard` pipeline that we will be using:

<p align="center">
<img src="../img/Picard_pipeline.png" width="800">
</p>

### Creating our `sbatch` script

Let's go ahead and start making a new `sbatch` within `vim`:

```
vim picard_alignment_processing_normal.sbatch
```

Start the `sbatch` script with our shebang line, description of the script and our `sbatch` directives. 

```
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o picard_alignment_processing_normal_%j.out
#SBATCH -e picard_alignment_processing_normal_%j.err
```

Next, let's define some variables that we will be using:

```
# Assign file paths to variables 
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}remove_duplicates.bam`
METRICS_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}remove_duplicates_metrics.txt`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
```

Load the `Picard` module: 

```
module load picard/2.8.0
```

**Note: `Picard` is one of the pieces of software that does NOT require gcc/6.2.0 to also be loaded** 

### Sorting and Removing Duplicates

In order to appropriately flag and remove duplicates, we first need to ***query*** sort our SAM file. Oftentimes, when people discuss sorted BAM/SAM files, they are refering to **coordinate**-sorted BAM/SAM files. 

- **Query**-sorted BAM/SAM files are sorted based upon their read names and order lexiographically
- **Coordinate**-sorted BAM/SAM files are sorted by their sequence name (chromosome/linkage group/scaffold) and position 

<p align="center">
<img src="../img/SAM_sorting.png" width="800">
</p>

`Picard` can mark and remove duplicates in either coordinate-sorted or query-sorted BAM/SAM files, however, if the alignments are query-sorted it can test secondary alignments for duplicates. A brief discussion of this nuance is discussed in the [`MarkDuplicates` manual of `Picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-). As a result, we will first **query**-sort our SAM file and convert it to a BAM file:

#### Query-sort the Alignment File

```
java -jar $PICARD/picard-2.8.0.jar SortSam \
INPUT=$SAM_FILE \
OUTPUT=$QUERY_SORTED_BAM_FILE \
SORT_ORDER=queryname
```

The components of this command are:

- `java -jar $PICARD/picard-2.8.0.jar SortSam ` Calls `Picard`'s `SortSam` software package

- `INPUT=$SAM_FILE` This is where we provide the SAM input file

- `OUTPUT=$QUERY_SORTED_BAM_FILE` This is the BAM output file. Because the extension is `.bam` rather than `.sam`, `Picard` will recognize this and create the output as a BAM file rather than the SAM inoput we have provided it.

- `SORT_ORDER=queryname` The method with which we would like the file to be sorted. The options here are either `queryname` or `coordinate`.

#### Mark and Remove Duplicates

Now we will add the command to our script that allows us to mark and remove duplicates:

```
java -jar $PICARD/picard-2.8.0.jar MarkDuplicates \
INPUT=$QUERY_SORTED_BAM_FILE \
OUTPUT=$REMOVE_DUPLICATES_BAM_FILE \
METRICS_FILE=$METRICS_FILE \
REMOVE_DUPLICATES=true
```

The componetns of this command are:

- `java -jar $PICARD/picard-2.8.0.jar MarkDuplicates` Calls `Picard`'s `MarkDuplicates` program

- `INPUT=$QUERY_SORTED_BAM_FILE` Uses our query-sorted BAM file as input

- `OUTPUT=$REMOVE_DUPLICATES_BAM_FILE` Write the output to a BAM file

- `METRICS_FILE=$METRICS_FILE` Creates a metrics file (required by `Picard MarkDuplicates`)

- `REMOVE_DUPLICATES=true` Not only are we going to mark/flag our duplicates, we can also remove them. By setting the `REMOVE_DUPLICATES` parameter equal to `true` to can remove the duplicates.

#### Coordinate-sort the Alignment File

For most downstream processes, coordinate-sorted alignment files are required. As a result, we will need to change our alignemnt file from being **query**-sorted to being **coordinate**-sorted and we will once again use the `SortSam` command within `Picard` to accomplish this.

```
java -jar $PICARD/picard-2.8.0.jar SortSam \
INPUT=$REMOVE_DUPLICATES_BAM_FILE \
OUTPUT=COORDINATE_SORTED_BAM_FILE \
SORT_ORDER=coordinate \
CREATE_INDEX=true
```

The components of this command are:

- `java -jar $PICARD/picard-2.8.0.jar SortSam` Calls `Picard`'s `SortSam` program

- `INPUT=$REMOVE_DUPLICATES_BAM_FILE` Our BAM file once we have removed our duplicates. **NOTE: Even though the software is called `SortSam`, it can use BAM or SAM files as input and also BAM or SAM files as output.**

- `OUTPUT=COORDINATE_SORTED_BAM_FILE` Our BAM output file sorted by coordinates.

- `CREATE_INDEX=true` Since this BAM file will be the final BAM file that we make and will use for downstream analyses, we will need to create an index for it. Setting the `CREATE_INDEX` equal to `true` will create an index of the final BAM output. The task can also be accomplished by using the `BuildBamIndex` command within `Picard`, but this `CREATE_INDEX` functionality is built into many `Picard` function, so you can often use it at the last stage of processing your BAM file to save having to run `BuildBamIndex` after.

---

Your final `sbatch` script for `Picard` should look like:

```
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o picard_alignment_processing_normal_%j.out
#SBATCH -e picard_alignment_processing_normal_%j.err

# Assign file paths to variables 
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}remove_duplicates.bam`
METRICS_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}remove_duplicates_metrics.txt`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`

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

<details>
  <summary>Click here for alignment file processing using <code>Samtools</code></summary>
<br><ol><li><details>
    <summary>Click here for setting up a <code>sbatch</code> script BAM/SAM Processing for the <code>Samtools</code> pipeline</summary>
<br>Below is the pipeline and explanation for how you would carry out the similar SAM/BAM processing steps within <code>Samtools</code>.<br>
<h2>Setting up <code>sbatch</code> Script</h2>
First we are going to need to set-up our <code>sbatch</code> submission script with our shebang line, <code>sbatch</code> directives, modules to load and file variables.
<pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Samtools<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o samtools_processing_normal_%j.out
#SBATCH -e samtools_processing_normal_%j.err<br>
# Load modules
module load gcc/6.2.0
module load samtools/1.15.1<br>
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
FINAL_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}final.bam`<br>
</pre>
</details></li>

<li><details>
    <summary>Click here for <b>Query</b>-sorting a SAM file and converting it to BAM for the <code>Samtools</code> pipeline</summary>
Similarly to <code>Picard</code>, we are going to need to initally <b>query</b>-sort our alignment. We are also going to be converting the SAM file into a BAM file at this step. Also similarly to <code>Picard</code>, we don't need to specify that our input or output files are BAM or SAM files. <code>Samtools</code> will use the extensions you provide it in your file names as guidance for whether you are providing it a BAM/SAM and whether you want the output to be a BAM/SAM file. Below is the code we will use to <b>query</b>-sort our SAM file and convert it into a BAM file:<br>
    
<pre>
# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
-@ 8 \
-n \
-o $QUERY_SORTED_BAM_FILE \
$SAM_FILE
</pre>
    
The components of this line of code are:
    
<ul><li><code>samtools sort</code> This calls the sort function within <code>samtools</code>.</li>

<li><code>-@ 8</code> This tells <code>samtools</code> to use 8 threads when it multithreads this task. Since we requested 8 cores for this <code>sbatch</code> submission, let's go ahead and use them all.</li>

<li><code>-n</code> This argument tells <code>samtools sort</code> to sort by read name as opposed the the default sorting which is done by coordinate.</li>

<li><code>-O bam</code> This is declaring the output format of `bam`.</li>

<li><code>-o $QUERY_SORTED_BAM_FILE</code> This is a <code>bash</code> variable that holds the path to the output file of the <code>samtools sort</code> command.</li>

<li><code>$SAM_FILE</code> This is a <code>bash</code> variable holding the path to the input SAM file.</li></ul>
</details></li>

<li><details>    
<summary>Click here for fixing mate information for the <code>Samtools</code> pipeline</summary>
Next, we are going to add more mate-pair information to the alignments including the insert size and mate pair coordinates. It is important to note with this command that <code>samtools</code> relies on positional parameters for assigning the the input and output BAM files. In this case the input BAM file (<code>$QUERY_SORTED_BAM_FILE</code>) needs to come before the output file (<code>$FIXMATE_BAM_FILE</code>):
    
<pre>
# Score mates
samtools fixmate \
-m \
$QUERY_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE   
</pre>

The parts of this command are:

<ul><li><code>samtools fixmate</code> This calls the <code>fixmate</code> command in <code>samtools</code></li>

<li><code>-m</code> This will add the mate score tag that will be critically important later for <code>samtools markdup</code></li>

<li><code>$QUERY_SORTED_BAM_FILE</code> Bash variable that holds the path to the input file</li>

<li><code>$FIXMATE_BAM_FILE</code> Bash variable that holds the path to the input file</li></ul>
</details></li>

<li><details>
<summary>Click here for <b>coordinate</b>-sorting a BAM file for the <code>Samtools</code> pipeline</summary>
    
Now that we have added the <code>fixmate</code> information, we need to <b>coordinate</b>-sort the BAM file. We can do that by: 

<pre>
# Sort BAM file by coordinate   
samtools sort \
-@ 8 \
-o $COORDINATE_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE
</pre>

We have gone through all of the these paramters already in the previous <code>samtools sort</code> command. The only difference in this command is that we are not using the <code>-n</code> option, which tell <code>samtools</code> to sort by read name. Now, we are sorting by coordinates, the default setting.
</details></li>
    
<li><details>
<summary>Click here for marking and removing duplicates for the <code>Samtools</code> pipeline</summary>

<pre>
# Mark and remove duplicates and then index the output file
samtools markdup \
-r \
--write-index \
-@ 8 \
$COORDINATE_SORTED_BAM_FILE \
${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai
</pre> 

The components of this command are:    
    
<ul><li><code>samtools markdup</code> calls the mark duplicates software in <code>samtools</code></li>
    
<li><code>-r</code> removes the duplicate reads</li>
    
<li><code>--write-index</code> writes an index file (see next section for more information on BAM index files) of the output</li>
    
<li><code>-@ 8</code> sets that we will be using 8 threads</li>
    
<li><code>$BAM_FILE</code> this is out BAM input file</li>
    
<li><code>${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai</code>This has two parts:
<ol><li>The first part (<code>${REMOVED_DUPLICATES_BAM_FILE}</code>) is our BAM output file with the duplicates removed from it</li>
<li>The second part (<code>##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai</code>) is a shortcut to creating a <code>.bai</code> index of the BAM file. If we use the <code>--write-index</code> option without this second part, it will create a <code>.csi</code> index file. <code>.bai</code> index files are a specific type of <code>.csi</code> files, so we need to specify it with the second part of this command to ensure that a <code>.bai</code> index file is created rather than a <code>.csi</code> index file.</li></ol></li></ul>
</details></li>

<li><details>
<summary>Click here for the final <code>sbatch</code> script to do the BAM/SAM processing for the <code>Samtools</code> pipeline</summary>

The final script should look like:
    
<pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Samtools<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o samtools_processing_normal_%j.out
#SBATCH -e samtools_processing_normal_%j.err<br>
# Load modules
module load gcc/6.2.0
module load samtools/1.15.1<br>
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
FINAL_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}final.bam`<br>
# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
-@ 8 \
-n \
-o $QUERY_SORTED_BAM_FILE \
$SAM_FILE<br>
# Score mates
samtools fixmate \
-m \
$QUERY_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE<br>
# Sort BAM file by coordinate   
samtools sort \
-@ 8 \
-o $COORDINATE_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE<br>
# Mark and remove duplicates and then index the output file
samtools markdup \
-r \
--write-index \
-@ 8 \
$COORDINATE_SORTED_BAM_FILE \
${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai<br>
</pre>
</details></li>
</details>


## Query-name Sorted Reads and SAM to BAM Conversion

Now that we have aligned our reads to the reference, we would see if we opened up the SAM file that the alignments are in the order they were processed by `bwa` and not in any particular order that would be useful for downstream analyses. So, we are going to ***sort them into order by read name*** for the `fixmates` tool in `samtools` downstream. It should be noted that we are going to sort our BAM file by coordinate later in the processing and when people discuss a "sorted SAM/BAM" file they are usually referring to a BAM file that is coordinate sorted.

Additionally, while SAM files are nice due to their human readability, they are typically quite large files and it is not an efficient use of space on the cluster. Fortunately, there is a binary compression version of SAM called BAM. While we sort the reads, we are going to use a very common tool, `samtools`, to convert our SAM files to BAM files. 


## Marking and Removing Duplicates

An important step in processing a BAM file is to mark and remove PCR duplicates. These PCR duplicates can introduce artifacts because regions that have preferential PCR amplification could be over-represented. These reads are flagged by having identical mapping locations in the BAM file. Importantly, it is impossible to distinguish between PCR duplicates and identical fragments. However, one can reduce the latter by doing paired-end sequencing and providing appropriate amounts of input material. We can mark and remove duplicates in `samtools` as well:



## Indexing the `BAM` File

Similar to the index of a book. Many software packages want an index of your BAM file in order to facilitate fast look-ups of a BAM file. While not all software packages that use a BAM file will require this, many will and thus it is a good practice to just index your BAM file. In our previous line of `samtools` command, we provided it with the `--write-index` option, so it automatically created an index for us after marking and removing duplicates. This option exists in the `markdup` command because this is often the last step of BAM file processing that people carry out, so it makes sense to offer the ability to index the BAM file at this point.

<h2>BAM-Indexing within <code>Samtools</code></h2>
    
In the previous step, we have already indexed our BAM file. But if for some reason we needed to index a BAM file, the command to index a BAM file with <code>samtools</code> would be:

<pre>
#### SKIP THIS STEP
# Index the BAM file
samtools index \
$BAM_FILE
</pre>

<code>samtools index</code> Calls the <code>index</code> function within <code>samtools</code>

<code>$BAM_FILE</code> This is a `bash` variable that holds the path to the BAM file that we want to index.

We don't need to provide an output file for <code>samtools index</code>, by default it will generate a new file using the same path and filename as the BAM file, but add `bai` as the extension to denote that it is a BAM-index file.


The tool we will be using for variant calling is called `GATK` and it was developed and maintained by the Broad Institute. The Broad Institute also maintains a tool that does many of the functions that `samtools` does and it is called `Picard`. In the dropdown below, we show the command you can use for sorting a SAM file in `Picard`. 

<details>
  <summary>SAM to BAM conversion in <code>Picard</code></summary>
  
  One benefit of using Picard for this task is that it will convert the SAM file to BAM and do the indexing of the BAM file in the same step. However, Picard doesn't support multithreading, so it is slower than `samtools`. Additionally, given that we requested 8 cores, we would be wasting 7 cores while running `Picard`, so the best practice be to put the below `Picard` command in a separate `sbatch` script. 
  
  <pre>
  module load picard/2.8.0
  
  java -jar $PICARD/picard-2.8.0.jar SortSam INPUT=$SAM_FILE OUTPUT=$BAM_FILE SORT_ORDER=coordinate CREATE_INDEX=true
  </pre>
</details>

---

The final script for your normal sample should look like:

```
#!/bin/bash

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_samtools_sorting_index_normal_%j.out
#SBATCH -e bwa_alignment_samtools_sorting_index_normal_%j.err

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1 

# Assign files to bash variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
FINAL_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}final.bam`

# Align reads with bwa
bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-normal\tPL:illumina\tPU:syn3-normal\tSM:syn3-normal' \
$REFERENCE_SEQUENCE \
$LEFT_READS \
$RIGHT_READS \
-o $SAM_FILE

# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
-@ 8 \
-n \
-o $QUERY_SORTED_BAM_FILE \
$SAM_FILE

# Score mates
samtools fixmate \
-m \
$QUERY_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE  

# Sort BAM file by coordinate   
samtools sort \
-@ 8 \
-o $COORDINATE_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE

# Mark and remove duplicates and then index the output file
samtools markdup \
-r \
--write-index \
-@ 8 \
$COORDINATE_SORTED_BAM_FILE \
${FINAL_BAM_FILE}##idx##${FINAL_BAM_FILE}.bai
```

## Creating Tumor `sbatch` script

Now that we have created the `sbatch` script for our normal samples, we need to repeat the process for our tumor samples. All of the parameters will stay the same, we just need to edit the SBATCH error file, SBATCH output file, LEFT_READS variable, RIGHT_READS variable, SAM_FILE variable and the read group information within the `bwa` command. You could very well do this by hand and it would be just fine. However, to cut down on typos we are going to use `sed`. `sed` is a powerful tool within `bash` and [has a wide variety of applications](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/sed.html). However, one of the most common uses for `sed` is as a "find-and-replace" tool. The syntax for this type of task is:

```
sed 's/pattern/replacement/g' file.txt 
```

The `s` before `/pattern/replacement/` is telling `sed` that we are going to use its **substittion** function and the `g` after `/pattern/replacement/` is telling `sed` that we want to apply that change **globally**, or every instance in the file. `pattern` represents the pattern that we are looking for and `replacement` is what we wish to replace the `pattern` with. Lastly, we need to provide `sed` some text source to apply this "find-and-replace" function, so we have provided it with `file.txt`, but you can also pipe in a string or file to apply this function to. 

In our case, we are hoping to replace each instance of "normal" with "tumor". Therefore, we could call `sed` to do this using:

```
sed 's/normal/tumor/g' bwa_alignment_samtools_processing_normal.sbatch
```

We can see that all instances of "normal" have been replaced with "tumor". Now we would like to redirect this output to a file called `bwa_alignment_samtools_sorting_index_tumor.sbatch` rather than standard output, se we need to add redirection to the end of out command:

```
sed 's/normal/tumor/g' bwa_alignment_samtools_processing_normal.sbatch >  bwa_alignment_samtools_processing_tumor.sbatch 
```

If we look at the output with:

```
cat bwa_alignment_samtools_processing_tumor.sbatch 
```

Then it should look like this:

```
#!/bin/bash

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_samtools_sorting_index_tumor_%j.out
#SBATCH -e bwa_alignment_samtools_sorting_index_tumor_%j.err

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1 

# Assign files to bash variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_tumor_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/tumor_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
FINAL_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}final.bam`

# Align reads with bwa
bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-tumor\tPL:illumina\tPU:syn3-tumor\tSM:syn3-tumor' \
$REFERENCE_SEQUENCE \
$LEFT_READS \
$RIGHT_READS \
-o $SAM_FILE

# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
-@ 8 \
-n \
-o $QUERY_SORTED_BAM_FILE \
$SAM_FILE

# Score mates
samtools fixmate \
-m \
$QUERY_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE  

# Sort BAM file by coordinate   
samtools sort \
-@ 8 \
-o $COORDINATE_SORTED_BAM_FILE \
$FIXMATE_BAM_FILE

# Mark and remove duplicates and then index the output file
samtools markdup \
-r \
--write-index \
-@ 8 \
$COORDINATE_SORTED_BAM_FILE \
${FINAL_BAM_FILE}##idx##${FINAL_BAM_FILE}.bai
```

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

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
