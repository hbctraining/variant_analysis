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

### Setting up our `sbatch` Submission Script

Now that we have indexed the reference sequence, we can align sequence reads to our indexed reference sequence. To align reads to the reference sequence, we will need to create and `sbatch` submission script in `vim`.

```
cd ~/variant_calling/scripts/

vim bwa_alignment_samtools_sorting_index_normal.sbatch
```

Now that we have opened up `vim` we need to enter `insert-mode` by pressing `i`. Once in `insert-mode`, we can copy and paste the following shebang line and `sbatch` directives:

```
#!/bin/bash

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_samtools_sorting_index_normal_%j.out
#SBATCH -e bwa_alignment_samtools_sorting_index_normal_%j.err
```
Next, we need to add the modules that we will be using in this command:

```
# Load modules
module load gcc/6.2.0
module load bwa/0.7.17
module load samtools/1.15.1 
```

> Remember that on O2, many of the common tools were compiled using `GCC` version 6.2.0, so to be able to access them, we first need to load the `GCC` module.

> Also, you can see that we also loaded the `samtools` version 1.15.1 module. That is not required for alignment, however, this script will have a few steps after alignment that will require `samtools`. It is common to load all modules at the top of a `sbtach` submission script. The only exception to this practice happen in the rare case that they are dependency conflicts between modules. In this case, you may see people loading/unloading/swapping modules within their `sbatch` scripts.

Now that we have the module load command for `bwa` in our SBATCH script, we are going to declare some bash variables that we are going to use:

```
# Assign files to bash variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/n/data1/cores/bcbio/gammerdinger/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/data1/cores/bcbio/gammerdinger/variant_calling/alignments/normal_GRCh38.p7.sam
BAM_FILE=`echo ${SAM_FILE%sam}bam`
REMOVED_DUPLICATES_BAM_FILE=`echo ${BAM_FILE%bam}removed_duplicates.bam`
```

Some of these variable assignment are straightforward and are simply assigning paths to known files to `bash` variables. However, `RIGHT_READS`, `BAM_FILE` and `REMOVED_DUPLICATES_BAM_FILE` all use a little `bash` trick in it in order to swap the last parts of their name. Instead, we could have simply written:

```
RIGHT_READS=/n/data1/cores/bcbio/gammerdinger/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_2.fq.gz
BAM_FILE=~/variant_calling/alignments/normal_GRCh38.p7.bam
REMOVED_DUPLICATES_BAM_FILE=~/variant_calling/alignments/normal_GRCh38.p7.removed_duplicates.bam
```

However, we chose to do it this way for two reasons:
1. Each time we use this script moving forward, we will never need to edit the `RIGHT_READS`, `BAM_FILE` and `REMOVED_DUPLICATES_BAM_FILE` variables
2. Reduces the chance for typos


#### Brieft `bash` Text Manipulation 

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

**ID**: This is the identification for a given batch of reads. This ***MUST*** be unique to your experiment. 

**PL**: This is the platform that the sequencing was run on. For aligning Illumina reads, you should use "illumina" here.

**PU**: This is the platform unit and it is ideally supposed to hold `<FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>`, where `<FLOWCELL_BARCODE>` is the barcode of the flowcell, `<LANE>` is the lane the data was run on and `<SAMPLE_BARCODE>` is supposed to be a library/sample specific identifer. In some software packages **PU** can take precedence over the **ID** field. If you don't happen to have the `<FLOWCELL_BARCODE>.<LANE>.<SAMPLE_BARCODE>`, just make this field something useful that will help identify the sample. In this case, we didn't have that information so we are re-using the **ID** field here. 

**SM**: This is to mark which *sample* your reads are coming from. Note, this does not need to be unique like the **ID** field, since you may have multiple Read Group **ID**s coming from a single sample. For example, you may be merging alignment (BAM/SAM) files that originated from the same individual, but we sequenced on different lanes or machines. In this case, those would have different Read Group **ID**s, but the same **SM** value.

More information about Read Groups and some fields we didn't discuss can be found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups).
  
- `$REFERENCE_SEQUENCE` This calls the variable that holds the path to our reference sequence (`/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa`).
  
- `$LEFT_READS` and `$RIGHT_READS` These are variables with paths to the FASTQ sequence read files that we wish to align (`~/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_1.fq.gz` and `~/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_2.fq.gz`). In this case you can see that the reads are paired-end reads because 1) we have provided two sets of reads and 2) they have `_1`/`_2` right before the `.fq`/`.fastq` extension. This is a very common annotation for pair-end reads. Alternatively, sometimes they might also be annotated as `_R1`/`_R2`, but once again, they will almost always be placed right before the `.fq`/`.fastq` extension. We can also, note that these reads are currently compressed with gzip as annotated by the `.gz` extension. This is not a problem as `bwa` can read both gzip compressed and uncompressed FASTQ files.

- `-o $SAM_FILE` This calls a variable that holds the path to your output file (`~/variant_calling/alignments/normal_GRCh38.p7.sam`) where your alignments written. Notice this is a SAM file (an uncompressed alignment file).

Now you have written your command to run `bwa` you are ready to run alignment. However, there are a few steps that we are going to add to the script so that they run immediately after the alignemnt finishes.

## Sorting Reads and SAM to BAM Conversion

Now that we have aligned our reads to the reference, we would see if we opened up the SAM file that the alignments are in the order they were processed by `bwa` and not in any particular order that would be useful for downstream analyses. So, we are going to sort them into a coordinate ordered format that downstream tools can use. 

Additionally, while SAM files are nice due to their human readability, they are typically quite large files and it is not an efficient use of space on the cluster. Fortunately, there is a binary compression version of SAM called BAM. While we sort the reads, we are going to use a very common tool, `samtools`, to convert our SAM files to BAM files. 

```
# Sort SAM file and convert it to a BAM file
samtools sort \
-@ 8 \
-O bam \
-o $BAM_FILE \
$SAM_FILE
```

Let's go ahead and break down this command:

`samtools sort` This calls the sort function within `samtools`.

`-@ 8` This tells `samtools` to use 8 threads when it multithreads this task. Since we requested 8 cores for this `sbatch` submission, let's go ahead and use them all.

`-O bam` This is declaring the output format of `bam`.

`-o $BAM_FILE` This is a `bash` variable that holds the path to the output file of the `samtools sort` command.

`$SAM_FILE` This is a `bash` variable holding the path to the input SAM file.

## Indexing the `BAM` File

Similar to the index of a book. Many software packages want an index of your BAM file in order to facilitate fast look-ups of a BAM file. While not all software packages that use a BAM file will require this, many will and thus it is a good practice to just go ahead and index your BAM file at this point. The command to index our BAM files with `samtools` is:

```
# Index the BAM file
samtools index \
-@ 8 \
$BAM_FILE
```

`samtools index` Calls the `index` function within `samtools`

`-@ 8` This tells `samtools` to use 8 threads when it multithreads this task. This process is usually quite fast, so it is a bit overkill to multithread it, but we had requested the cores for the rest of the job, so we might as well use them here as well.

`$BAM_FILE` This is a `bash` variable that holds the path to the BAM file that we want to index.

We don't need to provide an output file for `samtools index`, by default it will generate a new file using the same path and filename as the BAM file, but add `bai` as the extension to denote that it is a BAM-index file.


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
LEFT_READS=/n/data1/cores/bcbio/gammerdinger/variant_calling/fastq_files/synthetic_challenge_set3_normal_NGv3_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/data1/cores/bcbio/gammerdinger/variant_calling/alignments/normal_GRCh38.p7.sam
BAM_FILE=`echo ${SAM_FILE%sam}bam`
REMOVED_DUPLICATES_BAM_FILE=`echo ${BAM_FILE%bam}removed_duplicates.bam`

# Align reads with bwa
bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-normal\tPL:illumina\tPU:syn3-normal\tSM:syn3-normal' \
$REFERENCE_SEQUENCE \
$LEFT_READS \
$RIGHT_READS \
-o $SAM_FILE

# Sort SAM file and convert it to a BAM file
samtools sort \
-@ 8 \
-O bam \
-o $BAM_FILE \
$SAM_FILE

# Index the BAM file
samtools index \
-@ 8 \
$BAM_FILE
```

## Creating Tumor `sbatch` script

Now that we have created the `sbatch` script for our normal samples, we need to repeat the process for our tumor samples. All of the parameters will stay the same, we just need to edit the SBATCH error file, SBATCH output file, LEFT_READS variable, RIGHT_READS variable, SAM_FILE variable and the read group information within the `bwa` command. You could very well do this by hand and it would be just fine. However, to cut down on typos we are going to use `sed`. `sed` is a powerful tool within `bash` and [has a wide variety of applications](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/sed.html). However, one of the most common uses for `sed` is as a "find-and-replace" tool. The syntax for this type of task is:

```
sed 's/pattern/replacement/g' file.txt 
```

The `s` before `/pattern/replacement/` is telling `sed` that we are going to use its **substittion** function and the `g` after `/pattern/replacement/` is telling `sed` that we want to apply that change **globally**, or every instance in the file. `pattern` represents the pattern that we are looking for and `replacement` is what we wish to replace the `pattern` with. Lastly, we need to provide `sed` some text source to apply this "find-and-replace" function, so we have provided it with `file.txt`, but you can also pipe in a string or file to apply this function to. 

In our case, we are hoping to replace each instance of "normal" with "tumor". Therefore, we could call `sed` to do this using:

```
sed 's/normal/tumor/g' bwa_alignment_samtools_sorting_index_normal.sbatch 
```

We can see that all instances of "normal" have been replaced with "tumor". Now we would like to redirect this output to a file called `bwa_alignment_samtools_sorting_index_tumor.sbatch` rather than standard output, se we need to add redirection to the end of out command:

```
sed 's/normal/tumor/g' bwa_alignment_samtools_sorting_index_normal.sbatch > bwa_alignment_samtools_sorting_index_tumor.sbatch
```

If we look at the out

```
module load gcc/6.2.0 bwa/0.7.17 samtools/1.15.1 

bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-tumor\tPL:illumina\tPU:syn3-tumor\tSM:syn3-tumor' \
/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa \
fastq_files/synthetic_challenge_set3_tumor_NGv3_1.fq.gz \
fastq_files/synthetic_challenge_set3_tumor_NGv3_2.fq.gz \
-o alignments/tumor_GRCh38.p7.sam

samtools sort \
-@ 8 \
-O bam \
-o alignments/tumor_sorted_GRCh38.p7.bam \
alignments/tumor_GRCh38.p7.sam

samtools index \
alignments/tumor_sorted_GRCh38.p7.bam
```

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
