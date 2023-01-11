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

<details> 
<summary>Click here for details on creating a <code>bwa</code> index</summary>

While we will be using an index that has already been made for us, if you need to create and index for a reference sequence using <code>bwa</code>, the steps for this are laid out below.

<ol><b><li>Navigate to your reference sequence directory</b>

<pre>
cd ~/path/to/reference/sequence/directory/
</pre></li>

<b><li>Create a `bwa` index</b>

<pre>
bwa index reference_sequence.fasta
</pre></li>

This process may take up to 30+ minutes to run depending on the reference sequence size. The output of this command will produce five files:

<ul><li><code>reference_sequence.fasta.sa</code> This is the binary version of the Suffix Array Index</li>
<li><code>reference_sequence.fasta.bwt</code> This is the binary version of the Burrows-Wheeler Transform of the reference sequence</li>
<li><code>reference_sequence.fasta.pac</code> A special binary compression of the reference sequence</li>
<li><code>reference_sequence.fasta.ann</code> Notations regarding the reference sequence</li>
<li><code>reference_sequence.fasta.amb</code> Notations regarding base ambiguities (mostly Ns, but also other base ambiguities) in the reference sequence</li></ul>
</details>

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

## Creating Tumor `sbatch` script

Now that we have created the `sbatch` script for our normal samples, we need to repeat the process for our tumor samples. All of the parameters will stay the same, we just need to edit the SBATCH error file, SBATCH output file, LEFT_READS variable, RIGHT_READS variable, SAM_FILE variable and the read group information within the `bwa` command. You could very well do this by hand and it would be just fine. However, to cut down on typos we are going to use `sed`. `sed` is a powerful tool within `bash` and [has a wide variety of applications](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/sed.html). However, one of the most common uses for `sed` is as a "find-and-replace" tool. The syntax for this type of task is:

```
sed 's/pattern/replacement/g' file.txt 
```

The `s` before `/pattern/replacement/` is telling `sed` that we are going to use its **substittion** function and the `g` after `/pattern/replacement/` is telling `sed` that we want to apply that change **globally**, or every instance in the file. `pattern` represents the pattern that we are looking for and `replacement` is what we wish to replace the `pattern` with. Lastly, we need to provide `sed` some text source to apply this "find-and-replace" function, so we have provided it with `file.txt`, but you can also pipe in a string or file to apply this function to. 

In our case, we are hoping to replace each instance of "normal" with "tumor". Therefore, we could call `sed` to do this using:

```
sed 's/normal/tumor/g' bwa_alignment_tumor.sbatch
```

We can see that all instances of "normal" have been replaced with "tumor". Now we would like to redirect this output to a file called `bwa_alignment_tumor.sbatch` rather than standard output, se we need to add redirection to the end of out command:

```
sed 's/normal/tumor/g' bwa_alignment_normal.sbatch >  bwa_alignment_tumor.sbatch
```

If we look at the output with:

```
cat bwa_alignment_tumor.sbatch 
```

Then it should look like this:

```
#!/bin/bash
# This script is for aligning sequencing reads against a reference genome using bwa

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o bwa_alignment_tumor_%j.out
#SBATCH -e bwa_alignment_tumor_%j.err

# Load modules
module load gcc/6.2.0
module load bwa/0.7.17

# Assign files to bash variables
REFERENCE_SEQUENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_tumor_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
SAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/tumor_GRCh38.p7.sam

# Align reads with bwa
bwa mem \
-M \
-t 8 \
-R '@RG\tID:syn3-tumor\tPL:illumina\tPU:syn3-tumor\tSM:syn3-tumor' \
$REFERENCE_SEQUENCE \
$LEFT_READS \
$RIGHT_READS \
-o $SAM_FILE
```

Once we have created this script we can go ahead and submit it for processing:

```
sbatch bwa_alignment_tumor.sbatch 
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

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
