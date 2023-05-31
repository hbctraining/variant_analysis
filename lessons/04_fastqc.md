# Evaluating Read Qualities with FastQC

## Learning objectives
- Implement `FastQC` to evaluate read qualities
- Manipulate strings of `bash` variable
- Utilize `sed` to find-and-replace text

## Importance of Evaluating Read Qualities

Before engaging in any next-generation sequencing project is it best practice to inspect your sequence reads to ensure that they are of high-quality. Sources of error are be numerous and include poor library construction and a malfunctioning sequencer. Therefore, it is critically important that you analyze your sequenced reads to ensure that they are high-quality before you devote time and resources to downstream analyses.

<p align="center">
<img src="../img/Read_QC_Pipeline.png" width="800">
</p>

## Unmapped read data (FASTQ)

The [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format is the *de facto* file format for sequence reads generated from next-generation sequencing technologies. This file format evolved from FASTA in that it contains sequence data, but also contains quality information. Similar to FASTA, the FASTQ file begins with a header line. The difference is that the FASTQ header is denoted by a `@` character. For a single record (sequence read), there are four lines, each of which are described below:

|Line|Description|
|----|-----------|
|1|Always begins with '@', followed by information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+', and sometimes the same info as in line 1|
|4|Has a string of characters representing the quality scores; must have same number of characters as line 2|

Let's use the following read as an example:

```
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
```

The line 4 has characters encoding the quality of each nucleotide in the read. The legend below provides the mapping of quality scores (Phred-33) to the quality encoding characters. *Different quality encoding scales exist (differing by offset in the ASCII table), but note the most commonly used one is fastqsanger, which is the scale output by Illumina since mid-2011.* 
 ```
 Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                   |         |         |         |         |
    Quality score: 0........10........20........30........40                                
```
 
Using the above quality encoding character legend, the first nucelotide in the read (C) is called with a quality score of 31 (corresponding to encoding character `@`), and our Ns are called with a score of 2 (corresponding to encoding character `#`). **As you can tell by now, this is a bad read.** 

Each PHRED quality score represents the probability that the corresponding nucleotide call is incorrect, with higher PHRED scores representing lower probabilities of incorrect base calls. This quality score is logarithmically based and is calculated as:

	Q = -10 x log10(P), where P is the probability that a base call is erroneous

These probabaility values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------:|:---------------------------------:|:-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|

Therefore, for the first nucleotide in the read (C), there is less than a 1 in 1000 chance that the base was called incorrectly. Whereas, for the the end of the read there is greater than 50% probabaility that the base is called incorrectly.

Before we discuss implementing `FastQC` we are going to introduce string manipulation in `bash` which we will be using throughout our work.

***

**Exercise**

**1.** If the probability of a incorrect base call is 1 in 3,981, what is the associated PHRED score?

***

## FastQC

Now we understand what information is stored in a FASTQ file, the next step is to generate quality metrics for our sequence data.

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular **tool for analyzing read quality for NGS data**. It can evaluate many aspects of your NGS data including:
- Read quality by position
- GC distribution
- Overrepresented sequences
- More

When working in a cluster environment, you will find that generally many tools and software are pre-installed for your use. On the O2 cluster, these tools are available through the LMOD system. Let's first check to see if the tool FASTQC exists as a module:

```
$ module avail fastqc
```

We can decide on the version we would like to use and go ahead and load the FastQC module to use:

```
$ module load fastqc/0.11.9
```

You should now see that the module is loaded:

Now that we have loaded the module, FASTQC is directly available to you like any other basic UNIX command. That is, at the command prompt we just need to provide the name of the tool to use it. _This is because the path to the executable file for FastQC has now been added to our $PATH variable._ 

> Check your $PATH variable to see whether or not you see a relevant path. Is it appended to the beginning or end? Do you see any additional paths added?

To run FastQC we need to specify two arguments: 
1. the file name(s) of our FASTQ input (one file for single-end; two files for paired-end)
2. the directory where the results (ouput) will be stored, which is indicated after the -o flag

**Example code is provided below. DO NOT RUN!**

```
## DO NOT RUN!

$ fastqc -o ~/variant_calling/results/fastqc/ \
      --threads 4 \
     ~/variant_calling/raw_data/syn3_normal_1.fq.gz variant_calling/raw_data/syn3_normal_2.fq.gz 
```

This command is pretty strightforward, but we will explain each part:

- `fastqc` This calls the `FastQC` software package
- `--outdir` or `-o`: This is the directory for the output files to be written to
- `--threads` This specifies the number of threads that `FastQC` can use to speed up the processing
-  This is the left and right read (R1 and R2) from the input FASTQ file (paired-end data)


This would be fine if we just had a one or two samples to run. But a typical dataset can contain 10s to 100s of samples, and you don't want to be sitting around running each sample interactively. **It would be much more efficient if we submitted this as a job.**

### Running FastQC for normal samples

In order to run `FastQC` we are going to develop an `sbatch` script to run on the cluster. First, we will need to change directories to where our scripts are stored and start a new submission script using `vim`:

```
$ cd ~/variant_calling/scripts/
$ vim fastqc_normal.sbatch
```

Now that we have opened up `vim`, we need to enter `insert-mode` by pressing `i`. Once in `insert-mode`, we will add our shebang line, description and `sbatch` directives. These directives indicate to SLURM what resources we need in order to run the job.

```
#!/bin/bash
# This sbatch script is for running FastQC to evaluate read qualities

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -o fastqc_normal_%j.out
#SBATCH -e fastqc_normal_%j.err
```

Next, we will add loading the `FastQC` module to our script:

```
# Load module
module load fastqc/0.11.9
```

Next, we will **define the variables that we will be using**. This step isn't necessary and you can alternatively just manually input full paths to files and parameters within the actual command. However, the use of variables instead has some distinct advantages:

1. If it is a variable, you only need to type it once to define the value instead of re-writing each occurance of the value.
2. The use of variable allows us to use `bash` string manipulation to assign variables rather than typing each out
3. Variables help keep your commands cleaner and more organized
4. Variables help make repurposing your script for a different project easier, since it might just be the variables changing and you don't need to even touch the actual command

Hopefully, the case made for assigning varibles outside of your command has been successful. The **variables we will be adding to this script** are:

```
# Assign variables
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/reports/fastqc/syn3_normal/
THREADS=4
```

Notice here how **we are using string manipulation in bash to assign the `$RIGHT_READS` variable.** The left and right reads from paired-end sequencing will oftentimes be in the same directory. So this is a case where string manipulation in `bash` will be really helpful in saving time and reducing typos. 

> #### String manipulation in bash using `%` 
> Here, we introduce the % tool for text manipulation. % is placed within the {} of a variable and tells bash to remove the shortest match from the end that contains the text that follows the %. So in our case, the substring "1.fq.gz" will get removed. Then outside of the `{}` we have added the extension for the second read file. For more practice with string manipulation in bash [check out these materials]().

Before we run `FastQC` we need to make sure that the creation of the output directory exists in our script to accept the output:

```
# Create directory to hold output
mkdir -p $OUTPUT_DIRECTORY
```

The `-p` option for `mkdir` does two things:

1. It will make any parent directories necessary
2. If the directory already exists it will not return an warning message.

Now we can add the code to run `FastQC`, an use the variables we had created:

```
# Run FastQC
fastqc \
$LEFT_READS \
$RIGHT_READS \
--outdir $OUTPUT_DIRECTORY \
--threads $THREADS
```

All together, the final `sbatch` script should look like:

```
#!/bin/bash
# This sbatch script is for running FastQC to evaluate read qualities

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -o fastqc_normal_%j.out
#SBATCH -e fastqc_normal_%j.err

# Load module
module load fastqc/0.11.9

# Assign variables
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/reports/fastqc/syn3_normal/
THREADS=4

# Create directory to hold output
mkdir -p $OUTPUT_DIRECTORY

# Run FastQC
fastqc \
$LEFT_READS \
$RIGHT_READS \
--outdir $OUTPUT_DIRECTORY \
--threads $THREADS
```

Now we can submit this script to the cluster:

```
sbatch fastqc_normal.sbatch
```

### Running FastQC for the tumor sample

Now that we have created the `sbatch` script for our normal samples, we need to **repeat the process for our tumor samples**. All of the parameters will stay the same, we just need to edit the SBATCH error file, SBATCH output file and LEFT_READS variable. You could very well do this by hand and it would be just fine. However, to cut down on typos we are going to use `sed`. `sed` is a powerful tool within `bash` and [has a wide variety of applications](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/sed.html). 

We will create the new script by **using `sed` as a "find-and-replace" tool**. The syntax for this type of task is:

```
sed 's/pattern/replacement/g' file.txt 
```

* The `s` before `/pattern/replacement/` is telling `sed` that we are going to use its **substitution** function 
* The `g` after `/pattern/replacement/` is telling `sed` that we want to apply that change **globally**, or every instance in the file. 
* `pattern` represents the pattern that we are looking for 
* `replacement` is what we wish to replace the `pattern` with. 

In our case, **we are replacing each instance of "normal" with "tumor"**. Therefore, we could call `sed` to do this using:

```
$ sed 's/normal/tumor/g' fastqc_normal.sbatch
```

We can see that all instances of "normal" have been replaced with "tumor". Now we would like to **redirect this output to a file called `fastqc_tumor.sbatch`** rather than standard output, so we need to add redirection to the end of out command:

```
$ sed 's/normal/tumor/g' fastqc_normal.sbatch >  fastqc_tumor.sbatch
```

Now, we have the tumor sample to analyze with `FastQC`, so let's inspect the file to make sure it is right:

```
$ cat fastqc_tumor.sbatch
```

It should looke like:

```
#!/bin/bash
# This sbatch script is for running FastQC to evaluate read qualities

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 4
#SBATCH --mem 8G
#SBATCH -o fastqc_tumor_%j.out
#SBATCH -e fastqc_tumor_%j.err

# Load module
module load fastqc/0.11.9

# Assign variables
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_tumor_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/reports/fastqc/syn3_tumor/
THREADS=4

# Create directory to hold output
mkdir -p $OUTPUT_DIRECTORY

# Run FastQC
fastqc \
$LEFT_READS \
$RIGHT_READS \
--outdir $OUTPUT_DIRECTORY \
--threads $THREADS
```

Once we have looked it over and it matches what we have, we can go ahead and **submit this `SBATCH` submission script to the cluster**:

```
sbatch fastqc_tumor.sbatch
```

Traditionally, most people inspect their `FastQC` reports before continuing on with their analysis. However, we are going to merge all of our QC together using `MultiQC` and analyze it all after alignment.


***

[Next Lesson >>](05_sequence_alignment_theory.md)

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
