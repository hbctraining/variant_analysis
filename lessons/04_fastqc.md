# Evaluating Read Qualities with `FastQC`

## Learning objectives
- Implement `FastQC` to evaluate read qualities
- Manipulate strings of `bash` variable
- Evaluate `FastQC` output
- Utilize `sed` to find-and-replace text

## Importance of Evaluating Read Qualities

Before engaging in any next-generation sequencing project is it best practice to inspect your sequence reads to ensure that they are of good quality. Sources of error are be numerous and include poor library construction and a malfunctioning sequencer. Therefore, it is critically important that you analyze your sequenced reads to ensure that they are high-quality before you devote time and resources to downstream analyses.

<p align="center">
<img src="../img/Read_QC_Pipeline.png" width="800">
</p>

Before we discuss implementing `FastQC` we are going to introduce string manipulation in `bash` which we will be using throughout our work.

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
 
Using the quality encoding character legend, the first nucelotide in the read (C) is called with a quality score of 31 (corresponding to encoding character `@`), and our Ns are called with a score of 2 (corresponding to encoding character `#`). **As you can tell by now, this is a bad read.** 

Each quality score represents the probability that the corresponding nucleotide call is incorrect. This quality score is logarithmically based and is calculated as:

	Q = -10 x log10(P), where P is the probability that a base call is erroneous

These probabaility values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------:|:---------------------------------:|:-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|

Therefore, for the first nucleotide in the read (C), there is less than a 1 in 1000 chance that the base was called incorrectly. Whereas, for the the end of the read there is greater than 50% probabaility that the base is called incorrectly.

## `bash` String Manipulation

String manipulation in `bash` can be a very helpful tool in minimizing typos whenever evaluating a script that uses `bash`. Before we discuss manipulating strings, we should first define a string. A string is a data type that is used to represent text rather than integers. Examples of strings include:

- "Hello World"
- "/path/to/file.txt"
- "TPS_Report_2019"

With this understanding of what strings are, we might start to see some value in them, particularly with respect to filepaths. Imagine you have two filepaths:

```
/This/is/my/extremely/super/duper/long/and/annoying/filepath/file_1.txt
/This/is/my/extremely/super/duper/long/and/annoying/filepath/file_2.txt
``` 

Alternatively, we could use `bash` variables and string manipulation:

```
FILE_1=/This/is/my/extremely/super/duper/long/and/annoying/filepath/file_1.txt
FILE_2=${FILE_1%1.txt}2.txt
```

We will explain the syntax in a bit, but we might chose to do it this way for three reasons:
1. As long as File 1 is always in the same directory as File 2 and they both end with `1.txt` and `2.txt`, respectively, then each time we use this script moving forward, we will never need to edit the `FILE_2` variable. 
2. Reduces the chance for typos
3. Can also help keep filename nonmenclature consistent across files

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
  <hr />
</details>

## `FastQC`

[`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular tool for analyzing read quality for NGS data. It can evaluate many aspects of your NGS data including:
- Read quality by position
- GC distribution
- Overrepresented sequences
- More

### Running `FastQC` for normal samples

In order to run `FastQC` we are going to develop an `sbatch` script to run on the cluster. First, we will need to change directories to where our scripts are stored and start a new submission script using `vim`:

```
cd ~/variant_calling/scripts/
vim fastqc_normal.sbatch
```

First, we will add our shebang line, description and `sbatch` directives.

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

Next, we will load the `FastQC` module:

```
module load fastqc/0.11.9
```

Next, we will define the variables that we will be using. This step isn't necessary and you can alternatively just manually input full paths to files and parameters within the actual command. However, the use of variables instead has some distinct advantages:

1) Anything used more than once could have more easily been a variable and less prone to typos. Perhaps this is a common filepath or sample name that you need to write many times. If it is a variable, you only need to type it once instead of each occurance.
2) The use of variable allows us to use `bash` string manipulation to assign variables rather than typing each out
3) Can help keep your commands cleaner and more organized
4) Can help make repurposing your script for a different project easier, since it might just be the variables changing and you don't need to even touch the actual command

Hopefully, the case made for assigning varibles outside of your command has been successful. The variables we will be setting for this script are:

```
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/reports/fastqc/syn3_normal/
THREADS=4
```

Notice here how we are using string manipulation in `bash` to assign `$RIGHT_READS`. The left and right reads from paired-end sequencing will oftentimes be in the same directory. So this is a case where string manipulation in `bash` will be really helpful in saving time and reducing typos. 

Now that we have assigned parameters to variables, we are going to get ready to run `FastQC` and before we run `FastQC` we need to make sure that the output directory exists to accept the output:

```
mkdir -p $OUTPUT_DIRECTORY
```

The `-p` option for `mkdir` does two things:

1) It will make any parent directories necessary
2) If the directory already exits it will not return an warning message.

Now we can run `FastQC`:

```
fastqc \
$LEFT_READS \
$RIGHT_READS \
-o $OUTPUT_DIRECTORY \
-t $THREADS
```

This command is pretty strightforward, but we will explain each part:

- `fastqc` This calls the `FastQC` software package
- `$LEFT_READS` This is the left read (R1 or Read 1) input FASTQ file
- `$Right_READS` This is the right read (R2 or Read 2) input FASTQ file
- `-o $OUTPUT_DIRECTORY` This is the directory for the output files to be written to
- `-t $THREADS` This specifies the number of threads that `FastQC` can use to speed up the processing

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

module load fastqc/0.11.9

LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/reports/fastqc/syn3_normal/
THREADS=4

mkdir -p $OUTPUT_DIRECTORY

fastqc \
$LEFT_READS \
$RIGHT_READS \
-o $OUTPUT_DIRECTORY \
-t $THREADS
```

Now we can submit this script to the cluster:

```
sbatch fastqc_normal.sbatch
```

### Running `FastQC` for the tumor sample

Now that we have created the `sbatch` script for our normal samples, we need to repeat the process for our tumor samples. All of the parameters will stay the same, we just need to edit the SBATCH error file, SBATCH output file and LEFT_READS variable. You could very well do this by hand and it would be just fine. However, to cut down on typos we are going to use `sed`. `sed` is a powerful tool within `bash` and [has a wide variety of applications](https://hbctraining.github.io/Training-modules/Intermediate_shell/lessons/sed.html). However, one of the most common uses for `sed` is as a "find-and-replace" tool. The syntax for this type of task is:

```
sed 's/pattern/replacement/g' file.txt 
```

The `s` before `/pattern/replacement/` is telling `sed` that we are going to use its **substittion** function and the `g` after `/pattern/replacement/` is telling `sed` that we want to apply that change **globally**, or every instance in the file. `pattern` represents the pattern that we are looking for and `replacement` is what we wish to replace the `pattern` with. Lastly, we need to provide `sed` some text source to apply this "find-and-replace" function, so we have provided it with `file.txt`, but you can also pipe in a string or file to apply this function to. 

In our case, we are hoping to replace each instance of "normal" with "tumor". Therefore, we could call `sed` to do this using:

```
sed 's/normal/tumor/g' fastqc_normal.sbatch
```

We can see that all instances of "normal" have been replaced with "tumor". Now we would like to redirect this output to a file called `fastqc_tumor.sbatch` rather than standard output, se we need to add redirection to the end of out command:

```
sed 's/normal/tumor/g' fastqc_normal.sbatch >  fastqc_tumor.sbatch
```

Now, we have the tumor sample to analyze with FastQC:

```
sbatch fastqc_tumor.sbatch
```

Traditionally, most people inspect their `FastQC` reports before continuing on with their analysis. However, we are going to merge all of our QC together using `MultiQC` and analyze it all after alignment.

[Next Lesson >>>](05_sequence_alignment_theory.md)

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
