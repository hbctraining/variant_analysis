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

Next, we will define our variables:

```
LEFT_READS=/home/$USER/variant_calling/raw_data/syn3_normal_1.fq.gz
RIGHT_READS=`echo ${LEFT_READS%1.fq.gz}2.fq.gz`
OUTPUT_DIRECTORY=~/variant_calling/fastqc/normal/
```

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
-t 4
```

This command is pretty strightforward, but we will explain each part:

- `fastqc` This calls the `FastQC` software package
- `$LEFT_READS` This is the left read (R1 or Read 1) input FASTQ file
- `$Right_READS` This is the right read (R2 or Read 2) input FASTQ file
- `-o $OUTPUT_DIRECTORY` This is the directory for the output files to be written to
- `-t 4` This specifies the number of threads that `FastQC` can use to speed up the processing

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
OUTPUT_DIRECTORY=~/variant_calling/fastqc/normal/

mkdir -p $OUTPUT_DIRECTORY

fastqc \
$LEFT_READS \
$RIGHT_READS \
-o $OUTPUT_DIRECTORY \
-t 4
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

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
```
