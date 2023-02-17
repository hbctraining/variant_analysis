# Workshop Schedule

## Pre-reading

**Before the workshop:**

Please read the following page to learn about the dataset we will be using:

[ICGC-TCGA DREAM Mutation Calling Challenge Synthetic Dataset](../lessons/01_syn3_dataset.md)

## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:10 | [Workshop Introduction](../lectures/Variant_calling_Intro_to_workshop_all.pdf) | Will |
| 10:00 - 11:30 | [Introduction to Variant Calling]() | Sergey |
| 11:30 - 11:50 | [Project Organization](../lessons/02_project_organization.md) | Will |
| 11:50 - 12:00 | Overview of self-learning materials and homework submission | Will |

### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:

  1. [File Formats](../lessons/03_file_formats.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Before we dive too deeply into calling variants, there are a few file formats that we will see during our analysis. Understanding how these files are formatted will allow you to inspect them to ensure that the software programs that we are employing are working correctly.
         <br><br>This lesson will cover:<br>
             <ul><li>Describe the difference between 0-based and 1-based indexing</li>
             <li>Decode a FLAG in a SAM file in order to reveal information about the nature of the read's alignment</li>
             <li>Create a CIGAR string for an alignment</li>
             <li>Parse out variant information from a VCF file</li>
             <li>Create a BED file</li></ul>
             <hr />
        </details>
        
  2. [Evaluating Read Quality with `FastQC`](../lessons/04_fastqc.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>The first step in many NGS studies is first to evaluate the read qualites that you received from the sequencing facility. A common tool used for handling this analysis is <code>FastQC</code>. 
         <br><br>This lesson will cover:<br>
          <ul><li>Implement FastQC to evaluate read qualities</li>
          <li>Manipulate strings of bash variable</li>
          <li>Evaluate FastQC output</li>
          <li>Utilize sed to find-and-replace text</li></ul>
          <hr />
        </details>

  3. [Sequence Read Alignment](../lessons/05_sequence_alignment_theory.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Once we have completed our QC on sequence reads we will be aligning the reads to a reference sequence. This alignment step places each read in genomic space and creates the bedrock for calling variants.
         <br><br>This lesson will cover:<br>
             <ul><li>Enumerate difficulties with alignment</li>
             <li>Create an <code>sbatch</code> script to align reads</li></ul>
             <hr />
        </details>

  4. [Alignment File Processing ](../lessons/06_alignment_file_processing.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Before we can call variants from our alignment files, we need to do some processing to clean them up. The two major concerns here are organizing (sorting) our alignment files for our analyses and removing duplicates.
         <br><br>This lesson will cover:<br>
             <ul><li>Differentiate between query-sorted and coordinate-sorted alignment files</li>
             <li>Describe and remove duplicate reads</li>
             <li>Process a raw SAM file for input into a BAM for GATK</li></ul>
             <hr />
        </details>

  5. [Alignment File Quality Control](../lessons/07_alignment_QC.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Once we have our alignment files processed, we want to evaluate them to ensure that the data is of high-quality before proceeding into variant calling. We also need to merge our read quality QC from <code>FastQC</code> into a report with these alignment QC metrics using <code>MultiQC</code>.
         <br><br>This lesson will cover:<br>
             <ul><li>Estimate alignment rates using <code>Picard</code></li>
             <li>Merge <code>Picard</code> QC metrics with <code>FastQC</code> metrics using <code>MultiQC</code></li></ul>
             <hr />
        </details>
        
  6. [Evaluating Quality Control Metrics](../lessons/08_evaluate_QC.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Many high-performance computing clusters are not designed to render the HTML reports produced by <code>MultiQC</code>. Thus, we will use <code>FileZilla</code> to dowload our <code>MultiQC</code> HTML report and interpret the results within it.
         <br><br>This lesson will cover:<br>
             <ul><li>Evaluating alignment rates</li>
             <li>Intepretting read QC metrics within <code>MultiQC</code> HTML report</li></ul>
             <hr />
        </details>
        


> **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it).
> 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (replace the "XX" in the username with the number you were [assigned in class](https://docs.google.com/spreadsheets/d/1kBlYowhjjHJC9ZovmbBULmbqozKpprM17vZ2wPlhNg0/edit?usp=sharing)). 
> 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:30 --mem 1G /bin/bash` to get on a compute node.
> 3. Proceed once your command prompt has the word `compute` in it.
> 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

II. **Complete the exercises**:
   * Each lesson above contains exercises; please go through each of them.
   * **Copy over** your solutions into the [Google Form]() the **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:15 | Self-learning lessons review | All |
| 10:15 - 11:00 | [Variant Calling](../lessons/09_variant_calling.md) | Will |
| 11:00 - 11:30 | [Variant Filtering](../lessons/10_variant_filtering.md) | Will |
| 11:30 - 12:00 | [Variant Annotation with SnpEff](../lessons/11_variant_annotation.md) | Will |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:

1. [Automation of Variant Calling Pipeline](../lessons/12_automation_of_variant_calling.md)

      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Now that we have executed much of the standard workflow for variant calling, we might want to automate our workflow to make future analyses more streamlined and reproducible with a consistent workflow. We will need to adapt our current scripts to allow for a more streamlined workflow and also discuss some intricacies <code>bash</code> and <code>SLURM</code> that will help us create this automated workflow.
         <br><br>This lesson will cover:<br>
          <ul><li>Construct a flexible pipeline for automating variant calling</li>
          <li>Integrate the <code>--dependency</code> option for <code>sbatch</code> into workflows</li></ul>
          <hr />
        </details>

> **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it). For login instructions, please see above.

II. **Complete the exercises**:
   * Each lesson above contains exercises; please go through each of them.
   * **Copy over** your solutions into the [Google Form]() the **day before the next class**.
   
### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:00 | Self-learning lessons review | All |
| 10:00 - 10:45 | [Variant Prioritization with SnpSift](../lessons/13_variant_prioritization.md) | Will |
| 10:45 - 11:25 | [Visualization in IGV](../lessons/14_IGV.md) | Will |
| 11:25 - 11:45 | Q & A | Will |
| 11:45 - 12:00 | [Wrap up](../lectures/Variant_calling_wrapup_all.pdf) | Will |

***

## Answer key
* [Exercises Answer Key](../answer_key/Answer_key.md)

***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
