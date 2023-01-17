# Alignment Quality Control

## Learning Objectives

- Verify alginment rates using `Picard`

## Factors Impacting Alignment

One of the most important metrics for your alignment file is the alignment rate. Alignment rates can vary based upon many factors, including:

- **Quality of reference assembly** - A high-quality assembly like GRCh38.p7 will provide a more than adequate reference geneome for alignment. However, if you were studying a organism with a poorly assembled genome, parts of the reference genome could be missing from the assembly. Therefore, high-quality reads might not align because they there is missing reference sequence to align to that corresponds to their sequence.
- **Quality of libraries** - If the library generation was poor and there wasn't enough input DNA, then your sequencing could be filled with low-quality reads
- **Quality of the reads** - If the reads are poor quality, then it can make alignment more uncertain. If your `FASTQC` report shows any anomalous signs, contact your sequencing center for support.
- **Contamination** - If your samples are contaminated, then it can also skew your alignment. For example, if your samples were heavily contaminated with some bacteria, then much of what you will sequence will be bacteria DNA and not your sample DNA. As a result, most of the sequence reads will not align to your target sequence. If you suspect contamination might be the source of a poor alignment, you could consider running [Kraken](https://ccb.jhu.edu/software/kraken/) to evaluate the levels of contamination in your samples.
- **Evolutionary distance between the sampled organism and reference genome** If a reference genome doesn't exist for your species of interest, you are able to align reads to a closely-related organism. However, it does come at the cost of lowering the alignment rate. 
- **Aligner and alignment parameters** Different aligners work differently and are specialized for different types of data. Additionally, many aligners have a variety of parameters that are able to be adjusted. As a result, different aligners or different parameters for the same aligner will give different alignment rates, but they usually should be within the same approximate alignment rate. Generally speaking, the default parameters for most alignment tools are usually fine and they shouldn't need to be manually adjusted/optimized unless there is a specific reason to do so.

When aligning high-quality reads to a high quality reference genome, one should expect to see alignment rates at 90% or better. If alignment rates dipped below 80-85%, then there could be reason for furtheer inspection. 

## Collecting Alignment Statistics

We are going to use `Picard` once again in order to collect our alignment statistics. `Picard` has many packages for collecting different types of data, but the one we will be using is `CollectAlignmentSummaryMetrics`. Let's start creating an `sbatch` script that can utilize this package:

```
cd ~/variant_calling/scripts/
vim picard_CollectAlignmentMetrics_normal.sbatch
```

First, we need to add our shebang line, description and `sbatch` directives to the script:

```
#!/bin/bash
# This sbatch script is for collecting alignment metrics using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 1
#SBATCH --mem 4G
#SBATCH -o picard_CollectAlignmentMetrics_normal_%j.out
#SBATCH -e picard_CollectAlignmentMetrics_normal_%j.err
```

Next we need to load `Picard`:

```
# Load picard
module load picard/2.8.0
```

Next, let's assign our files to variables:

```
# Assign variables
COORDINATE_SORTED_BAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
METRICS_OUTPUT_FILE=/home/${USER}/variant_calling/reports/normal__GRCh38.p7.CollectAlignmentSummaryMetrics.txt
```

Next, we can add the `Picard` command to gather the alignment metrics:

```
# Run Picard CollectAlignmentSummaryMetrics
picard CollectAlignmentSummaryMetrics \
INPUT=$COORDINATE_SORTED_BAM_FILE
OUTPUT=$METRICS_OUTPUT_FILE
REFERENCE_SEQUENCE=$REFERENCE
```

We can breakdown this command into each of it's components:

- `picard CollectAlignmentSummaryMetrics` Calls the `CollectAlignmentSummaryMetrics` from within `Picard`

- `INPUT=$COORDINATE_SORTED_BAM_FILE` This is the output from our previous `Picard` alignment processing steps.

- `OUTPUT=$METRICS_OUTPUT_FILE` This is the file to write the output metrics to.

- `REFERENCE_SEQUENCE=$REFERENCE` This isn't a required parameter, but `picard` can do a subset of mismatch-related metrics if this is provided.

The `sbatch` submission script for collecting the alignment metrics should look like:

```
#!/bin/bash
# This sbatch script is for collecting alignment metrics using Picard 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:30:00
#SBATCH -c 1
#SBATCH --mem 4G
#SBATCH -o picard_CollectAlignmentMetrics_normal_%j.out
#SBATCH -e picard_CollectAlignmentMetrics_normal_%j.err

# Load picard
module load picard/2.8.0

# Assign variables
COORDINATE_SORTED_BAM_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/alignments/normal_GRCh38.p7.coordinate_sorted.bam
REFERENCE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7_genomic.fa
METRICS_OUTPUT_FILE=/home/${USER}/variant_calling/reports/normal_GRCh38.p7.CollectAlignmentSummaryMetrics.txt

# Run Picard CollectAlignmentSummaryMetrics
picard CollectAlignmentSummaryMetrics \
INPUT=$COORDINATE_SORTED_BAM_FILE
OUTPUT=$METRICS_OUTPUT_FILE
REFERENCE_SEQUENCE=$REFERENCE
```

Create the tumor version of this submission script using `sed`:

```
sed 's/normal/tumor/g' picard_CollectAlignmentMetrics_normal.sbatch > picard_CollectAlignmentMetrics_tumor.sbatch
```

Now that we have created these files to submit, let's check the status of our previous `Picard` alignment processing steps:

```
squeue -u $USER
```

**If your `Picard` alignment processing steps are not completed yet**, wait until they have finished before submitting these jobs to collect alignment metrics.

**If your `Picard` alignment processing steps are completed**, then submit these jobs to collect alignment metrics:

```
sbatch picard_CollectAlignmentMetrics_normal.sbatch
sbatch picard_CollectAlignmentMetrics_tumor.sbatch
```

## Options for Inspecting `Picard` Alignment Metrics

Once the job has finished we would inspect the output files. This could be done in one of a few ways:

1) View each metrics file in a `less` buffer
  **Pros:**
    - Simple to do
  **Cons:**
    - Hard to compare across samples
    - Tedious to parse the columns
    
2) Download each metrics from the O2 cluster and import them into Excel/Excel-like program that puts tab-delmited files into a grid
  **Pros:**
    - Easier to interpret the data than the `less` buffer approach
  **Cons:**
    - Hard to compare across samples
    - Have to download the metrics files from the O2 cluster
    
3) Collate metrics files using `MultiQC` and download the `MultiQC` HTML report from the O2 cluster
  **Pros:**
    - Alignment metrics are easy to compare across samples
    - Easy to interpret results
    - We can combine the `picard CollectAlignmentSummaryMetrics` output with our `FASTQC` reports
  **Cons:**
    - Have to download a file from the O2 cluster
    - Have to run samples through an extra `MultiQC` step

None of the above methods are wrong, but some are more elegant than others. One might use **Method 1)** if they only had a handful of samples (~<5) to analyze and only wanted a single statistic, like alignment rate, from each. Then, it might be fastest just to open them up in a `less` buffer. However, if one has lots of samples (>5) then the advantages of `MultiQC` collating the results starts to become really helpful and one might choose **Method 3)**. Unfortunately, `MultiQC` doesn't display *ALL* of the data contained in the metrics file, so one may be inclined to do **Method 2** and downloard the directory full of metrics files in order to view the metrics not included in the `MultiQC` report in a program like Microsoft Excel.

However, for this workshop, we are going to collate our results in `MultiQC` and download the HTML report to our local computers.

## Inspecting `Picard` Alignment Metrics

One nice feature of `MultiQC` is that it accepts many different file formats. It figures out which format was submitted and tailors the report to that type of analysis. Collating our `MultiQC` results would be relatively quick to just run from the command-line, but it's best practice to write our steps to scripts so that we always have a record of what we did and how we created our reports. We will start by writing a `sbatch` script in `vim` for submission:

```
vim multiqc_alignment_metrics.sbatch
```

First, we will add our sheband line, description and `sbatch` directives

```
#!/bin/bash
# This sbatch script is for collating alignment metrics from Picard using MultiQC 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:10:00
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -o multiqc_alignment_metrics_%j.out
#SBATCH -e multiqc_alignment_metrics_%j.err
```

Next, we will load out modules:

```
# Load modules
module load gcc/9.2.0
module load multiqc/1.12
```
> NOTE: `MultiQC` version 1.12 requires `gcc/9.2.0` on the O2 cluster.

Next, we will assign our variables:

```
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/
```

Then, we will add the command to run `MulktiQC`:

```
mulitqc $REPORTS_DIRECTORY
```

So our final `sbatch` script should look like:

```
#!/bin/bash
# This sbatch script is for collating alignment metrics from Picard using MultiQC 

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-00:10:00
#SBATCH -c 1
#SBATCH --mem 1G
#SBATCH -o multiqc_alignment_metrics_%j.out
#SBATCH -e multiqc_alignment_metrics_%j.err

# Load modules
module load gcc/9.2.0
module load multiqc/1.12

REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/

mulitqc $REPORTS_DIRECTORY
```

Like the previous step, we will need to check to ensure that the previous `Picard` step for collecting metrics for each sample is down before we can submit this script. To do this, we will check out `squeue`:

```
squeue -u $USER
```

**If your `Picard` collect alignment metric steps are not completed yet**, wait until they have finished before submitting these jobs to `MultiQC`.

**If your `Picard` collect alignment metric steps are completed**, then submit this `MultiQC` job to collate the alignment metrics:

```
sbatch multiqc_alignment_metrics.sbatch
```

This job should finish fairly quickly and then we can proceed to downloading it with `FileZilla`.

## Downloading `MultiQC` HTML Report with `FileZilla`

While the O2 cluster cluster is fantastic at many things, it is not designed to render HTML files. For that we will need a browser, such as Safari, Chrome, Firefox, etc., on our local computer. Thus, we will need to download the HTML report from the cluster to our local computers. There are ways to do this from the command line using tools like `scp` and `rsync`, however, we are going to use `FileZilla` which has an easy-to-use GUI to help us.

### Filezilla - Step 1

Open up *FileZilla*, and click on the File tab. Choose 'Site Manager'.

<p align="center">
<img src="../img/filezilla_setup.png" width="500">
</p>

### Filezilla - Step 2

Within the 'Site Manager' window, do the following: 

1. Click on 'New Site', and name it something intuitive (e.g. O2)
2. Host: transfer.rc.hms.harvard.edu 
3. Protocol: SFTP - SSH File Transfer Protocol
4. Logon Type: Normal
5. User: Username (i.e rc_trainingXX) 
6. Password: O2 password
7. Click 'Connect'

> NOTE: While using the temporary training accounts on the O2 cluster, two-factor authentication ***IS NOT*** required. However, if you explore this lesson when using your personal account, two-factor authentication ***IS*** required. 
> 
> In order to connect your laptop using FileZilla to the O2 cluster, follow steps 1-7 as outlined above. Once you have clicked 'Connect', you will receive a Duo push notification (but no indication in Filezilla) which you must approve within the short time window. Following Duo approval, FileZilla will connect to the O2 cluster.

<p align="center">
<img src="../img/filezilla_login.png" width="500">
</p>

### Filezilla Interface

You will see messages printed in the message window in the top window pane, giving a you an indication of whether or not you have successfully connected to O2. Next, if this if your first time using Filezilla we recommend that you take some time to get familiar withe the basics of the interface. This [tutorial](https://wiki.filezilla-project.org/FileZilla_Client_Tutorial_(en)) is a helpful resource.

You will see two panels in the interface. On the left hand side you will see your the files in your laptop and on the right hand side you have your home directory on O2. Both panels have a directory tree at the top and a detailed listing of the selected directory's contents underneath. In the right hand panel, navigate to where the HTML files are located on O2 `~/variant_calling/reports/`. Then decide where you would like to copy those files to on your computer and move to that directory on the left hand panel.

Once you have found the HTML output for `MultiQC` **copy it over** by double clicking it or drag it over to right hand side panel. Once you have the HTML file copied over to your computer, you can leave the `Filezilla` interface. You can then locate the HTML file on your computer and open it up in a browser. 


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
