# Variant Annotation

## Learning Objectives

- Annotate a VCF file with for functional impacts with `SnpEff`
- Differentiate between an unannotated and annotated VCF file

## Introduction

Now that we have a VCF file for our samples, we want to annotate those variants to figure out what impacts they could have on our samples. For example, we could get interested to know if a given mutation is in the 5' UTR or alters a stop codon for a given gene model. In order to do these types of analyses, we have to merge our variants with annotated transcript information and we will use `SnpEff` to do this. `SnpEff` is extremely fast and also is bundled with `SnpSift`, which we can use downstream to prioritize our variants. 

## SnpEff

[SnpEff](http://pcingola.github.io/SnpEff/) uses transcriptome annotations to predict the functional impacts of mutations in a VCF file and will modify the INFO file of the VCF file to carry the predicted functional impacts.

### Databases

Thie first step in annotating your VCF file is finding the appropriate SnpEff database to use so that the annotations are consistent with your version of the reference genome. These databases hold the gene model information that is critical for annotating a variant. `SnpEff` has tens of thousands of genome databases pre-built that you can use and they represent many publicly availible genomes. In the rare case that your genome is not represented, you can build one for to use by following the instructions [here](http://pcingola.github.io/SnpEff/se_build_db/). However, unless you are working with a rare organism or an unpublished reference genome, it is highly unlikely that you will need to do build the database yourself. To see if your genome of interest is in the SnpEff database, we first need to load the the SnpEff module:

```
module load snpEff/4.3g
```

Now that you had loaded the `SnpEff` module, you can use the following command to display all of the current genomes availible and pipe them into a `less` buffer page:

```
java -jar $SNPEFF/snpEff.jar databases | less
```

The first column is the database name and the second column in the `Genus_species` for the organism. There is also a link to where the database can be downloaded at but it can mostly be ignored as SnpEff will download it if it needs it automatically. As you can see there are tens of thousands of these. So let's exit the `less` buffer page and see which GRCh databases are availible:

```
java -jar $SNPEFF/snpEff.jar databases | grep "GRCh" 
```

We can see that this build of SnpEff has five possible GRCh databases that we can use for annotation, including one for GRCh38.p7 called GRCh38.p7.RefSeq. Now that we have found the database that we would like to use for our analysis, we can run `SnpEff`.

### Running SnpEff

Move to your `scripts` directory and create a new script named `run_SnpEff.sbatch` using `vim`.:

```
cd ~/variant_calling/scripts/
vim run_SnpEff.sbatch
```

Once inside insert mode, can can enter the shebang line, description and `SBATCH` directives:

```
#!/bin/bash
# Using SnpEff to annotate our variants

#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 1
#SBATCH --mem 24G
#SBATCH -o run_SnpEff_GRCh38.p7_%j.out
#SBATCH -e run_SnpEff_GRCh38.p7_%j.err
```

Next, we will added the line to load the `snpEff` module: 

```
module load snpEff/4.3g
```

Also, we will add our variables:

```
REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/snpeff/
SAMPLE_NAME=mutect2_syn3_normal_syn3_tumor
REFERENCE_SEQUENCE_NAME=GRCh38.p7
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=GRCh38.p7.RefSeq
DATADIR=/n/groups/hbctraining/variant_calling/reference/snpeff/data/
FILTERED_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-LCR-filt.vcf
ANNOTATED_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-LCR-filt.snpeff.vcf
```

We need to create a directory to hold out reports:

```
mkdir -p $REPORTS_DIRECTORY
```

Lastly, we need to add out `SnpEff` command:

```
java -jar -Xmx4g $SNPEFF/snpEff.jar  eff \
-dataDir $DATADIR \
-cancer \
-noLog \
-csvStats $CSV_STATS \
-s $HTML_REPORT \
$REFERENCE_DATABASE \
$FILTERED_VCF_FILE > $ANNOTATED_VCF_FILE
```

This script should look like:

```
#!/bin/bash
# Using SnpEff to annotate our variants

#SBATCH -p short
#SBATCH -t 0-1:00
#SBATCH -c 1
#SBATCH --mem 24G
#SBATCH -o run_SnpEff_GRCh38.p7_%j.out
#SBATCH -e run_SnpEff_GRCh38.p7_%j.err

module load snpEff/4.3g

REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/snpeff/
SAMPLE_NAME=mutect2_syn3_normal_syn3_tumor
REFERENCE_SEQUENCE_NAME=GRCh38.p7
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=GRCh38.p7.RefSeq
DATADIR=/n/groups/hbctraining/variant_calling/reference/snpeff/data/
FILTERED_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-LCR-filt.vcf
ANNOTATED_VCF_FILE=/n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-LCR-filt.snpeff.vcf

mkdir -p $REPORTS_DIRECTORY

java -jar -Xmx4g $SNPEFF/snpEff.jar  eff \
-dataDir $DATADIR \
-cancer \
-noLog \
-csvStats $CSV_STATS \
-s $HTML_REPORT \
$REFERENCE_DATABASE \
$FILTERED_VCF_FILE > $ANNOTATED_VCF_FILE
```

Submit this script using:

```
sbatch run_SnpEff.sbatch
```

Let's breakdown this command and discuss each argument:

`java -jar $SNPEFF/snpEff.jar  eff` `Snpeff` is a `java` packaged program, so it needs to be called with `java -jar` followed by the path where the JAR file is located on the cluster. `$SNPEFF` is just a bash variable that contains the path to JAR file. `eff` is the command within `SnpEff` to annotate variants.

`-dataDir $DATADIR` This is the path to the data directory that holds the `SnpEff` annotations

`-cancer` Performs "cancer" comparisons

`-noLog` Does not report usage statistics

`-csvStats $CSV_STATS` This produces a flat-text file with summary statistics regarding the variants annotated. (Optional)

`-s $HTML_REPORT` This creates an HTML file with summary statistics regarding the variants annotated. This HTML file is mostly just an HTML stylized version of the CSV file above. (Optional)

`$REFERENCE_DATABASE` This is the `SnpEff` database we are going to use for the annotation.

`$FILTERED_VCF_FILE` This is the input VCF file to be annotated

`> $ANNOTATED_VCF_FILE` The output of `SnpEff` will be redirected into this file.

### Output

Let's take a look at our output from `SnpEff` now to see how our VCF file has been modified to include variant annotations. Let's open up our SnpEff annotated VCF file in `less`:

```
less /n/scratch3/users/${USER:0:1}/${USER}/variant_calling/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7--LCR-filt.snpeff.vcf
```

Scroll down to right before your variants are and you will notice that `SnpEff` has inserted two lines into your VCF file. First, it has inserted:

`##SnpEffVersion="4.3g (build 2016-11-28 08:31), by Pablo Cingolani"`

This is just a bit of extra metadata on your VCF file letting you know the version of `SnpEff` that was used to generate this VCF file. While it is best practice to have your code for annotating a VCF file archived in a script  like we have done, you could encounter a situation where you've either lost this script or you are analyzing someone else's work and need to know which version of `SnpEff` they used. 

```
##SnpEffCmd="SnpEff  -cancer -csvStats /n/data1/cores/bcbio/gammerdinger/cancer-dream-syn3/replicate/vcf_files/syn3_GRCh38.p7-effects-stats.csv -s /n/data1/cores/bcbio/gammerdinger/cancer-dream-syn3/replicate/vcf_files/syn3_GRCh38.p7-effects-stats.html GRCh38.p7.RefSeq /n/data1/cores/bcbio/gammerdinger/cancer-dream-syn3/replicate/vcf_files/syn3_GRCh38.p7-raw-filt.vcf.gz "
```

This tells you the command that was used by SnpEff (without the redirect) to create this SnpEff annotated file. Thus, between the aforementioned version line and this line, you should be able to reproduce the annotation of any SnpEff annotation as long as the annotation database is part of the buit-in databases.

Let's go ahead and compare the first entry of the filtered VCF file with the annotated VCF file:

#### Filtered VCF

```
chr1    16206   .       T       A       .       map_qual;normal_artifact;strand_bias    AS_FilterStatus=map_qual,strand_bias;AS_SB_TABLE=176,11|25,0;ClippingRankSum=-5.376;DP=215;ECNT=1;FS=4.331;GERMQ=93;MBQ=35,32;MFRL=347,345;MMQ=39,22;MPOS=29;MQ=34.22;MQ0=0;MQRankSum=-3.652;NALOD=-1.359e+01;NLOD=3.05;POPAF=6.00;ReadPosRankSum=-3.033;TLOD=17.71       GT:AD:AF:DP:F1R2:F2R1:SB        0/0:100,11:0.105:111:37,4:61,7:96,4,11,0        0/1:87,14:0.140:101:39,4:46,9:80,7,14,0
```

#### Annotated VCF

```
chr1    16206   .       T       A       .       map_qual;normal_artifact;strand_bias    AS_FilterStatus=map_qual,strand_bias;AS_SB_TABLE=176,11|25,0;ClippingRankSum=-5.376;DP=215;ECNT=1;FS=4.331;GERMQ=93;MBQ=35,32;MFRL=347,345;MMQ=39,22;MPOS=29;MQ=34.22;MQ0=0;MQRankSum=-3.652;NALOD=-1.359e+01;NLOD=3.05;POPAF=6.00;ReadPosRankSum=-3.033;TLOD=17.71;ANN=A|downstream_gene_variant|MODIFIER|DDX11L1|DDX11L1|transcript|NR_046018.2|pseudogene||n.*1797T>A|||||1797|,A|downstream_gene_variant|MODIFIER|MIR6859-1|MIR6859-1|transcript|NR_106918.1|pseudogene||n.*1163A>T|||||1163|,A|downstream_gene_variant|MODIFIER|MIR6859-2|MIR6859-2|transcript|NR_107062.1|pseudogene||n.*1163A>T|||||1163|,A|downstream_gene_variant|MODIFIER|MIR6859-3|MIR6859-3|transcript|NR_107063.1|pseudogene||n.*1163A>T|||||1163|,A|downstream_gene_variant|MODIFIER|MIR6859-4|MIR6859-4|transcript|NR_128720.1|pseudogene||n.*1163A>T|||||1163|,A|intron_variant|MODIFIER|WASH7P|WASH7P|transcript|NR_024540.1|pseudogene|8/10|n.1081-259A>T||||||       GT:AD:AF:DP:F1R2:F2R1:SB        0/0:100,11:0.105:111:37,4:61,7:96,4,11,0        0/1:87,14:0.140:101:39,4:46,9:80,7,14,0
```

You can see that the only difference between these files is that `SnpEff` has appended on an `ANN` field within the `INFO` field. This `ANN` field contains lots of information such as the type of modification it is, the gene symbol of the gene modified and the accession number of the transcript modified. Because gene models overlap, you may have a single variant that alters multiple transcript models. In this case, each alteration will be separated by a `,`.

Now that we have successfully annotateed our variants, let's talk about filtering our variants in the next lesson!

[Next lesson >>](variant_filtering.md)

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
