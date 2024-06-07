# Oncoprint Creation

## Learning Objectives

- Create oncoprint on cBioPortal using data analyzed in class 

## Oncoprints



## Wrangling our data for an Oncoprint

```
cd ~/variant_calling/scripts/
vim VCF_to_oncoprint.sh
```

```
#!/bin/bash
# This script was wrttien by the Training Team at the Harvard Chan Bioinformatics Core on June 6th, 2024  as part of training materials for the Introduction to Variant Analysis workshop.
# USAGE: sh VCF_to_oncoprint.sh <INPUT_VCF_FILE> <SAMPLE_NAME>

# Assign variable for input and output
INPUT_VCF_FILE=$1
SAMPLE_NAME=$2
OUTPUT_FILE=${INPUT_VCF_FILE%vcf}oncoprint.txt

# Load SnpEff
module load snpEff/4.3g

java -jar $SNPEFF/SnpSift.jar filter \
  -noLog \
  "( ANN[*].EFFECT has 'missense_variant' ) | ( ANN[*].EFFECT has 'stop_gain' ) | ( ANN[*].EFFECT has 'frameshift_variant' ) | ( ANN[*].EFFECT has 'inframe_insertion' ) | ( ANN[*].EFFECT has 'inframe_deletion' ) " \
  $INPUT_VCF_FILE  | \
  $SNPEFF/scripts/vcfEffOnePerLine.pl | \
  java -jar $SNPEFF/SnpSift.jar extractFields \
  - \
  "ANN[*].GENE"  "ANN[*].HGVS_P"  "ANN[*].EFFECT" |  \
  awk -v sample_name=$SAMPLE_NAME  'NR>1 {print sample_name,$1,$2,$3}' | \
  grep -E "missense_variant|stop_gain|frameshift_variant|inframe_insertion|inframe_deletion" | \
  awk '{sub(/p\./, "", $3); print}' | \
  awk '{sub(/Phe/, "F", $3); print}' | \
  awk '{sub(/Leu/, "L", $3); print}' | \
  awk '{sub(/Ile/, "I", $3); print}' | \
  awk '{sub(/Met/, "M", $3); print}' | \
  awk '{sub(/Val/, "V", $3); print}' | \
  awk '{sub(/Ser/, "S", $3); print}' | \
  awk '{sub(/Pro/, "P", $3); print}' | \
  awk '{sub(/Thr/, "T", $3); print}' | \
  awk '{sub(/Ala/, "A", $3); print}' | \
  awk '{sub(/Tyr/, "Y", $3); print}' | \
  awk '{sub(/His/, "H", $3); print}' | \
  awk '{sub(/Gln/, "Q", $3); print}' | \
  awk '{sub(/Asn/, "N", $3); print}' | \
  awk '{sub(/Lys/, "K", $3); print}' | \
  awk '{sub(/Asp/, "D", $3); print}' | \
  awk '{sub(/Glu/, "E", $3); print}' | \
  awk '{sub(/Cys/, "C", $3); print}' | \
  awk '{sub(/Trp/, "W", $3); print}' | \
  awk '{sub(/Arg/, "R", $3); print}' | \
  awk '{sub(/Gly/, "G", $3); print}' | \
  awk 'BEGIN { OFS="\t" } {print $1,$2,$3,$4}' | \
  sed 's/frameshift_variant.*/TRUNC/g' | \
  sed 's/stop_gain.*/TRUNC/g' | \
  sed 's/inframe_insertion.*/INFRAME/g' | \
  sed 's/inframe_deletion.*/INFRAME/g' | \
  sed 's/missense_variant.*/MISSENSE/g' > $OUTPUT_FILE
```

```
sh VCF_to_oncoprint.sh /n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/mutect2_syn3_normal_syn3_tumor_hg38-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf  syn3
```



***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
