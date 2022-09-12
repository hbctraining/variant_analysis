# File Formats

## Learning Objectives

-
-
-

## Alignment Files

### SAM

### BAM

(Briefly)

## VCF

## BED

**B**rowser **E**xtensible **D**ata (BED) is a tab-delimited file format that contains information on genomic features. A BED file's first three columns (Chromosome, Starting Position and Ending Position) are required fields. Some BED files have additional columns but these are not required.

<p align="center">
<img src="../img/bed.png" width="600">
</p>

It is important to note that BED files positioning have ***zero-based indexing***. What does this mean? There are two ways to express genomic coordinates:

- **Zero-based** is shown at the top of the image
- **One-based** is shown at the bottom of the image

<p align="center">
<img src="../img/Interbase.png" width="300">
</p>

The benefits to having a **zero-based** system is the ease of calculating distance or length of sequences. We can easily determine the length of the `ATG` sequence using the zero-based coordinates by subtracting the start from the end, whereas for one-based coordinates we would need to add one after the subtraction. Therefore, many file formats used in computation, including the BED file format, use zero-based coordinates.

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
