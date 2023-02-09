# Integrative Genomics Viewer

## Learning Objectives

- Open VCF files in IGV
- Load additional tracks in IGV
- Perform basic tasks in IGV, such as renaming tracks, resizing tracks and navigating to different regions
- Saving an IGV session


## Downloading our VCF file with `FileZilla`

The first thing we need to do to visualize our variants in IGV is to download our VCF to our local computer from the O2 cluster using `FileZilla`. Below is a refresher on the steps that we need to do in order to connect to O2 using `FileZilla`.

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
> In order to connect your computer using FileZilla to the O2 cluster, follow steps 1-7 as outlined above. Once you have clicked 'Connect', you will receive a Duo push notification (but no indication in Filezilla) which you must approve within the short time window. Following Duo approval, FileZilla will connect to the O2 cluster.

<p align="center">
<img src="../img/filezilla_login.png" width="500">
</p>

### Transferring the VCF file

Once you are connected to O2, navigate the O2 window to where your VCF files are stored:

```
/n/scratch3/users/r/rctrainingXX/variant_calling/vcf_files
```

From here you should see you your annotated VCF file:

```
mutect2_syn3_normal_syn3_tumor_GRCh38.p7-LCR-filt.pedigree_header.snpeff.dbSNP.vcf
```

Drag this file over to your your computer and places it in the appropriate directory. We are going to place ours in Desktop.

<p align="center">
<img src="../img/Move_VCF_file.png" width="1000">
</p>

>**NOTE:** If you do this transfer on your own O2 account, it will have two-factor authentication and you will need approve a Duo push notification before the download starts. 

## Viewing Variants in IGV

[Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/) is an application developed by the Broad Institute that we can use to visualize a wide variety of genomics tracks, including:

- BAM/SAM files
- BED files
- BEDGraph files
- VCF files
- Wiggle files
- BigWig files
- Many more

### Selecting a Genome

To get started, we are going to open up IGV on our computers. Once open, we are going to navigate to the reference genome that we are insterested in using fromt he dropdown menu in the top left of the IGV window. In this case, we are going to select "Human(GRCh38/hg38)":

<p align="center">
<img src="../img/Select_genome_IGV.png" width="1000">
</p>

### Load our VCF file

In order to load a file, like our VCF file, into IGV, we need to go to the top of our screen and left-click <kbd>File</kbd> -> <kbd>Load from File...</kbd>:

<p align="center">
<img src="../img/Load_from_file_IGV.png" width="600">
</p>

Then select the file we wish to open and left-click <kbd>Open</kbd>:

<p align="center">
<img src="../img/Open_file_IGV.png" width="600">
</p>

Now our IGV window should display the loaded file:

<p align="center">
<img src="../img/Loaded_file_IGV.png" width="600">
</p>

This loaded file is referred to as a *track* in IGV. We can have multipple tracks loaded as once and in the next section we will demonstrate how to load tracks provided by IGV.

### Load IGV provided tracks

For a select few genomes, like human, IGV comes with a few annotation tracks such as RefSeq genes. While Refseq genes is loaded by default in the lower panel of the IGV window, the other tracks are not automatically loaded. Let's see how to load these track by left-clicking on <kbd>File</kbd> -> <kbd>Load from Server...</kbd>:

<p align="center">
<img src="../img/Load_from_server_IGV.png" width="600">
</p>

Left-click the dropdown arrow on the left side of the window to expand all of the possible provided annotation tracks:

<p align="center">
<img src="../img/Annotation_dropdown_IGV.png" width="600">
</p>

Select the checkboxes next to <kbd>CpG Islands</kbd>, <kbd>GC %</kbd> and <kbd>Phastcons (20 way)</kbd> and left-click <kbd>Ok</kbd>:

<p align="center">
<img src="../img/Load_datasets_IGV.png" width="600">
</p>

Now our IGV window should have additional tracks and look like:

<p align="center">
<img src="../img/Loaded_tracks_IGV.png" width="800">
</p>

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
