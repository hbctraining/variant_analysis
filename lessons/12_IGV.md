# Integrative Genomics Viewer

## Learning Objectives

- Open VCF files in IGV
- Load additional tracks in IGV
- Perform basic tasks in IGV, such as renaming tracks, resizing tracks and navigating to different regions
- Saving an IGV session and opening a saved IGV session


## Downloading our VCF file with FileZilla

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
> In order to connect your computer using `FileZilla` to the O2 cluster, follow steps 1-7 as outlined above. Once you have clicked 'Connect', you will receive a Duo push notification (but no indication in Filezilla) which you must approve within the short time window. Following Duo approval, `FileZilla` will connect to the O2 cluster.

<p align="center">
<img src="../img/filezilla_login.png" width="500">
</p>

### Transferring the VCF file

Once you are connected to O2, navigate the O2 window to where your VCF files are stored:

```
/n/scratch/users/r/rctrainingXX/variant_calling/vcf_files
```

From here you should see you your annotated VCF file:

```
mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf
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

To get started, we are going to open up IGV on our computers. Once open, we are going to navigate to the reference genome that we are interested in using from the dropdown menu in the top left of the IGV window. In this case, we are going to select "Human(GRCh38/hg38)":

<p align="center">
<img src="../img/Select_genome_IGV.png" width="1000">
</p>

### Loading Our VCF File

In order to load a file, like our VCF file, into IGV, we need to go to the top of our screen and left-click <kbd>File</kbd> &#8594; <kbd>Load from File...</kbd>:

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

This loaded file is referred to as a *track* in IGV. We can have multiple tracks loaded as once and in the next section we will demonstrate how to load tracks provided by IGV.

### Load IGV Provided Tracks

For a select few genomes, like human, IGV comes with a few annotation tracks such as RefSeq genes. While Refseq genes is loaded by default in the lower panel of the IGV window, the other tracks are not automatically loaded. Let's see how to load these track by left-clicking on <kbd>File</kbd> &#8594; <kbd>Load from Server...</kbd>:

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

### Navigating IGV

Now that we have our tracks of interest loaded, we likely want to navigate them to view our variants.

#### Zooming in/out on Regions

The first way that we can zoom in on a region in IGV is to left-click and hold while dragging over the region we are interested in. This can be iteratively done as one narrows down the region that they are interested in viewing.

<p align="center">
<img src="../img/Zooming_click_and_drag_IGV.png" width="800">
</p>

Once we zoom, our IGV window might look like:

<p align="center">
<img src="../img/Zoomed_window_IGV.png" width="800">
</p>

We can zoom out and in using the <kbd>-</kbd> and <kbd>+</kbd>, respectively, on the right of the top bar:

<p align="center">
<img src="../img/Zoom_buttons_IGV.png" width="400">
</p>

>**NOTE:** If you can't see the <kbd>-</kbd> and/or <kbd>+</kbd> sign on the right side of the screen, you may need to widen your window. You can see in our image above this one that we couldn't see the <kbd>+</kbd> sign, but we made the IGV window wider and then we could see it.

#### Jumping to Regions

We can jump to a given region in the genome using the following syntax in the middle of the top bar:

```
chr<Chromosome_Number>:<Start_position>:<End_position>
```

Then left-clicking <kbd>Go</kbd>.

<p align="center">
<img src="../img/Coordinate_jump_IGV.png" width="400">
</p>

Alternatively, if there is a gene we are particularly interested in going to, we can also enter the gene's name in this same box and left-click <kbd>Go</kbd>:

<p align="center">
<img src="../img/Gene_jump_IGV.png" width="400">
</p>

If we zoom in enough we can get down to seeing individual variants. Underneath the top part of the variant track, we can see the tracks for the variant in the normal and tumor samples. The intensity of the variant marks on those tracks is proportional to the frequency of the variant in the sample:

<p align="center">
<img src="../img/Zoomed_variant_IGV.png" width="800">
</p>


### Modifying Tracks

Now that we have an idea of how we can navigate in the IGV window, we want to learn how to modify the tracks. IGV provides lots of ways in which you can modify the tracks.

#### Resizing Tracks

First, we might want to be interested in resizing the tracks. Currently, we have a lot of whitespace that we might want to eliminate. The easiest way eliminate lots of this whitespace is to hit the button on the top bar to "Resize tracks to fit in window":

<p align="center">
<img src="../img/Resize_button_IGV.png" width="400">
</p>

Now the tracks should fit into the window space a but better and look like:

<p align="center">
<img src="../img/Resized_to_window_IGV.png" width="800">
</p>

> **NOTE:** If you have too few tracks, like we have here, they likely won't encompass all of the white space or if you have too many tracks, IGV might not be able to cram them all into the window and you will have a scroll bar on the side.

We can also manually adjust the track height by right-clicking on the track that we want to alter in size and the left-clicking "Change Track Height.."

<p align="center">
<img src="../img/Change_height_IGV.png" width="600">
</p>

A window should pop-up and allow you to modify the height. Once you have selected a height you can left-click <kbd>OK</kbd> and the track will be resized.

<p align="center">
<img src="../img/Adjust_track_height_manual_IGV.png" width="400">
</p>


#### Modifying Track Data Range

We can also adjust the data range that we want displayed in the track. Similarly to resizing the track height, we start by right-clicking the track we want to adjust, but this time we will left-click "Set Data Range..."

<p align="center">
<img src="../img/Set_data_range_IGV.png" width="400">
</p>

A window should pop-up and allow you to select the minimum, maximum and whether you would like the data to be log-scaled. Once you have selected the parameters you want, you can left-click <kbd>OK</kbd>:

<p align="center">
<img src="../img/Set_data_range_manual_IGV.png" width="400">
</p>

#### Renaming Tracks

When tracks are loaded from files, they are named by the filename and sometimes these names can be long and unwieldy within IGV. Therefore, we are oftentimes interested in changing the track name to something that is more easy to understand. To change the track name, we need to right-click on the track and then left-click "Rename Track...":

<p align="center">
<img src="../img/Rename_track_IGV.png" width="600">
</p>

A window should pop-up and allow you to type the desired name of the track. Once you have typed the desired name, left-click <kbd>OK</kbd>:

<p align="center">
<img src="../img/Rename_track_manual_IGV.png" width="400">
</p>

#### Changing Track Color

Many tracks load into IGV as blue by default, but you do have options for which color you'd like the tracks to be. In order to change the color of a track, right click on the track and left-click "Change Track Color...":

<p align="center">
<img src="../img/Change_track_color_IGV.png" width="400">
</p>

A window will pop-up on the default "Swatches" tab on the top. You can pick from a wide array for pre-selected colors here. If you find one you like, left-click the color then click <kbd>OK</kbd>:

<p align="center">
<img src="../img/Swatches_color_IGV.png" width="600">
</p>

However, you may want finer control over your color selection and you can use some of the other tabs to do this. The "RGB" tab allows you to define the level of red, green and blue you want in the color. Of particular note, it also allows you to place the hexidemical code for the color you want in the "Color Code" text box. For instance, this could be of interest if you are trying to keep consistent colors from other figures where you defined a hexidecimal code for a given dataset. Once you have selected a color that you like, you can left-click <kbd>OK</kbd>:

<p align="center">
<img src="../img/RGB_color_IGV.png" width="600">
</p>

Now that we've changed a few features in our IGV window it should now look something like this:

<p align="center">
<img src="../img/Changed_color_IGV.png" width="800">
</p>

#### Changing Type of Graph

Different types of data might be visualized better in different formats. For example, we might think that "GC %" is better visualized as a line rather than as a barplot. In order to change this barplot to a line plot, we need to right-click on the track and then left-click on "Line Plot":

<p align="center">
<img src="../img/Line_plot_IGV.png" width="400">
</p>

Depending on your datatype, different types of plots might be more appropriate than others.

#### Remove Track

Lastly, you may want to remove a track from your IGV window. In order to remove a track, right-click on the track and then left-click on "Remove Track":

<p align="center">
<img src="../img/Remove_track_IGV.png" width="400">
</p>

### Saving and Loading IGV Sessions

Oftentimes, you'll want to save your IGV session or load up an IGV session that you've been previously working on. Below we will describe how to save and load IGV sessions.

#### Saving an IGV Session

Now that you have edited your tracks to get them just the way you want you them, you might want to save the IGV session so that you can easily reload it for when you want to revisit it. To save your IGV session, go to the top bar and left-click <kbd>File</kbd> &#8594; <kbd>Save Session...</kbd>:

<p align="center">
<img src="../img/Save_session_IGV.png" width="400">
</p>

Select a name and location for the IGV session to be saved under and left-click <kbd>Save</kbd>. It will then be saved as an XML file.

<p align="center">
<img src="../img/Saving_session_IGV.png" width="400">
</p>

Let's go ahead and close our IGV session now.

#### Loading an IGV Session

If we now open IGV back up, we will notice that it provides a fresh session. If we want to load a previous IGV session we will need to load it. To load an IGV session, go to the top bar and left-click <kbd>File</kbd> &#8594; <kbd>Open Session...</kbd>:

<p align="center">
<img src="../img/Load_session_IGV.png" width="400">
</p>

A window should pop-up and let you select the file you would like to load. Left-click the file you would like to load and then left-click <kbd>Open</kbd>:

<p align="center">
<img src="../img/Loading_session_IGV.png" width="400">
</p>

We can now see that we have loaded our previous IGV session! ***It is VERY IMPORTANT that if you move files that were loaded into IGV into a different location on your computer, IGV will not be able to find them and therefore not load your saved session!***

## Exercises

**1)** Download and load the SnpSift file that we created with "High Impact" mutations

**2)** Load the IGV provided "Variation and Repeats" track to your IGV session

**3)** Change the height of the CpG Islands track to 60

**4)** Navigate to your favorite gene. Do you see any high-impact variants there?

**5)** Find a high impact variant on Chromosome 4

**6)** Rename the "Refseq Genes" track to "Genes"

**7)** Save the IGV session as "Improved_IGV_session.xml"

**8)** *Bonus Challenge* Change the type of plot for "GC%" from a line plot to a heatmap.

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
