# Evaluating `MultiQC` Report

## Learning Objectives

- Interpret FastQC figures for read quality
- Evaluate read alignment

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

## Inspect `MultiQC` HTML Report 

Now, we will evalute all of our FASTQC and alignments metrics at once. 



Lastly, we can evaluate our alignments in the first table:

INSERT PICTURE OF TABLE

We can see that both of our alignment rates are 99% which is great (but also a bit expected since this is a synthetic dataset). 


[Next Lesson >>>](variant_calling.md)

[Back to Schedule](../schedule/README.md)


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
