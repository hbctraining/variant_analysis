# Oncoprint Creation

## Learning Objectives

-
-
-

Filter out TCGA file from cBioPortal data_mutations.txt

*Note: this is from GRCh37*

```
awk '{print $18,"\t",$1,"\t",$39,"\t"$9}' brca_tcga_pan_can_atlas_2018/data_mutations.txt | \
sed '1d' > brca_filtered.txt
```

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
