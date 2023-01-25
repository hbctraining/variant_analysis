# ICGC-TCGA DREAM Mutation Calling Challenge Synthetic Dataset

## Learning Objectives
- Describe the dataset we will be analyzing
- Explain limitations of our dataset

## The ICGC-TCGA DREAM Mutation Calling Challenge

While it may not be obvious at this point, calling variants is not a simple task. As a result, International Cancer Genome Consortium (ICGC) and The Cancer Genome Atlas (TCGA) co-sponsered the [DREAM Mutation Calling Challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/) in order to help develop methods to more accuarately predict cancer-associated mutations from whole genome resequencing data. The idea was to allow software designers to compete to see who could most accuarately identify cancer-associateed mutations in real tumor/normal datasets.

## The Dataset

There are a few datasets that the ICGC-TCGA DREAM Mutation Calling Challenge made availible:

1) 10 Real Tumor/Normal datasets from anonymous cancer patients
2) 5 Synthetic Tumor/Normal datasets developed *in silico*

As would be expected, due to ethical standards, using the real data requires approval from ICGC and would be difficult to use in a workshop like this where we need to able to distribute the datasets to participants. Fortunately, the synthetic datasets are freely availble for use and do not require ICGC approval, so this workshop will be using a single synthetic tumor/normal sample (synthetic dataset 3) which has multiple subclones, enabling detection of lower frequency variants. 

In order to expedite our methodologies and minimize resource usage in the O2 computing cluster, we will just be using the whole exome sequencing (WES) dataset rather than a whole genome sequencing (WGS) dataset, but all of the methods that we will be using will be applicable to both WES and WGS datasets. 

[Back to Schedule](../schedule/README.md)
