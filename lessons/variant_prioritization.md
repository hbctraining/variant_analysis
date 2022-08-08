# Variant Priorization

## Learning Objectives

- Filter records in a VCF file various effects/impacts
- Extract fields of interest from VCF file
- Create the input necessacary to create an OncoPrint

Now that we have annotated and filtered our variants, we are likely interested in subsetting our variants to find those of most interest to our study. Perhaps we are interested in finding variants that substantially disrupt a transcript, such as a variant causing a premature stop codon, or just find all of the missense mutations in the sample. `SnpSift` is part of the `SnpEff` suite and it is built explicitly for this purpose.

Let's start by discussing some of the ways you can filter your data with SnpSift.

## Filter

### Fields

First, you can filter your SnpEff annotated VCF file based upon the first seven fields of the VCF file:

- **CHROM**
- **POS**
- **ID**
- **REF** 
- **ALT** 
- **QUAL**
- **FILTER**

Let's go ahead and do our first `SnpSift` command to extract variants, but before we do we will need to load the `SnpEff` module:

```
module load snpEff/4.3g
```

Now that we have loaded the SnpEff module we can utilize the SnpSift to find all of the variants on `chr1` and pipe the output into `less`:

```
java -jar $SNPEFF/SnpSift.jar filter "( CHROM = 'chr1' )" syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | less
```

Let's break down the syntax a bit:

- `java -jar $SNPEFF/SnpSift.jar` This calls the `SnpSift` program.

- `filter` There are several useful commands within `SnpSift`, but `filter` is the one we will be using to extract variants meeting out criteria

- `"( CHROM = 'chr1' )"` This is the syntax needed to extract variants on `chr1`. The left side of the equals sign corresponds to the VCF field you wish to filter by and the right side if the string you would like to match.

- `syn3_GRCh38.p7-LCR-filt.snpeff.vcf` The input VCF file. Importantly, this need to go at the end.

- `| less` Piping the output into a `less` buffer page for inspection.

### Multiple Filters

It is likely that you could want to filter on multiple criteria. You can do that by separating the filter criteria with either and `&` (and) or `|` (or).

For example, let's consider a case where you want to filter your filter for any variant on `chr1` ***OR*** `chr2`. That might look like:

```
java -jar $SNPEFF/SnpSift.jar filter "( CHROM = 'chr1' ) | ( CHROM = 'chr2' )" syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | less
```

Note the `"( CHROM = 'chr1' ) | ( CHROM = 'chr2' )"` syntax allows us to filter for `chr1` or `chr2` by using the `|` to separate our criteria.

Alternatively, we could be interested in variants on `chr1` between positions `1000000` and `2000000`. It would look like:

```
java -jar $SNPEFF/SnpSift.jar filter "( CHROM = 'chr1' ) & ( POS > 1000000 ) & ( POS < 2000000 )" syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | less
```

### INFO Field

`SnpSift` also allows the user to filter based upon the `INFO` field. This is particularly helpful since `SnpEff`'s annotations are placed into the `INFO` field. There are many `INFO` field filters that one can apply but we will discuss some of the more popular ones. Filtering the `INFO` will importantly always begin with an `ANN` in the filtering criteria.

#### Gene

If you are interested in all of the variants corresponding to a single gene of interest, you can filter by the gene name in this case `CPSF3L`:

```
java -jar $SNPEFF/SnpSift.jar filter "( ANN[*].GENE = 'CPSF3L' )" syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | less
```

To filter by the gene name you will need `"( ANN[*].GENE = 'INSERT_GENE_NAME' )"`. 

> NOTE: When handling multiple valued fields (i.e. fields with commas), `SnpSift` uses a 0-based index to describing those elements. In the above example, the output have multiple annotations for `CPSF3L`, but if we wanted to specify that the first one needed to be `CPSF3L`, then we would need to filter by `"( ANN[0].GENE = 'CPSF3L' )"`. The `*` tell `SnpSift` to extract the record if "any" annotations corresponds to `CPSF3L`. For most cases, you will want to use the `*`, but you should know why it is there.

#### Effects

It is also quite common to want to filter your output by the effects the variants have on the annotated gene models. The syntax for this is quite similar to the example for genes:

```
java -jar $SNPEFF/SnpSift.jar filter "( ANN[*].EFFECT has 'missense_variant' )" syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | less
```

To filter by a variant effect, the filter syntax is "( ANN[\*].EFFECT has 'VARIANT_EFFECT' )"

> Note: Importantly, notice the use of `has` instead of `=` here. Sometimes effects will contain mutliple effects such as `missense_variant&splice_donor_variant`. Using `ANN[\*].EFFECT = missense_variant` here ***WILL NOT*** return this line, because the line is not equal to `missense_variant`, however `ANN[\*].EFFECT has missense_variant` ***WILL*** return this line. Oftentimes for effects, one would be interested in the `has` query as opposed to the `=` one.

There are many different variant effects and some of the more common ones are:

- `missense_variant` for missense/non-synonymous variants 

- `frameshift_variant` for frameshift variants

- `stop_gain` for nonsense variants

- `stop_lost` for variants that lose a stop codon

- `start_gain` for variants that gain a start codon

- `start_lost` for variants that lose a start codon

- `synonymous_variant` for synonymous/silent variants

- `splice_donor_variant` for a variant in the splice donor site

- `splice_acceptor_variant` for a variant in the splice acceptor site

- `5_prime_UTR_variant` for a variant in the 5' untranslated region 

- `3_prime_UTR_variant` for a variant in the 3' untranslated region 

Many more effects can be found [here](https://pcingola.github.io/SnpEff/se_inputoutput/#effect-prediction-details).

#### Impacts

`SnpEff` also predicts the deleterious nature of a variant by binning it into one of several categories:

- `HIGH` These are variants that will almost certainly have a deleterious impact on the transcript. Examples of this would be the loss or gain of a stop codon or a frameshift mutation. 

- `MODERATE` These are variants where the impact may have a deleterious impact on the transcript. Examples of this would be missense/non-synonymous variants and in-frame deletions/insertions.

- `LOW` These are variants that are unlikely to have a deleterious impact on the transcript. Examples of this would be silent/synonymous variants and alterations between different stop codons.

- `MODIFER` These variants are typically in non-coding regions and their impacts are difficult to assertain. 

More information on these categories can be found [here](https://pcingola.github.io/SnpEff/se_inputoutput/#impact-prediction) and a complete listing of the categories for each effect can be found [here]((https://pcingola.github.io/SnpEff/se_inputoutput/#effect-prediction-details). 

Let's go ahead and filter out all of our `HIGH` impact muations:

```
java -jar $SNPEFF/SnpSift.jar filter "( ANN[*].IMPACT has 'HIGH' )"  syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | less
```

> ***Note:*** Similarly to `EFFECT`, oftentimes you will want to use `has` rather than `=`.

#### Other ANN fields

In addition to `GENE`, `EFFECT` and `IMPACT`, there are a whole host of other `ANN` fields. Some of the others we may come across are:

- `TRID` - Transcript ID or NCBI accesssion number
- `HGVS_P` - The alteration in protein notation
- `HGVS_C` - The alteration in DNA notation

A full list of `ANN` fields can be found [here](http://pcingola.github.io/SnpEff/ss_filter/#snpeff-ann-fields).

## vcfEffOnePerLine

A useful tool within the `SnpSift` toolkit is the `perl` script named `vcfEffOnePerLine.pl`. This script allows the user to separate each effect onto its own line instead of having them lumped into a single line. In order to utilize this script we need to pipe the output of our `filter` command into `$SNPEFF/scripts/vcfEffOnePerLine.pl`. We can use it on our previous example to demonstrate:

```
java -jar $SNPEFF/SnpSift.jar filter "( ANN[*].IMPACT has 'HIGH' )"  syn3_GRCh38.p7-LCR-filt.snpeff.vcf  | $SNPEFF/scripts/vcfEffOnePerLine.pl | less
```

Now, we can see that each variant has a separate entry depending on its effect.

```
chrX    153285023       .       C       A       .       strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,strand_bias;AS_SB_TABLE=12,45|0,2;ClippingRankSum=1.288;DP=61;ECNT=1;FS=0.000;GERMQ=93;MBQ=26,35;MFRL=339,349;MMQ=60,60;MPOS=15;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.45;NLOD=8.07;POPAF=6.00;ReadPosRankSum=-1.027;TLOD=3.78;LOF=(IRAK1|IRAK1|3|1.00);NMD=(IRAK1|IRAK1|3|1.00);ANN=A|stop_gained|HIGH|IRAK1|IRAK1|transcript|NM_001569.3|protein_coding|2/14|c.163G>T|p.Glu55*|242/3571|163/2139|55/712||       GT:AD:AF:DP:F1R2:F2R1:SB        0/0:28,0:0.035:28:15,0:12,0:7,21,0,0    0/1:29,2:0.091:31:17,0:10,2:5,24,0,2
chrX    153285023       .       C       A       .       strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,strand_bias;AS_SB_TABLE=12,45|0,2;ClippingRankSum=1.288;DP=61;ECNT=1;FS=0.000;GERMQ=93;MBQ=26,35;MFRL=339,349;MMQ=60,60;MPOS=15;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.45;NLOD=8.07;POPAF=6.00;ReadPosRankSum=-1.027;TLOD=3.78;LOF=(IRAK1|IRAK1|3|1.00);NMD=(IRAK1|IRAK1|3|1.00);ANN=A|stop_gained|HIGH|IRAK1|IRAK1|transcript|NM_001025242.1|protein_coding|2/14|c.163G>T|p.Glu55*|242/3481|163/2049|55/682||    GT:AD:AF:DP:F1R2:F2R1:SB        0/0:28,0:0.035:28:15,0:12,0:7,21,0,0    0/1:29,2:0.091:31:17,0:10,2:5,24,0,2
chrX    153285023       .       C       A       .       strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,strand_bias;AS_SB_TABLE=12,45|0,2;ClippingRankSum=1.288;DP=61;ECNT=1;FS=0.000;GERMQ=93;MBQ=26,35;MFRL=339,349;MMQ=60,60;MPOS=15;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.45;NLOD=8.07;POPAF=6.00;ReadPosRankSum=-1.027;TLOD=3.78;LOF=(IRAK1|IRAK1|3|1.00);NMD=(IRAK1|IRAK1|3|1.00);ANN=A|stop_gained|HIGH|IRAK1|IRAK1|transcript|NM_001025243.1|protein_coding|2/13|c.163G>T|p.Glu55*|242/3334|163/1902|55/633||    GT:AD:AF:DP:F1R2:F2R1:SB        0/0:28,0:0.035:28:15,0:12,0:7,21,0,0    0/1:29,2:0.091:31:17,0:10,2:5,24,0,2
chrX    153285023       .       C       A       .       strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,strand_bias;AS_SB_TABLE=12,45|0,2;ClippingRankSum=1.288;DP=61;ECNT=1;FS=0.000;GERMQ=93;MBQ=26,35;MFRL=339,349;MMQ=60,60;MPOS=15;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.45;NLOD=8.07;POPAF=6.00;ReadPosRankSum=-1.027;TLOD=3.78;LOF=(IRAK1|IRAK1|3|1.00);NMD=(IRAK1|IRAK1|3|1.00);ANN=A|downstream_gene_variant|MODIFIER|MIR718|MIR718|transcript|NR_031757.1|pseudogene||n.*348G>T|||||348|       GT:AD:AF:DP:F1R2:F2R1:SB        0/0:28,0:0.035:28:15,0:12,0:7,21,0,0    0/1:29,2:0.091:31:17,0:10,2:5,24,0,2
chrX    153285023       .       C       A       .       strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,strand_bias;AS_SB_TABLE=12,45|0,2;ClippingRankSum=1.288;DP=61;ECNT=1;FS=0.000;GERMQ=93;MBQ=26,35;MFRL=339,349;MMQ=60,60;MPOS=15;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.45;NLOD=8.07;POPAF=6.00;ReadPosRankSum=-1.027;TLOD=3.78;LOF=(IRAK1|IRAK1|3|1.00);NMD=(IRAK1|IRAK1|3|1.00);ANN=A|downstream_gene_variant|MODIFIER|MECP2|MECP2|transcript|NM_004992.3|protein_coding||c.*10795G>T|||||2241|  GT:AD:AF:DP:F1R2:F2R1:SB        0/0:28,0:0.035:28:15,0:12,0:7,21,0,0    0/1:29,2:0.091:31:17,0:10,2:5,24,0,2
```

Which was previously:

```
chrX    153285023       .       C       A       .       strand_bias;weak_evidence       AS_FilterStatus=weak_evidence,strand_bias;AS_SB_TABLE=12,45|0,2;ClippingRankSum=1.288;DP=61;ECNT=1;FS=0.000;GERMQ=93;MBQ=26,35;MFRL=339,349;MMQ=60,60;MPOS=15;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.45;NLOD=8.07;POPAF=6.00;ReadPosRankSum=-1.027;TLOD=3.78;ANN=A|stop_gained|HIGH|IRAK1|IRAK1|transcript|NM_001569.3|protein_coding|2/14|c.163G>T|p.Glu55*|242/3571|163/2139|55/712||,A|stop_gained|HIGH|IRAK1|IRAK1|transcript|NM_001025242.1|protein_coding|2/14|c.163G>T|p.Glu55*|242/3481|163/2049|55/682||,A|stop_gained|HIGH|IRAK1|IRAK1|transcript|NM_001025243.1|protein_coding|2/13|c.163G>T|p.Glu55*|242/3334|163/1902|55/633||,A|downstream_gene_variant|MODIFIER|MIR718|MIR718|transcript|NR_031757.1|pseudogene||n.*348G>T|||||348|,A|downstream_gene_variant|MODIFIER|MECP2|MECP2|transcript|NM_004992.3|protein_coding||c.*10795G>T|||||2241|;LOF=(IRAK1|IRAK1|3|1.00);NMD=(IRAK1|IRAK1|3|1.00)        GT:AD:AF:DP:F1R2:F2R1:SB        0/0:28,0:0.035:28:15,0:12,0:7,21,0,0    0/1:29,2:0.091:31:17,0:10,2:5,24,0,2
```

This step is particularly helpful for cleaning up the files for use in the next step `extractFields`.

## extractFields

Lastly, we have another extremely useful feature of `SnpSift` and that is the `extractFields` command. This allows us to parse the VCF file and print only the fields we are interested in. 

If we wanted to parse out the missense mutations, create a single line for each effect, then extract the fields for chromosome, position, gene ID as well as the alteration in terms of protein and DNA space and also the effect.

```
java -jar $SNPEFF/SnpSift.jar \
filter \
"( ANN[*].EFFECT has 'missense_variant' )"  vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf| \
$SNPEFF/scripts/vcfEffOnePerLine.pl | \
java -jar $SNPEFF/SnpSift.jar \
extractFields \
-  \
"CHROM" "POS" "ANN[*].GENE" "ANN[*].TRID" "EFF[*].HGVS_P" "ANN[*].HGVS_C" "ANN[*].EFFECT" | less
```

Note the use of `-` within the `extractFields` fields command. `-` is very commonly used to define the input as coming from standard input, or in other words, the input is being piped into the command. 

`"CHROM" "POS" "ANN[*].GENE" "ANN[*].TRID" "EFF[*].HGVS_P" "ANN[*].HGVS_C" "ANN[*].EFFECT"` is defining the fields that we would like to filter. 

Notice however, that some of the fields don't correspond to a `missense_variant`. This is because when we initially extracted sites we filtered sites where at least one effect was a `missense_variant`, then we separated the variants into separate lines before extracting our fields of interest. At this point if we wanted to remove those non-missense variant lines, we can pipe the output into a simple `grep` command:

```
java -jar $SNPEFF/SnpSift.jar \
filter \
"( ANN[*].EFFECT has 'missense_variant' )"  vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf| \
$SNPEFF/scripts/vcfEffOnePerLine.pl | \
java -jar $SNPEFF/SnpSift.jar \
extractFields \
-  \
"CHROM" "POS" "ANN[*].GENE" "ANN[*].TRID" "EFF[*].HGVS_P" "ANN[*].HGVS_C" | \
grep 'missense_variant' | less
```

## Output for OncoPrint
```
java -jar $SNPEFF/SnpSift.jar \
filter \
"( ANN[*].EFFECT has 'missense_variant' )"  vcf_files/syn3_GRCh38.p7-LCR-filt.snpeff.vcf| \
$SNPEFF/scripts/vcfEffOnePerLine.pl | \
java -jar $SNPEFF/SnpSift.jar \
extractFields \
-  \
"ANN[*].GENE" "EFF[*].HGVS_P" "ANN[*].EFFECT"  | \
sort | \
uniq | \
awk '$3 == "missense_variant"' | \
sed 's/missense_variant/MISSENSE/g' | \
awk '{print "sample_1","\t",$0}' > vcf_files/syn3_GRCh38.p7.missense_variants.txt
```

[Next Lesson >>>](IGV.md)

[Back to Schedule](../schedule/README.md)


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
