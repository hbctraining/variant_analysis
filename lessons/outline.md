# Starting position

Raw VCF file with tumor and normal samples

# Filtering 

Tool: VCFtools 

Manual: http://vcftools.sourceforge.net/man_latest.html

What filtering settings do we want to do? Depth? Number of alleles? Remove entries without "PASS" filter status? Remove indels?

```
vcftools --vcf raw_variants.vcf --out filtered_variants.vcf
```

Options to consider:

```
--remove-indels
--remove-filtered-all
--max-alleles <string>
```

Filter for depth

Average across samples:

```
--min-meanDP <float> 
--max-meanDP <float>
```

OR

Across samples:

```
--minDP <float> 
--maxDP <float>
```
# Variant Annotation

Tool: SNPEff

Manual: http://pcingola.github.io/SnpEff/se_running/

```
java -Xmx8g -jar snpEff.jar Reference_genome filtered_variants.vcf > filtered_variants.ann.vcf
```

Option to consider:

```
-cancer
-cancerSamples <file>
``
The cancer samples file is a tab-delimited file with original<tab>derived.


# Variant Prioritization

Tool: SNPSift

Manual: http://pcingola.github.io/SnpEff/ss_filter/

```
java -jar SnpSift.jar filter "ANN[*].EFFECT == 'missense_variant'" filtered_variants.ann.vcf > missense_variants.vcf

java -jar SnpSift.jar filter  "ANN[*].IMPACT == 'HIGH'" filtered_variants.ann.vcf > high_impact_variants.vcf
```

We could save the filtering step until here as well. Not sure if we have thoughts either way on this.

**My biggest hang up is how to utilize our normal samples to. Not sure the bast way to do this.**

Optional or fork towards IGV:

Tool: dbNSFP

http://pcingola.github.io/SnpEff/ss_dbnsfp/

java -jar SnpSift.jar dbnsfp -v missense_variants.vcf > missense_variants.scored.vcf

I'm not sure exactly how this output looks or much about this tool. Input is welcome.


# Oncoprint Creation

Once we have the VCF file the way we want it. We can create the input file for creating an Oncoprint.

Oncoprints input format needs to be:

```
Sample  Gene  HGVS_P  Type
```

The Extract fields function in SNPSift will be useful for this step. 

Pull out a set of genes using SNPSift.

