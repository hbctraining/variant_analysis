# Starting position

Raw VCF file with tumor and normal samples

# Filtering 

Tool: FilterMutectCalls

Manual: https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls

```
 gatk FilterMutectCalls -V somatic.vcf -R reference_genome.fasta -O filtered.vcf.gz
```

Options to consider:

```
--max-alt-allele-count <integer>
```

Next?

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
``` 

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
 
Alternative to dbNSFP 
 
Tool: open-cravat

Website: https://opencravat.org/resources.html 
 
This seems like a web interface but supposedly seems give similar output as dbNSFP (maybe?).
 
# Oncoprint Creation

Once we have the VCF file the way we want it. We can create the input file for creating an Oncoprint.

Oncoprints input format needs to be:

```
Sample  Gene  HGVS_P  Type
```

The Extract fields function in SNPSift will be useful for this step. 

Pull out a set of genes using SNPSift.

Plug text file into OncoPrinter

https://www.cbioportal.org/oncoprinter
  
  
# Visualize in IGV

Pull up several tumor/normal bam files in IGV. Subset them down to like 2-3 genes so that participants can just download them.
  
# Mutational Signature Analysis
 
 ***Needs R Background***
 
 ***Also, not designed for WES, but rather WGS***
 
 https://github.com/Nik-Zainal-Group/signature.tools.lib
  
# Copy Number Variant Analysis
 
 I think we nixed this from the initial course, but Sergey linked a great paper comparing six different CNV detection algorrithms, so this could be a starting place if we want to do this. I am adding the paper here for posterity.
 
 https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07686-z
