## Different types of syntax for Picard

The syntax that `Picard` uses is quite particular and the syntax shown in the documentation is not always consistent. There are **two main ways** for providing input for `Picard`:

1. The **Traditional Syntax**
   * When providing inputs in this fashion you will provide the option followed immediately without any whitespace by an `=` sign followed once again immediately without whitespace by the argument for that option.
   * It will work and produce a valid output, but your error file might contain a warning

2. The **New (Barclay) Syntax**
   * Picard migrated to a new syntax several years back and this is the syntax we are demonstrating in this workshop. 
   * It uses a double hyphen (`--`) followed by the long-form name for the option _followed by whitespace_ followed by the argument for that option. 
 

For either syntax, you **can also provide a single/couple letter abbreviation**s for an option (as is the case with many other software). However, we have elected to write out the full names of options wherever possible to increase readabilty of the code. 

Commands written in either syntax are equally valid and produce the same output. However, we think you should be aware of the differences as Picard's documentation will sometimes use the New (Barclay) syntax is one place and the Traditional syntax in another place, sometimes even on the same page! 

***

**Exercise**

**1.** We are comparing our `SortSam` command with our colleague's command. Is there anything wrong with their syntax? Why or why not?

**Our syntax**
```
java -jar $PICARD/picard.jar SortSam \
--INPUT $SAM_FILE \
--OUTPUT $QUERY_SORTED_BAM_FILE \
--SORT_ORDER queryname
```

**Our colleague's syntax**
```
java -jar $PICARD/picard.jar SortSam \
I=$SAM_FILE \
O=$QUERY_SORTED_BAM_FILE \
SO=queryname
```

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
