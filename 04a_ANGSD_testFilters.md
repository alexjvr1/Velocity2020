# Test filters and overlap in loci across datasets

There's a big dropout in the number of loci kept in the datasets after applying ANGSD filters. 

This is particularly true for the modern data where only 2% of the data are kept! (although initial counts include invariant sites)

In the museum data ~63% of the loci are retained using the initial filters (inlcuding keeping SNPs with p-value 0.05)

### Aim:

1) How many loci are filtered with each option?

2) How does this affect the intersect between datasets? 

3) How much time does it add to the analysis? 


### Method: 

I'll use only the largest contig for the analysis. 
I'll start with the basic filters then add one each time. 

## 1. Basic filter set

Filters we have to include: 

-remove_bads 1 : remove reads with 255 flag (not primary, failure and duplicate reads) (1=default)

-uniqueOnly 1 : remove reads with multiple best hits

-minMapQ 20 : PHRED 20. This should already be in place during the mapping.

-minQ 20 : PHRED 20 for individual base score.

-only_proper_pairs 0 : NBNB THIS flag is changed to 0 because some of the reads in my final mus files are not properly paired! ***OLD FILTER include only properly paired reads (default) and should already have been applied to the museum reads prior to this.

-trim 0 : We're not trimming any data

-baq 1 : estimate base alignment quality using samtools method.

###ALLELE FREQUENCY ESTIMATION

-doMajorMinor 4 : Force Major allele based on reference. The minor allele is then inferred using doMajorMinor 1. This option needs to be used when calculating SFS for multiple populations as ANGSD otherwise determines a minor allele within each population. I.e. this may not be the same across all the populations.

-ref [..fasta] : For doMM 4 above we need to specify a reference genome.



