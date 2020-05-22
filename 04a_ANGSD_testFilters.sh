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


