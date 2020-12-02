# Filter ANGSD dataset using strict filters described by [Bi et. al. 2019](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008119)

## 1. Filter at the Indiv level

a. Remove indivs with very low and very high depth. 

I'll keep individuals with high depth, but remove indivs with very low depth. 

Calculate individual depths from bam files using this [script](https://github.com/alexjvr1/Velocity2020/blob/15e687209208aec5af0a1fe6a01616e31afc224a/LR75.depth.sh)

This gives us the 
```

```

b. Remove all individuals with excessively high sequencing error rate

This is tested by comparing the proportion of mismatched bases in all bases aligned to the mtDNA genome. Individuals with 3x the mean number of mismatches should be removed. 

I've asked Stefanos to look at this. 


## 2. Filtering at contig level

