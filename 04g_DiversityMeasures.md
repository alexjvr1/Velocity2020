# Estimate diversity

## Goal

We want to test if our measures of diversity (Watterson's theta and Nucleotide diversity) are stable to differences in 

1) sequencing depth and 2) variant calls vs Genotype Likelihood based methods 


### Genotype calls vs GLs 

We're estimating GLs with ANGSD and using the previous variant calling pipeline (samtools and bcfools call) to call variants. 
These methods are similar - both use the samtools approach to estimate GLs. But the variant calling pipeline then calls genotypes using a likelihood cut-off (p 0.05). 


```


```


