# Weekly Velocity Meetings


## 6 Nov 2020

### 1. Diversity issue

In early Oct we realised that it difficult to comparatively estimate genetic diversity in windows between the different populations. This is most likely due to different levels of missingness in the different datasets. 
While MODC and MODE have a high coverage and ~80% genotyping rate, the MUS dataset has a ~50% genotyping rate. 

Overall it seems that MUS has higher diversity than both modern datasets based on 1) the proportion of variants as estimated by the genotype calling pipeline (bcftools call), 2) Watterson's theta calculated in ANGSD. But we have two major issues: 

1) Nucleotide diversity estimates from the called genotypes and from ANGSD show the opposite pattern to Watterson's theta. We think this is due to missingness across individuals which effectively reduces the sample size to different degrees in the three datasets, but making them incomparable. 

2) If we filter the datasets for minimal missingness we lose the resolution to estimate diversity in windows across the chromosome. This window-based approach pared with the change in effective population size compared between populations will be a novel result, so I'm trying various approaches to avoid losing resolution in the MUS dataset. 

*ToDo:* 

- Estimate "ground truth" diversity from reference individual which has been sequenced to high coverage. 
AJvR will download these data from the CGR servers, map back to the reference genome, and estimate global and window-based diversity

- Apply [pixy](https://www.biorxiv.org/content/10.1101/2020.06.27.175091v1.full.pdf) which corrects for missingness within individual when calculating diversity

- Apply ROHan which finds runs of homozygosity across the genome


##### *ISSUE*

ANGSD reports a per-site value for Watterson's theta, but we're unsure how theta can be calculated on a per-site basis. 
After looking at the realSFS and thetaStat code Mark thinks there might be a default window size of 4096. 

*ToDo*

- Calculate Watterson's theta in ANGSD for windows of size 4096bp and step 1bp to check if we get the same result


##### *ISSUE*

Ilik spoke to Simon Martin about the difference in the proportion of variants we found between chromosomes using our genotype calling pipeline. Simon has a paper in prep (with Steve Montgomery from UoB) where they report a strong negative relationship between 4Dpi (i.e. 4 fold degenerate sites) and chromosome length. They attribute this to background selection in species with a large Ne. They find an average pi of 2-4% per chromosome. 

*ToDo*

- Estimate 4Dpi across LR75 chromosome


### 2. Attributing Significance to Fst outliers

Mark is working on a method to attribute a significance to the Fst outliers from our genome-wide comparisons between populations. 


*ToDo:* 

- Which loci are missing in the test saf files is using? 

See [here](https://github.com/alexjvr1/Velocity2020/blob/master/Missingness_Plots.md#missingness-in-saf-files-used-by-mark)



### 3. LD estimates

Sam is using LDhat to estimate linkage across each chromosome. 

##### *ISSUE*

Simon Martin suggests that if there is high background selection in a high Ne species (as they've found in Heliconius sp) the LDhat results will be nonsensical (see above). 


### 4. Other
