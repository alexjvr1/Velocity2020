# Runs of Homozygosity and IBD

Runs of homozygosity and increased Identity by Descent (IBD) indicates a loss of diversity. We expect higher genetic divrsity in the museum data compared to modern core. 

Example of using GLs in [this](https://www.biorxiv.org/content/10.1101/2020.04.08.031344v1.full.pdf) paper by Hooper et al. (bioRxiv; Orcas)  

This [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6614887/) developed a method for identifying ROH and calculate Watterson's theta for the regions between ROH. 

It's been specifically developed with low coverage ancient DNA sequencing in mind. 


## [ROHan](http://grenaud.github.io/ROHan/)

Download on server
```
module load tools/cmake-3.8.1
module load tools/git-2.22.0
module load tools/autoconf-2.69

git clone --depth 1 https://github.com/grenaud/rohan.git
cd rohan
make
```

The executable is in the rohan folder: 
```
module load languages/gcc-6.1

/newhome/aj18951/software/rohan/src/rohan
```

We need the Ti/Tv ratio for the analysis. We can obtain that from the called genotypes called using the samtools mplileup/bcftools call pipeline. 

#Change the submission script to output ALL sites rather than just variants. Change "-v 1" to "-v 0". 

I'm running this initially only for LR761675.1. Find the intersect in the raw bcf files and calculate the Ti/Tv ratio for each population
```
module load apps/vcftools-0.1.12b

for i in $(ls *vcf); do vcftools --vcf $i --TsTv-summary

for i in $(ls *summary); do ls $i && awk -F "\t" 'FNR==8{Ts=$2} FNR==9{Tv=$2} END {print Ts/Tv}' $i; done
MODC.TsTv.summary
1.13842
MODE.TsTv.summary
1.13964
MUS.TsTv.summary
1.13316


```
