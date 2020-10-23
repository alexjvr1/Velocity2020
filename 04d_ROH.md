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

