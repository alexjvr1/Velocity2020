# 4Dpi

Estimate pi from 4-fold degenerate sites across all individuals. 

These sites are used to estimate neutral diversity across the genome. 

It's been used on Lepidoptera by [Mackintosh et al. 2019](https://www.nature.com/articles/s41467-019-11308-4) from Konrad Lohse's group. They find a large range in 4Dpi: 

0.4 - 4% (mean 1.75%)

Simon Martin finds differences between chromosomes in Heliconius sp: 2-4% 4Dpi (pers comm to Ilik). 

## Pipeline

1. Extract genes using gff file

2. Translate nucleotides to Amino Acids. 

3. Identify the 4-fold degenerate sites

4. Calculate pi

We're using a script from the Lohse lab found [here](https://github.com/DRL/gIMble/blob/master/README.md)




