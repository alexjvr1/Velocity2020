# Plot Nucleotide diversity across the genome

See Fig 3 in [Feng et al. 2019](https://reader.elsevier.com/reader/sd/pii/S0960982218316099?token=65EE80B75A634FD527312349DBECB9E11DFF80648448692DC1E7DA07DD7E7A6552DFB5DD8B92B7420E4E045F526A4075#go_to_%22%C3%BE%C3%BF\u0000f\u0000l\u0000i\u0000n\u0000k\u00004%22)


Thetas calculated in ANGSD in windows (-win 50kb -step 10kb)

## Plots in R

Copy diversity estiates to mac
```
pwd
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/DiversityEstimates

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/M*/*pestPG .

```


R
```
library(reshape2)
library(ggplot2)



```
