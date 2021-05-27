# Variants called

Called variants using bcftools call on filtered bam files. SNPs were called together for MODC and MODE, and separately for MUS

## E3

MUS
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03.2_Variant_calling_museum

cat job0029.20210526-101529.report.txt

Calling variants for regions: LR761675.1 using samtools/bcftools...


Finished.

Merging results for all scaffolds/chromosomes...
63560 variants called for 48 samples/individuals.

```

MODC & MODE
```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/03.2_Variant_calling/PerChr
cat job0027.20210525-152341.report.txt

Calling variants for regions: LR761675.1 using samtools/bcftools...


Finished.

Merging results for all scaffolds/chromosomes...
165029 variants called for 78 samples/individuals.
```


# Modern samples: Why are we finding no singletons in the MODE dataset??

We're first having a look at the population structure and diversity estimated from the Modern dataset



Basic filters
```
#1. Remove loci with QUAL < 100 (i.e. Phred confidence in variant site)
#2. remove all multi-allelic SNPs

vcftools --bcf job0027.variants.raw.bcf --minQ 20 --max-alleles 2 --recode --recode-INFO-all --out LR75.step1
VCFtools - 0.1.17
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--bcf job0027.variants.raw.bcf
	--recode-INFO-all
	--max-alleles 2
	--minQ 20
	--out LR75.step1
	--recode

Using zlib version: 1.2.11
Warning: Expected at least 2 parts in INFO entry: ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias",IDX=21>
Warning: Expected at least 2 parts in INFO entry: ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias",IDX=21>
Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases",IDX=24>
Warning: Expected at least 2 parts in INFO entry: ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases",IDX=24>
After filtering, kept 78 out of 78 Individuals
Outputting VCF file...
BCF file contains IDX values in header. These are being removed for conversion to VCF.
After filtering, kept 157949 out of a possible 165029 Sites
```


Estimate the depth of the remaining loci and calculate the max depth filter
```
#calculates the mean depth per indiv
vcftools --vcf LR75.step1.recode.vcf --depth

#There are 6 indivs in MODE that sequenced badly. We're removing them from the dataset
cat toremove 
AH-02-2019-74_mod.exp_R1.fastq.gz.flt
AH-02-2019-75_mod.exp_R1.fastq.gz.flt
AH-02-2019-77_mod.exp_R1.fastq.gz.flt
AH-02-2019-78_mod.exp_R1.fastq.gz.flt
AH-02-2019-79_mod.exp_R1.fastq.gz.flt
AH-02-2019-80_mod.exp_R1.fastq.gz.flt

vcftools --vcf LR75.step1.recode.vcf --remove toremove --recode --recode-INFO-all --out LR75.step2

##Depth per locus



depth <- read.table("out.ldepth.mean", header=T)
summary(depth$MEAN_DEPTH)
    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  0.01389   2.94444   3.73611   3.65662   4.19444 224.76400 
  
library(plyr)
ddply(depth,.(CHROM), colwise(mean))
      CHROM     POS MEAN_DEPTH VAR_DEPTH
1 LR761675.1 3062144   3.656624   16.2482

ddply(depth,.(CHROM), colwise(sd))
       CHROM     POS MEAN_DEPTH VAR_DEPTH
1 LR761675.1 1854566   2.925705  232.2034
```





```
#3. Max mean depth of mean + 2xSD of meanDP (here 10x)

#4. Minimum depth of 2 (i.e. remove loci with lower than mean 2x depth)


vcftools --vcf LR75.step2.recode.vcf --max-meanDP 10 --minDP 2 --max-missing 0.5 --recode --recode-INFO-all --out LR75.step3

After filtering, kept 72 out of 72 Individuals
Outputting VCF file...
After filtering, kept 137088 out of a possible 157949 Sites

```


```
#5. Remove all loci genotyped in <50% of individuals in each population
#6. Remove individuals with >60% missingness

#Split the vcffile into MODC and MODE

bcftools query -l LR75.step3.recode.vcf 

##MODE

vcftools --vcf LR75.MODE.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out LR75.MODE.0.5miss
After filtering, kept 34 out of 34 Individuals
Outputting VCF file...
After filtering, kept 136735 out of a possible 137088 Sites
Run Time = 16.00 seconds

vcftools --vcf LR75.MODE.0.5miss.recode.vcf --missing-indv

vcftools --vcf LR75.MODE.0.5miss.recode.vcf --remove MODE.remove --recode --recode-INFO-all --out LR75.MODE.FINAL
Excluding individuals in 'exclude' list
After filtering, kept 32 out of 34 Individuals
After filtering, kept 136735 out of a possible 136735 Sites

##MODC

vcftools --vcf LR75.MODC.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out LR75.MODC.0.5miss
After filtering, kept 38 out of 38 Individuals
Outputting VCF file...
After filtering, kept 113387 out of a possible 137088 Sites

vcftools --vcf LR75.MODC.0.5miss.recode.vcf --remove MODC.remove --recode --recode-INFO-all --out LR75.MODC.FINAL

vcftools --vcf LR75.MODC.0.5miss.recode.vcf --missing-indv

#Remove 8 indivs for >60% missingness
```


Intersect and combine data
```
module load apps/samtools-1.9
module load apps/bcftools-1.8
module load apps/tabix-0.2.6

bcftools isec -n 2 LR75.MODC.FINAL.recode.vcf.gz LR75.MODE.FINAL.recode.vcf.gz -p E3.LR75.FINAL
cd E3.LR75.FINAL/

##Had to add header to vcf files to say they're vcf (for bcftools error)

bgzip 0000.vcf
tabix 0000.vcf.gz
bgzip 0001.vcf
tabix 0001.vcf.gz

vcftools --gzvcf 0000.vcf.gz --freq --out MODC.LR75.freq
vcftools --gzvcf 0001.vcf.gz --freq --out MODE.LR75.freq

awk -F ":" '{print $3}' MODE.LR75.freq.frq > MODE.frq
awk -F ":" '{print $3}' MODC.LR75.freq.frq > MODC.frq


gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "SFS: E3.MODE"
set ylabel "Number of Occurrences"
set xlabel "SFS"
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'MODE.frq' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title SFS: E3.MODC"
set ylabel "Number of Occurrences"
set xlabel "SFS"
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'MODC.frq' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

```

![alt_txt][MODC.sfs]

[MODC.sfs]:https://user-images.githubusercontent.com/12142475/119901911-6d720a80-bf3e-11eb-9c59-c4475a32df06.png



![alt_txt][MODE.sfs]

[MODE.sfs]:https://user-images.githubusercontent.com/12142475/119901916-706cfb00-bf3e-11eb-9628-8dcde07c8f7d.png

