# Summary of ANGSD filtering

Details of the pre-processing and ANGSD steps and scripts can be found in the [README.md](https://github.com/alexjvr1/Velocity2020/edit/master/README.md) file. 

I've summarised them here and linked specifically to the pages where I've investigated the effect of depth on ANGSD estimates, as well as the estimates of genetic diversity. 

Finally I link to the scripts used to generate the data I've sent Mark

## Note on file naming convention

A typical file will be named as follows: 

```
MUS.LR761675.1.minDP20.MinIND10.HumanReadable.10k.saf.gz
```

pop.chromosome.filters.any additional info.file type

Files start with the population (MUS/MODC/MODE), followed by the chromosome name (LR7616xx.1). We're initially working with LR761675.1 which is the shortest chromosome (6.2Mb). 




## 1. Pre-processing 

1a. Trim adapter sequence using cutadapt

1b. Concatenate raw museum data for samples that have been sequenced twice.

1c. Repair problems in museum paired end (PE) data for data from 1b. (BBrepair)

1d. Merge overlapping PE reads in museum data (BBmerge)

2. Map and process

2a. Map museum and modern data to Sanger genome using BWA mem

2b. Correct museum data for possible deamination (MapDamage -> output = corrected bam file)


## 2. ANGSD

Depth estimated in ANGSD [here](https://github.com/alexjvr1/Velocity2020/blob/1dfabf272901cc963ac2f59174e409ffada37244/03_DepthEstimate.md#angsd) 

Observed heterozygosity vs depth [here](https://github.com/alexjvr1/Velocity2020/blob/1dfabf272901cc963ac2f59174e409ffada37244/03_DepthEstimate.md#obs-het-vs-depth)

Missingness in the raw datasets [here](https://github.com/alexjvr1/Velocity2020/blob/master/Missingness_Plots.md)

PCA after removing missing data [here](https://github.com/alexjvr1/Velocity2020/blob/master/04c_PCAngsd.md)

Fst estimates using Bhatia method [here](https://github.com/alexjvr1/Velocity2020/blob/1dfabf272901cc963ac2f59174e409ffada37244/03_DepthEstimate.md#1-fst)

Using some basic ANGSD filters I compared theta between datasets [here](https://github.com/alexjvr1/Velocity2020/blob/01406f6b6dba49a140c7499dd703bd8e44b63998/NucelotideDiversityPlot.md#plot-nucleotide-diversity-across-the-genome)

Proportion of variants from called genotypes (vcf files) [here](https://github.com/alexjvr1/Velocity2020/blob/01406f6b6dba49a140c7499dd703bd8e44b63998/NucelotideDiversityPlot.md#vcf-files-allsites)

Missingness in the saf files I generated for Mark [here](https://github.com/alexjvr1/Velocity2020/blob/master/Missingness_Plots.md#missingness-in-saf-files-used-by-mark)


## 3. Script

The saf files I've sent were created using the following options: 

** Note I used a global minDP 20 for MUS and MODC, but minIndDP 2 for MODE. This should equate to the same filter (2x minInd depth x 10 indivs = 20x Global depth). 

```
time $angsd -b $POP.$N.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 \
-remove_bads 1 -only_proper_pairs 0 -r LR761675.1 \
-GL 1 -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna \
-doCounts 1 -setMinDepth 20 -setMaxDepth $maxDP -doMajorMinor 4\
 -out $SPECIESDIR/04b_ANGSD_FINAL/SFS_and_Fst/$POP/$POP.$REGION.minDP20.MinIND10 \
 -C 50 -baq 1 -dumpCounts 2 -doDepth 1 -doGlf 2 -minInd 10
```

minQ 20 : minimum base quality PHRED score of 20

minMapQ 20 : minimum mapping PHRED score of 20

uniqueOnly 1 : only uniquely mapped reads

-remove_bads 1 : remove reads with bad flags in bam file

-only_proper_pairs 0 : use all reads. MUS has merged reads, so I'm using flag 0 here

-r LR761675.1 : region or chromosome to analyse

-GL 1 : use Samtools Genotype likelihood model

-doSaf 1 : perform multi-sample GL estimation

-anc and -ref : use the reference genome as the ancestral and modern reference. The ancestral reference is used to identify the reference and alternate alleles. I've used this method to be able to compare between populations. 

-doCounts 1 : used in depth estimates

-setMinDepth 20 : minimum depth per site

-setMaxDepth : max depth set at 2xmean depth + SD

-baq 1 : estimate base alignment quality using samtools method. BAQ1 and BAQ2 give me very different final variant counts. BAQ1 is more stringent, but could be removing more loci. BAQ2 has a higher number of false positives.

-doDepth 1 : write depth per site

-doGlf 2 : write beagle output 

-minInd 10 : include only sites with data in at least 10 individuals. 




Create human readable saf files: 
```
module load languages/gcc-6.1
~/bin/angsd/misc/realSFS print MUS.LR761675.1.minDP20.MinIND10.saf.idx > MUS.LR761675.1.minDP20.MinIND10.HumanReadable.saf
```


Write the first 10k positions to file:
```
sed -n 1,10000p file.HumanReadable.saf > file.HumanReadable.10k.saf
```



Find the overlap between datasets by extracting sites in MODC.sites from MODExx10k.saf
```
awk -F "\t" 'NR==FNR {id[$1]; next} $2 in id' MODC.sites MODE.LR761675.1.minInd10.minDP2.HumanReadable.10k.saf > MODE.subset.saf

##write sites to MODE.sites

awk -F "\t" '{print $2}' > MODE.sites
wc -l MODE.sites 
    4522 MODE.sites
    
wc -l *4522.saf
    4522 MODC.subset4522.saf
    4522 MODE.subset4522.saf
    4522 MUS.subset4522.saf
```
