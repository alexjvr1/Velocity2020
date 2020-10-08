# Plot Space vs time

Plotting for three comparisons: 

Filtered Fst plots: 

1. [MODC:MODE](https://github.com/alexjvr1/Velocity2020/blob/master/04_FSTplots.md#modcmode)

2. [MODC:MUS](https://github.com/alexjvr1/Velocity2020/blob/master/04_FSTplots.md#modcmus)

3. [MODE:MUS](https://github.com/alexjvr1/Velocity2020/blob/master/04_FSTplots.md#modemus)

See [Initial Test plots](https://github.com/alexjvr1/Velocity2020/blob/master/04_FSTplots.md#initial-test-plots) below



## Bhatia's Fst: 3pops

```
/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst

##calculate folded 2D sfs
~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.saf.idx MODE/MODE.LR761675.1.saf.idx -fold 1 > MODC.MODE.LR75.ml
~/bin/angsd/misc/realSFS MODC/MODC.LR761675.1.saf.idx MUS/MUS.LR761675.1.saf.idx -fold 1 > MODC.MUS.LR75.ml
~/bin/angsd/misc/realSFS MUS/MUS.LR761675.1.saf.idx MODE/MODE.LR761675.1.saf.idx -fold 1 > MUS.MODE.LR75.ml

# Index and calculate sliding window Fst

~/bin/angsd/misc/realSFS fst index MODC/MODC.LR761675.1.saf.idx MODE/MODE.LR761675.1.saf.idx -sfs MODC.MODE.LR75.ml -whichFst 1 -fold 1 -fstout MODC.MODE.LR75.Bhatia.here

~/bin/angsd/misc/realSFS fst index MODC/MODC.LR761675.1.saf.idx MUS/MUS.LR761675.1.saf.idx -sfs MODC.MUS.LR75.ml -whichFst 1 -fold 1 -fstout MODC.MUS.LR75.Bhatia.here

~/bin/angsd/misc/realSFS fst index MUS/MUS.LR761675.1.saf.idx MODE/MODE.LR761675.1.saf.idx -sfs MUS.MODE.LR75.ml -whichFst 1 -fold 1 -fstout MUS.MODE.LR75.Bhatia.here

~/bin/angsd/misc/realSFS fst stats2 MODC.MODE.LR75.Bhatia.here.fst.idx -win 50000 -step 10000 > MODC.MODE.LR75.Bhatia.win50k.step10k.fst
~/bin/angsd/misc/realSFS fst stats2 MODC.MUS.LR75.Bhatia.here.fst.idx -win 50000 -step 10000 > MODC.MUS.LR75.Bhatia.win50k.step10k.fst
~/bin/angsd/misc/realSFS fst stats2 MUS.MODE.LR75.Bhatia.here.fst.idx -win 50000 -step 10000 > MUS.MODE.LR75.Bhatia.win50k.step10k.fst


#Or Index all three pops for window-based analysis. Make sure the populations are in the correct order
~/bin/angsd/misc/realSFS fst index MODC/MODC.LR761675.1.saf.idx MODE/MODE.LR761675.1.saf.idx MUS/MUS.LR761675.1.saf.idx -sfs MODC.MODE.LR75.ml -sfs MODC.MUS.LR75.ml -sfs MUS.MODE.LR75.ml -whichFst 1 -fold 1 -fstout 3pop.LR75.Bhatia.here

~/bin/angsd/misc/realSFS fst stats2 3pop.LR75.Bhatia.here.fst.idx -win 50000 -step 10000 > 3pop.LR75.win50k.step10k.fst

```



## Filtered Fst plots

Plots based on min global depth 20X and minInd 10X (ie min 2X each)

On mac
```
/Users/alexjvr/2020.postdoc/Velocity/E3/ANGSD_FINAL/FstPlots


##copy sliding window fst from the server 

scp bluecp3:/newhome/aj18951/E3_Aphantopus_hyperantus_2020/04b_ANGSD_FINAL/SFS_and_Fst/*Bhatia.win50k.step10k.fst .

#Add "fst" to header in each file
```


### MODC:MODE

```
library(ggplot2)

fstMODC.MODE <- read.table("MODC.MODE.LR75.Bhatia.win50k.step10k.fst", header=T)

ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761675.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("Mod Core vs Expanding LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
```

### MODC:MUS
```
library(ggplot2)

fstMODC.MUS <- read.table("MODC.MUS.LR75.Bhatia.win50k.step10k.fst", header=T)

ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761675.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("Mod Core vs Museum LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
```


### MODE:MUS
```
library(ggplot2)

fstMODE.MUS <- read.table("MUS.MODE.LR75.Bhatia.win50k.step10k.fst", header=T)

ggplot(fstMODE.MUS[which(fstMODE.MUS$chr=="LR761675.1"),], aes(x=midPos, y=fst, colour=pop)) + geom_point() + scale_fill_manual(values=c("#2E8B57", "#DAA520"))+ ggtitle("Mod Exp vs Museum LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
```

Comparison with previous plots

Top = new, filtered dataset. Bottom = previous plots (minDP 2x)

![alt_txt][MODCvsMODE.old]

[MODCvsMODE.old]:https://user-images.githubusercontent.com/12142475/94155559-f9dfa100-fe76-11ea-8371-255c2e0c0a48.png



![alt_txt][MODCvsMUS.old]

[MODCvsMUS.old]:https://user-images.githubusercontent.com/12142475/95355334-57371180-08bd-11eb-8204-23ea943dc64d.png



### E3ALL
```
library(dplyr)

fstMODC.MODE$pop <- "MODC.MODE"
fstMODC.MUS$pop <- "MODC.MUS"
fstMODE.MUS$pop <- "MODE.MUS"

fstALL <- bind_rows(fstMODC.MODE, fstMODC.MUS)
fstALL <- bind_rows(fstALL, fstMODE.MUS)

ggplot(fstALL[which(fstALL$chr=="LR761675.1"),], aes(x=midPos, y=fst, colour=pop)) + geom_line() + scale_color_manual(values=setNames(c("#2E8B57", "#DAA520", "#B0C4DE"), c("MODC.MODE", "MODC.MUS", "MODE.MUS")))+ggtitle("Ringlet population pairwise Fst: LR761675.1 - PX11") +  theme(axis.title.x = element_blank())

```

Bhatia et al. 2013 Fst estimator is more robust to differences in sample sizes between populations than the standard W&C estimator which is the default used in ANGSD. Here I compare the two: 

Top=Bhatia

Bottom=W&C

![alt_txt][Fst.all2]

[Fst.all2]:https://user-images.githubusercontent.com/12142475/95356013-1c81a900-08be-11eb-85d0-fa1a9e2e9934.png



# Compare Fst catagories: 

For all the significant peaks in each population, calculate deltaFst between datasets. 

e.g. for MODC-MODE vs MODC-MUS (ie. space vs time), we're using the top 1% of loci for the initial test: 

```
par(mfrow=c(2,1))
hist(fstMODC.MODE$fst)
hist(fstMODC.MUS$fst)


#Find the top 1% of loci in each dataset

subset(fstMODC.MODE, fst > quantile(fst, prob=1-1/100))
                                                 region        chr  midPos
138 (1029003,1060753)(1380000,1429999)(1380000,1430000) LR761675.1 1405000
139 (1035200,1070650)(1390000,1439999)(1390000,1440000) LR761675.1 1415000
387 (2943407,2984723)(3870000,3919999)(3870000,3920000) LR761675.1 3895000
388 (2950041,2987771)(3880000,3923199)(3880000,3930000) LR761675.1 3905000
389 (2959881,2993662)(3890000,3939640)(3890000,3940000) LR761675.1 3915000
390 (2969198,3000162)(3900027,3949935)(3900000,3950000) LR761675.1 3925000
391 (2976795,3005600)(3910236,3959999)(3910000,3960000) LR761675.1 3935000
    Nsites      fst       pop deltafst
138  31752 0.292400 MODC.MODE 0.204360
139  35452 0.275554 MODC.MODE 0.205668
387  41318 0.269946 MODC.MODE 0.211235
388  37732 0.287049 MODC.MODE 0.223936
389  33783 0.331218 MODC.MODE 0.274496
390  30966 0.372452 MODC.MODE 0.315517
391  28807 0.325008 MODC.MODE 0.308476


subset(fstMODC.MUS, fst > quantile(fst, prob=1-1/100))
                                             region        chr  midPos Nsites
50        (7879,8541)(500818,548345)(500000,550000) LR761675.1  525000    664
95     (13027,13915)(953530,995214)(950000,1000000) LR761675.1  975000    890
96    (13208,13988)(961825,1008492)(960000,1010000) LR761675.1  985000    782
97    (13433,14192)(971141,1019488)(970000,1020000) LR761675.1  995000    761
98    (13747,14192)(981408,1019488)(980000,1030000) LR761675.1 1005000    447
99    (13779,14259)(993001,1037140)(990000,1040000) LR761675.1 1015000    482
584 (92147,92384)(5841226,5878539)(5840000,5890000) LR761675.1 5865000    239
         fst      pop
50  0.168106 MODC.MUS
95  0.225792 MODC.MUS
96  0.243193 MODC.MUS
97  0.261846 MODC.MUS
98  0.302631 MODC.MUS
99  0.320239 MODC.MUS
584 0.159993 MODC.MUS

#Write the position of each to file 

subsetfstpos <- bind_rows(subset(fstMODC.MODE, fst > quantile(fst, prob=1-1/100),midPos), subset(fstMODC.MUS, fst > quantile(fst, prob=1-1/100),midPos))

#Extract these rows from both files

TimeSpace.fst.top0.01 <- left_join(fstMODC.MUS[fstMODC.MUS$midPos %in% subsetfstpos$midPos,], fstMODC.MODE[fstMODC.MODE$midPos %in% subsetfstpos$midPos,], by="midPos", suffix=c(".time", ".space"))

#Calculate the difference in fst Time-Space
TimeSpace.fst.top0.01$deltafst.timeminusspace <- TimeSpace.fst.top0.01$fst.time-TimeSpace.fst.top0.01$fst.space

#plot
ggplot(TimeSpace.fst.top0.01, aes(x=midPos, y=deltafst.timeminusspace))+geom_point()

ggplot(fstMODC.MODE, aes(midPos, deltafst))+ggtitle("delta fst: Time - Space")+ylab("delta fst")+geom_point(aes(colour=cut(deltafst, c(-Inf, -0.1, 0.2, Inf))))+scale_color_manual(name="deltafst", values=c("(-Inf,-0.1]"="green" ,"(-0.1,0.2]"="black","(0.2, Inf]"="orange"), labels=c("Space not Time", "Neutral/Both", "Time not Space"))


ggplot(TimeSpace.fst.top0.01, aes(fst.time, fst.space))+ggtitle("Top 1% Fst loci: Time - Space")+ylab("Fst: Space")+ xlab("Fst: Time")+geom_point(aes(colour=cut(deltafst.timeminusspace, c(-Inf, -0.1, 0.2, Inf))))+scale_color_manual(name="deltafst", values=c("(-Inf,-0.1]"="green" ,"(-0.1,0.2]"="black","(0.2, Inf]"="orange"), labels=c("Space not Time", "Both", "Time not Space"))
```

![alt_txt][hist.fst]

[hist.fst]:https://user-images.githubusercontent.com/12142475/95358065-8ac76b00-08c0-11eb-8f19-0467f9b2d13e.png


DeltaFst

![alt_txt][deltaFst]

[deltaFst]:https://user-images.githubusercontent.com/12142475/95445187-a54c2380-0956-11eb-9ca3-c31faf952a13.png


Outlier Fsts Time vs Space

![alt_txt][FstvsFst]

[FstvsFst]:https://user-images.githubusercontent.com/12142475/95445078-7d5cc000-0956-11eb-823f-f76bb58780be.png




## Blast outliers

Sam is doing this part: 

Use "bedtools closest" to find the genes in the region closest to the two fst peaks using the A.hyperantus annotation on [NCBI]()

```
LR761675.1      3900000 3950000 LR761675.1      XP_034838826.1
LR761675.1      3900000 3950000 LR761675.1      XP_034838827.1
LR761675.1      3900000 3950000 LR761675.1      XP_034838828.1
LR761675.1      3900000 3950000 LR761675.1      XP_034839076.1
LR761675.1      3900000 3950000 LR761675.1      XP_034839078.1
LR761675.1      990000  1040000 LR761675.1      XP_034838757.1
LR761675.1      990000  1040000 LR761675.1      XP_034838758.1
LR761675.1      990000  1040000 LR761675.1      XP_034839043.1
LR761675.1      990000  1040000 LR761675.1      XP_034839047.1
```


And these are: 

ruvB = Time candidate potentially associated with [immunity](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0117200). 
```
>XP_034838757.1 uncharacterized protein LOC117994891 [Aphantopus hyperantus]

MVPRATQAAGTTQIQKVVGALHGLLPKFGGPREGVRRLYAGVVRSIALYGAPVWSQRLPGVQRLRKKLNSVQRKIAIRVARLRGGYYSGSFSASGYIGGHGCETVYPDPRSLREECRRVKLVNWRKRLDEERNARQRTVGAILPNFEAWLGRKHGSVTFRLTQVLTGHGCFGEYLCRIGREMTPRCHHCEEDRDSAQHTLEECPTWESERCTLVGCIGRDLSPPAVIAAMLEDNQAWRTVVSFCETVMFKKEAAERVRERADPARRRRPDQSNQ

>XP_034838758.1 xaa-Pro aminopeptidase ApepP-like [Aphantopus hyperantus]

MVAAGSAITKANRQTSLERLMAVREVMSRTHVHAFVVPTADAHNSQYIAPSDARREWLSGLSGSSGTALVTQDHALVWTDGRYFTQFEIEVDTTVWTLMRQGTDLSLQAWLAQNMRRDSVVGIDPTTYTRSAWTSLENALLPNNVTLLAVSENIVDEARVNIGDVPPPRPNNPLIPLGVEFTGRVSSDKILELLRQIRARDAAALVLTALDDIAYTLNLRGTDIPYNPVFFAYLVIRADVQTKNVILFWGNGQLSSDIQQHLQSEGTQVEIQSYGQIFPYLRNMANELQAGSTIWLSRDGSHAIYLAAESNNAINTLSTVSPVALMKCVKNDVELQGFRSAHIKDGIAVVRFLRWVHERVAAGVNVTEIDVVDKLDELRGEENLYRGPSFATIAGAGENGAIIHYKPSREGVQKVILRDHMLLVDSGGQYMDGTTDITRTRHMSSSPTPAQRLAFTRVLKGQILLGTAVFPRGTTGNILETLARKALWDVGLNYAHGTGHGVGHFLNVHEGPSGIGAALMSGDPGIVPGMIFSNEPGYYEVGEYGIRHEDLVEVIEINNRTEDHIFAEGMVGDFLGLGAVGFRTIALAPHQTACLDVKLLNNDEIRYLNSYHARVLATLGPILKARNLNADYDWLEKECAPIRSGAIFVASSPFLVLSVVYSLWFGA

--

>XP_034838826.1 zinc finger protein 43-like [Aphantopus hyperantus]

MYKMLGTDHSRLWQMVLFAEGSEPSSMEFKCYPHNLTFSKIEQYSTINSQPIEVQTVDHTIKSEPVEADIEIKKEHNNISDFDNDIPLVQVTDIYRYATSVAPSGPLARAKSNTKTDESLPQITKRKRRSKKTESSLTQSKKKPKSKEKKVVRKKVDRCEDKVKIIVLSEQEMLEDRKLNASKSRYLRLPYKCEKCITSFDHEVNLVDHINKKHKKEEGSQVCNICESVLNTKISFEEHYKRHYRRYECIKCGRRDNNVYSVLQHYKETHGKVDMTFTCQLCDFTTESIHKLSNRTYSCEPCGKVYRNKSGLSTHMALKHSDMKPAYCVPCNQSFPSESSLRHHLRKASKHISEEDKKFTCSDCNAKFLAKRQLQEHIDWVHLKCTNHTCNKCSKSFKCERNLKRHLMYVHDKIRPPRNKICDHCGRAFQTSAILRCHIRTHTGERPLQCTHCAATFAHSAALYTHKKLVHNKQKK

>XP_034838827.1 uncharacterized protein LOC117994961 [Aphantopus hyperantus]

MDQSLLAILHSMKTKLQNLSTNASQIIETLDDPTTQSPTIIQEAEIISSKLQTLLNRLSSELTNYFRLSNDPAPDDISEISTVQLEAEDILCELKVKIKEKIKVIERENQVNRDTNPNLNSRLPKLSLPEFDGDILQWSTFWDQFSSNIDKRNLTDVDKLLYLKASLKGEAKKIVEGLETTNKNYSIALVTLKNRYGKENHIIDAHYSALYRVKASTAMNVTEVRQTLNEIERHLRVLKSLGEDINHNHLRFLIMEKFPQEIIYEMKMKINTESIEEIRKCLEVIITAREDAERVIQGKSEDNYTTQSLHANASGDNNMKLKSKRPTSKEMSTKKQDKNRGLFQGRQRSFRKQFKRKWEEPVGNQERQHEKKRKLHCIFCQENHFNDQCKAFTTLRERKNKLINRCYACLRIGHTVKTCVRKPKCAHCNEIGSHNRALCPRTLSKEIEPTNTLMQLNAKGTTILQTAVVNAKSEKEGHVETKCRILLDSGSQRTYITREVAKQLNLPIEEESRLSIFTFGSKTPQIIDSPIVKIKLLTRTNETLLLYANVVPYISQCVPYPDTELGHWENKTVLADDGSLSSRIDILIGNDYYHNVMKTGKLKIRENLYLVNSKLGWILSGETNLKPVDELSVITYFLSCGETKLNKPDLPLNNVDIKSLWDLESIGITDSPKISADEEAVKTFNESTEVIDNRYTVSWPWTEYPPSLSTNFGLSFGRLVSLLKRLDSNTLLEYDNILKEQISKGVIELVSKHTNQEKHPVHYLPHHCVRQKEKPIRIVYDASAKTKDSNSLNECLYCGPLMLEDLTGLLVKFRYHQIALSADVEKAFLQIALHERDRDVTRFLWIKDLNKPPTEDNIIYMRFCRVPFGIISSPFLLNATIKLHLSNDDKPKVRDLANDIYVDNLVTGSSSVAEALELYNDSREAFKQMSMNLREWSSNSYEFTEQIPDSSKEIRVKILGLDWDLKEDLLQLRYKIDSNANNKREVLRVIASIYDPCGFVAPHILPAKLFLQELWKTKIKWDTTFSKQMKEEWSTIREQLDKIKEISIPRCYMANAQNNDVQLHCFTDASLKAYAATVYIVNENKISFVIGKSRLVPIKDQDNLKIPRLELLGVLIGSRLIKFIRNFLPAKVTHQFLWTDSQIVIGWCKSSKLLPPFVARRIQEIKRNKELIIRYVPSELNPADVATKPYSTCEDKVKWLAGPQYLLQACDTWPQKTTNTNSLLSREDLPNESEDEEVEMDETISANCNKRLDEINNQCKINEKEKEQLQETYLEEIKKLQEEYFKDEVNGKDTNLSRNLGLFRDIDGILRCKGRMRNTNWSFDKRYPILLPKDSQLTNRIIMETHNKNYHVGVNHTLSIIRETYWILQGKRQVQKLLKKCPTCVKHGGGPYKLPPTPALPPERVNYSSPFTFTGIDYLGPVLVNNGNGLEKRWICLFTCLAIRAIHLEVVQDLSAEEGLRALRRMIATRGLPTLITSDNALQFKLMSEILAKPYCIENKIRWRFIPALAPWFGGFYERLMGIVKHCMKRTLQKHCLKDCELSTVMKEIEAVVNTRPLTHIDAELDHILKPADFLAVGKCITVENTRKEPLAEGTVTKRELVKGWKRALKIQEEFREMFSNRYLLSLRERYQHSHRDPRVTSKLEPQEGQIVQIKGEHPNRESWRVGKITKLMKGKDGFSRTAKVKVDNSEYIRSISHLYPLEMDEQEAHALADKPRPSYNYIPSERITREVHEETTGETSEIRSGEIDAPEKENDHTIQDHNETMLEDQLLPEEITNENEIDLSSEEPRPKREAAVRALEKIKEWTQNLIVLL

>XP_034838828.1 zinc finger protein 721-like [Aphantopus hyperantus]

METKIAITCSGCLNNGRRMVEISEINLQKCFSEIISENMISRYSVGRLPLQLCFECAALLLKFTRFKLQVAASYKTFQNLNSSLDRKTIIRLQTNRLPNMDMLINCHSIEVQDVAEEQVKEESGTDVVLDERLDGTEENLKEELSSHEDEIPLIKIKKKKKKLNVKLEKDGELEFIEVLLDEKEIKQERVMLAMRDDYVNAMFTCDRCILTFPNEDDLNDHILVKHQQNASLYKCSICTCSFSSEVSLTYHTHKHTRRYQCTVCADRFASKRSAVKHYNLHHCPGTAIEYQFEENNLEVNDPNAVKDPEENNVSTAEENSFPCDFCAKIFKWKTSLRKHLEKHRIETGQKRKPYCAPCRLSFTTTPNLQKHVRTSSKHKIQLKLRKLNGNSESNDTGATQKQSKQRTKLMDEIKSSVNNSQPKYLCHQCDKVFLWRGNLYRHLESHAAKAKGDLVCKPCNRTFSSIATYQQHMKISKKHASENDFKYMCSECGKRFPTKSHLKDHINWEHLKNYVYTCDECQKVFKTSNSVYLHKQAVHKKDCMEHLCDHCGKGFPNNTKLRNHILGIHSGGALRQCPWCDARFSWQSCLSRHVRQKHRKIASKMS

--

>XP_034839043.1 ruvB-like 2 [Aphantopus hyperantus]

MASLAAAQVQEVRSITRIERIGAHSHIRGLGLDDALEPRAVSQGMVGQKMARKAAGVVLQMIREGKIAGRAVLLAGQPGTGKTAIAMGLAQALGPDTPFTSMAGSEIYSLEMSKTEALTQAIRKSIGIRIKEESEIIEGEVVEVIVERAAGGGGARTGRLTLKTTDMETNYDMGTKMIDSLLKEKVQAGDVITIDKATGKINKLGRSFARARDYDATGQQARFVQCPEGELQKRKEVVHTVTLHEVDVINSRTHGFLALFSGDTGEIKAEIREQISGKLQEWREEGKAEMIPGVLFIDEAHMLDIECFSFLNRALESETAPVVIMATNRGITRIRGTDYRSPHGIPLDLLDRMIIVPTSPYSHQELREILNIRCEEEDCQMSGDALTVLTRVATETSLRYAIQLITTASLVAKRRKAAEVSMEDVKKVYSLFLDEHRSEQFLKEYQDEFMFSDGSGDGSQFMEVSQ

--

>XP_034839047.1 uncharacterized protein LOC117995157 [Aphantopus hyperantus]

MDGIQMEPIGKYKADRNGSAKAAAWRNMVDQDLDDVAFLADTDKAREARELGDAPVAVCSQRRAIFIAALVFIAMLTTALIIVYATPQPECPCLAESDLMPGQPPTGPEENATAANKERIANNGMVFPWRGARLPTFLIPKHYSLWLHPNLTTGELRGEVSIDFKVDRDTIFVVLNVRDMNVTERALFKTGGSLGPKVVKTLDYPPADQTYIEFKEKLRRKFNYTLSLRFITKLDRSDKQRGFFLAGTHRHRCAVSRFWLTHARSAFPCLDEPHLRATFKLTVVRDRFHVSLTNMPIVATEEAGFYLGHRLLQDEFAVSPPMAPHMMALAVCRLQRRAAPPLPTNTSDVTEVPDELVLPPPEISLYSDQQVILDESGPLLEWVQKTIQLFGYELNTSYPLPKFDIVVVEGGNQYSEGWGLITLSPSMLSDTKVIARLLAQQWFGGLVSPRWWSAQWLLEALTSVLGERAPALGGAPAEREQDAILLDHVLPALRLDSSTSVRAVASPRLERADIESAADELSLHKGAAIVSMAIEAAGQEAGRAALARLLRDHRSASADARDLWRALQHDGTDEAPAHAWDGWCERPGYPLLCAKTSGTDVVLKQERFVMSADPPDPEPPMLNHLLTLELRTELDELFFEPENFTEIMEEDINATTTTSVPTTTTKKPTVKTTKAPPTPKWVIPVTFSVGPLERDNEELEKFTKFWKNVSENVHVETWYDIVNTNETKKNSRWSENVTYLVWMNETEMVIPELGKHTWVRYNVGARGLYRVAPQDRSAGEAAAQRAAALYESGAPAERALLLDDAFVLSRAGRLPASRAIAAAAHLTTERHWAVWRVVVGHLSWWRELLQVSASGPRVLALLSTLPPRLPLYSGQDLQEAAVNEDQLWLSGALLTAGVEWGNENVTEQALELFDAWLNDNETIPEIYQEATFTAGVRKHGELAWQACWGALLASHAAARPLYSHAALLAALSSPVDDWLLYRFAYTVLSSEAQRGREWKTWVTALCAGVCRWRGAAGIWRILQMAPLPPSALQAAARCLHQPGDYYRFKELFGEQKGALLALDTISLNAAWVSRADADLLAYFNSVHKKSG

--

>XP_034839076.1 uncharacterized protein LOC117995182 [Aphantopus hyperantus]

MHESQVTQLNLLQPQYKHLQPVQEQNVPQNTNPPEEQQPVKPRVENKKKKEAIDGQAREIIFKVIKFFESEKQNRGYHFPVENVVKRACAATGLSESTIKRIKREGLRAEANQTKMTGPKKKRVRKTKVQLDYFQLCALRGIVNAYSMRKEVPTLGKILTAAKHELNYRGGKESLRLILLNKLGIKFKKCEKKNKKPPEPVAPTQQMVPVMSHIQMQPMKPDNQCIYTNMMPQVPPVSY

--

>XP_034839078.1 uncharacterized protein LOC117995184 [Aphantopus hyperantus]

MICLDADSKLFLMNKHKLDEAYEKLTGHPLCDQGNLKQTVCVQCAQTLINFSRFRDKSLRARTLMMDLVNKHELITRRHIKMINRTKHKLTSNMVVTTLGPDHCNLYILEHPSEDKQTELEETSHQVLVKTEGSDDSMSVDEDVEVINEDHNNIDIVKDEFVTSNDEDISDYSIMMEAKAMDEDLYKALKMKHPYMQHTSAKFEASRPQPFCTFCSDLSRSSSRDKQRNKSSAN
```






###############################
###############################
## Initial test plots

These plots were made using the previous filters (min global depth 2X and no minInd filter)

```
alexjvr/2020.postdoc/Velocity/E3/Test.ANGSDstats/
```



MOD core vs exp
```
p2 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761674.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("Mod Core vs Mus LR761674.1 - PX02") +  theme(axis.title.x = element_blank())
p3 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761666.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761666.1 - PX03") +  theme(axis.title.x = element_blank())
p4 <-ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761647.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761647.1 - PX04") +  theme(axis.title.x = element_blank())
p5 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761648.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761648.1 - PX05") +  theme(axis.title.x = element_blank())
p6 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761657.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761657.1 - PX06") +  theme(axis.title.x = element_blank())
p7 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761668.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761668.1 + PX07") +  theme(axis.title.x = element_blank())
p8 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761662.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761662.1 - PX08") +  theme(axis.title.x = element_blank())
p9 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761659.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761659.1 - PX09") +  theme(axis.title.x = element_blank())
p10 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761649.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761649.1 - PX10") +  theme(axis.title.x = element_blank())
p11 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761675.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
p12 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761655.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761655.1 - PX12") +  theme(axis.title.x = element_blank())
p13 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761653.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761653.1 - PX13") +  theme(axis.title.x = element_blank())
p14 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761669.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761669.1 - PX14") +  theme(axis.title.x = element_blank())
p15 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761651.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761651.1 - PX15") +  theme(axis.title.x = element_blank())
p16 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761656.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761656.1 - PX30/16") +  theme(axis.title.x = element_blank())
p17 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761660.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761660.1 - PX17") +  theme(axis.title.x = element_blank())
p18 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761663.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761663.1 - PX18") +  theme(axis.title.x = element_blank())
p19 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761667.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761667.1 - PX19") +  theme(axis.title.x = element_blank())
p20 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761654.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761654.1 - PX20/28") +  theme(axis.title.x = element_blank())
p21 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761664.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761664.1 - PX21") +  theme(axis.title.x = element_blank())
p22 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761652.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761652.1 - PX22") +  theme(axis.title.x = element_blank())
p23 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761658.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761658.1 - PX23") +  theme(axis.title.x = element_blank())
p24 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761673.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761673.1 - PX24") +  theme(axis.title.x = element_blank())
p25 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761665.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761665.1 - PX25") +  theme(axis.title.x = element_blank())
p26 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761671.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761671.1 - PX05/26") +  theme(axis.title.x = element_blank())
p27 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761670.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761670.1 - PX27") +  theme(axis.title.x = element_blank())

p29 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761661.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761661.1 - PX29") +  theme(axis.title.x = element_blank())
p31 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761672.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761672.1 - PX31") +  theme(axis.title.x = element_blank())

p32 <- ggplot(fstMODC.MODE[which(fstMODC.MODE$chr=="LR761650.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761650.1 - PX1Z") +  theme(axis.title.x = element_blank())


```




Mod vs Mus
```
p2mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761674.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("Mod Core vs Mus LR761674.1 - PX02") +  theme(axis.title.x = element_blank())
p3mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761666.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761666.1 - PX03") +  theme(axis.title.x = element_blank())
p4mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761647.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761647.1 - PX04") +  theme(axis.title.x = element_blank())
p5mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761648.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761648.1 - PX05") +  theme(axis.title.x = element_blank())
p6mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761657.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761657.1 - PX06") +  theme(axis.title.x = element_blank())
p7mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761668.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761668.1 + PX07") +  theme(axis.title.x = element_blank())
p8mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761662.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761662.1 - PX08") +  theme(axis.title.x = element_blank())
p9mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761659.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761659.1 - PX09") +  theme(axis.title.x = element_blank())
p10mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761649.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761649.1 - PX10") +  theme(axis.title.x = element_blank())
p11mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761675.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761675.1 - PX11") +  theme(axis.title.x = element_blank())
p12mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761655.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761655.1 - PX12") +  theme(axis.title.x = element_blank())
p13mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761653.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761653.1 - PX13") +  theme(axis.title.x = element_blank())
p14mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761669.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761669.1 - PX14") +  theme(axis.title.x = element_blank())
p15mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761651.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761651.1 - PX15") +  theme(axis.title.x = element_blank())
p16mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761656.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761656.1 - PX30/16") +  theme(axis.title.x = element_blank())
p17mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761660.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761660.1 - PX17") +  theme(axis.title.x = element_blank())
p18mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761663.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761663.1 - PX18") +  theme(axis.title.x = element_blank())
p19mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761667.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761667.1 - PX19") +  theme(axis.title.x = element_blank())
p20mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761654.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761654.1 - PX20/28") +  theme(axis.title.x = element_blank())
p21mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761664.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761664.1 - PX21") +  theme(axis.title.x = element_blank())
p22mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761652.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761652.1 - PX22") +  theme(axis.title.x = element_blank())
p23mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761658.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761658.1 - PX23") +  theme(axis.title.x = element_blank())
p24mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761673.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761673.1 - PX24") +  theme(axis.title.x = element_blank())
p25mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761665.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761665.1 - PX25") +  theme(axis.title.x = element_blank())
p26mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761671.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761671.1 - PX05/26") +  theme(axis.title.x = element_blank())
p27mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761670.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761670.1 - PX27") +  theme(axis.title.x = element_blank())

p29mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761661.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761661.1 - PX29") +  theme(axis.title.x = element_blank())
p31mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761672.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761672.1 - PX31") +  theme(axis.title.x = element_blank())

p32mus <- ggplot(fstMODC.MUS[which(fstMODC.MUS$chr=="LR761650.1"),], aes(x=midPos, y=fst)) + geom_point() + ggtitle("LR761650.1 - PX1Z") +  theme(axis.title.x = element_blank())
```


Write all plots to pdf and combine using an online pdf combiner. 
```
pdf("FST1.pdf")

multiplot(p2,p2mus, p3, p3mus, cols=1)

dev.off()



pdf("FST2.pdf")

multiplot(p4,p4mus, p5, p5mus, cols=1)

dev.off()



pdf("FST3.pdf")

multiplot(p6,p6mus, p7, p7mus, cols=1)

dev.off()



pdf("FST4.pdf")

multiplot(p8,p8mus, p9, p9mus, cols=1)

dev.off()



pdf("FST5.pdf")

multiplot(p10,p10mus, p11, p11mus, cols=1)

dev.off()



pdf("FST6.pdf")

multiplot(p12,p12mus, p13, p13mus, cols=1)

dev.off()



pdf("FST7.pdf")

multiplot(p14,p14mus, p15, p15mus, cols=1)

dev.off()


pdf("FST8.pdf")

multiplot(p16,p16mus, p17, p17mus, cols=1)

dev.off()



pdf("FST9.pdf")

multiplot(p18,p18mus, p19, p19mus, cols=1)

dev.off()



pdf("FST10.pdf")

multiplot(p20,p20mus, p21, p21mus, cols=1)

dev.off()



pdf("FST11.pdf")

multiplot(p22,p22mus, p23, p23mus, cols=1)

dev.off()



pdf("FST12.pdf")

multiplot(p24,p24mus, p25, p25mus, cols=1)

dev.off()



pdf("FST13.pdf")

multiplot(p26,p26mus, p27, p27mus, cols=1)

dev.off()



pdf("FST14.pdf")

multiplot(p29, p29mus, p31, p31mus, cols=1)

dev.off()

pdf("FST15.pdf")
multiplot(p32,p32mus, cols=1)
dev.off()

```
