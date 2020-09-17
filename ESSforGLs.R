################
##Mark Beaumont
########################################################
##Score for GLs as proxy for ESS
##To be used to test for bias in correlations with Fst
########################################################
#!/usr/bin/env Rscript
#arg 1= input file *beagle.gz
#arg 2= output file

args = commandArgs(trailingOnly=TRUE)

geno.lik.raw = read.table(gzfile(arg[1]),header=T)
geno.lik = data.matrix(geno.lik.raw[,-c(1:3)])
scoremat = matrix(nrow=nrow(geno.lik),ncol=(ncol(geno.lik)/3))
for(j in 1:nrow(geno.lik)){
	for(k in 1:(ncol(geno.lik)/3)){
		i.1 = (k-1)*3 + 1
		i.2 = i.1+1
		i.3 = i.2+1
		l1 = geno.lik[j,i.1]
		l2 = geno.lik[j,i.2]
		l3 = geno.lik[j,i.3]
		vec1 = sort(c(l1,l2,l3),decr=T)
		#This is purely a made-up, intuitive score,designed to go from 0 when 1/3, 1/3, 1/3
		#to 1 when just one of l1,l2,l3 is equal to 1.000
		#There is no theory to equate to a N_eff, but I guess there must be a monotonic relationship
            #between the two
		#The idea is that when the largest genotype likelihood is 1, the others must be 0, so the 
		#score is 1. If it's 1/3, the others must be a 1/3 also, and the score is 0
		scoremat[j,k] = 1-0.5*(vec1[2] + vec1[3])/vec1[1]
	}
	scorevec = apply(scoremat,1,sum)

}

plot(scorevec,ylim=c(0,ncol(scoremat)))
abline(h=ncol(scoremat),col=2,lwd=2)
abline(h=0,col=2,lwd=2)
write.table(scorevec, (args[2]),quote=F, sep="\t")
print(paste("number of individuals is ",ncol(scoremat)))
print(paste("maximum 'N_eff' is ",max(scorevec)))
print(paste("minimum 'N_eff' is ",min(scorevec)))
