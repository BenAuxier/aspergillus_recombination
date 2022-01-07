### this are script is used to filter variants from a freebayes variant call file 
### In the files the first two columns represent parents, with column 1 being P0 and 
### column 2 being P1. P0 is similar to the reference genome, P1 is the alternative
### parent. All other columns are potential offspring. 

## filtering steps are the following
## (1) make the list of those that are differentiated between the parents
## (2) select those samples that have a mean allele frequecy >0.2 and <0.8
## (3) select those samples that have loci that are similar to one of the parents
## (4) select those loci that are biallelic
## (5) samples remain if 95% or higher  is within 0 - 0.05 or a parent 
## (6) remove loci that often have a frequency that is different from homozygosity

library(vcfR)

vcf<-read.vcfR("final.annotated.vcf")
 
## extract alternative depth
AO<-extract.gt(vcf,element='AO', as.numeric=TRUE)
## extract reference depth
RO<-extract.gt(vcf,element='RO', as.numeric=TRUE)

## calculate coverage
DP<-AO+RO

## create a alternate read frequency matrix 
AF<- AO / DP

## select those variants that are similar to reference in P0 and alternative in P1, with a minimum depth of >15, maximum of <150
W01<- which(AF[,1]<0.05 & AF[,2]>0.95 & is.na(DP[,1])==FALSE & is.na(DP[,2])==FALSE &
DP[,1]>15 & DP[,2]>15 & DP[,1]<150 & DP[,2]<150)

## select those variants that are similar to reference in P1 and alternative in P0, with a minimum depth of >15, maximum of <150
W10<- which(AF[,1]>0.95 & AF[,2]<0.05 & is.na(DP[,1])==FALSE & is.na(DP[,2])==FALSE &
DP[,1]>15 & DP[,2]>15 & DP[,1]<150 & DP[,2]<150)

## check how many are reference in P0 (which is the reference genome) and alternate in P1. 
length(W01)

## check how many are alternate in P0 (which is the reference genome) and reference in P1. (none in our case) 
length(W10)

## we make a list of chromosomes and list them as a table
table(vcf@fix[W01,1])

## get chromosome names
CHROS<- unique(vcf@fix[,1])


## (1) make the list of those that are differentiated between the parents
L<-W01

## check length
length(L)


## number of samples
D<-dim(AF)[[2]]
## mean allele frequency
meanFREQ<-colMeans(round(AF[L,1:D]),na.rm=TRUE)

## (2) select those samples that have a mean allele frequecy >0.2 and <0.8

GS<-which(meanFREQ>0.2 & meanFREQ<0.8)



## there might be samples that are not offspring of the parents. 
## these samples have a SNP frequency that is different from both partens,
## i.e., the mean allele frequency difference between the mean parents
## and offspring is 1 for many snps

MP<- rowMeans(AF[,1:2])
WO<-c()
MO<-c()
for (i in 1:length(GS))
{
MO<- c(MO,sum(abs(AF[,GS[i]]-MP)>0.99,na.rm=TRUE))
WO<-c(WO,which(abs(AF[,GS[i]]-MP)>0.99))

}

## (3) select those have loci that are similar to one of the parents

GS<-GS[which(MO<1000)]


## (4) select those loci that are biallelic

N<-c()
for (i in 1:length(W01)){
N<-c(N,length(strsplit(vcf@fix[W01[i],5],",")[[1]]))
}

Ls<- L[which(N==1)]


## now we check the proportion of SNP between 0 - 0.05 and 0.95 - 1 or alternative read frequency for each sample

propper05<-c()
for (i in 1:length(GS))
{
propper05<-c(propper05, (sum(AF[Ls,GS[i]]<0.05 & DP[Ls,GS[i]]>15,na.rm=TRUE)+sum(AF[Ls,GS[i]]>0.95 & DP[Ls,GS[i]]>15,na.rm=TRUE))/sum(AF[Ls,GS[i]]<2,na.rm=TRUE))
}

### (5) samples remain if 95% or higher  is within 0 - 0.05 or a parent 

GS<- GS[which(propper05>0.95)] ## leaves us with 195 samples

### (6) remove loci that often have a frequency that is different from homozygosity

Proximity2_0.5<-rowMeans(abs(AF[Ls,GS]-0.5),na.rm=TRUE)

## removes roughly 2.8%

Lss<- Ls[which(Proximity2_0.5>0.49)]


### write out genotype table using read frequencies

write.table(cbind(vcf@fix[Lss,c(1,2,4,5)],AF[Lss,c(1,2,GS)]), 
file="genotypes as frequencies.txt",col.names=TRUE, row.names=FALSE, 
quote=FALSE, sep='\t')


### make loci NA is read frequency is between 10 - 90%

for (i in 1:length(GS))
{
AF[Lss[which(AF[Lss,GS[i]]>0.1 & AF[Lss,GS[i]]<0.9)],GS[i]]<- NA
}

### write out alleles are rounded frequencies (so all become 0 or 1 [or NA]). 

write.table(cbind(vcf@fix[Lss,c(1,2,4,5)],round(AF[Lss,c(1,2,GS)])), 
file="genotypes as frequencies rounded.txt",col.names=TRUE, row.names=FALSE, 
quote=FALSE, sep='\t')



