library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)
library(data.table)

setwd("/mnt/scratch/auxie001/afum_genomes/barber.data")

dfr_full <- read.delim("barber.data_LD_out.ld.summary")
colnames(dfr_full) <- c("dist","rsq")
rsq_means_f<- aggregate(rsq~dist, data=dfr_full, FUN=function(x) c(mean=mean(x), var=var(x),count=length(x)))
print("writing full")
write.table(rsq_means_f, "rsq_means.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
rm(dfr_full)

dfr_s1 <- read.delim("barber.sample1_LD_out.ld.summary")
colnames(dfr_s1) <- c("dist","rsq")
rsq_means_s1<- aggregate(rsq~dist, data=dfr_s1, FUN=function(x) c(mean=mean(x), var=var(x),count=length(x)))
head(rsq_means_s1)
rm(dfr_s1)

dfr_s2 <- read.delim("barber.sample2_LD_out.ld.summary")
colnames(dfr_s2) <- c("dist","rsq")
rsq_means_s2<- aggregate(rsq~dist, data=dfr_s2, FUN=function(x) c(mean=mean(x), var=var(x),count=length(x)))
head(rsq_means_s2)
rm(dfr_s2)

dfr_s3 <- read.delim("barber.sample3_LD_out.ld.summary")
colnames(dfr_s3) <- c("dist","rsq")
rsq_means_s3<- aggregate(rsq~dist, data=dfr_s3, FUN=function(x) c(mean=mean(x), var=var(x),count=length(x)))
head(rsq_means_s3)
print(rsq_means_s3[1:5,3])
rm(dfr_s3)

df <- data.frame(sample=c(1,2,3),
                 mean_rsq=c(mean(rsq_means_s1$rsq.mean),2,3))

write.table(df, "rsq_samples.tab", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
