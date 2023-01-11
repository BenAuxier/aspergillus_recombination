library(tidyverse)
library(dplyr)
library(ggrepel)
library(cowplot)
library(qtl2)

if(Sys.getenv("USER")=="ben"){setwd("/Users/ben/Library/CloudStorage/OneDrive-WageningenUniversity&Research/aspergillus_recombination/Figure1and2/Fig1.lengths")}

set.seed(167436)
#first we rarify the markers as a different comparison----
#progeny2 is the map after the first estimation of RF based on haldane calculation
progeny2 <- read_cross2("asp_June9_RFadjust.yaml",quiet=F)

reduced_maps_means <- read.csv("reduced_maps_means.csv",header=T)
reduced_maps_means <- reduced_maps_means[reduced_maps_means$dist != 1,]
reduced_maps_means$plotting_heights <- reduced_maps_means$mean_len+300
reduced_maps_means$plotting_heights[1] <- reduced_maps_means$plotting_heights[[1]] +800
reduced_maps_means$plotting_heights[1] <- reduced_maps_means$plotting_heights[[2]] +400

rarified <- ggplot(reduced_maps_means) + theme_classic() + 
  scale_y_continuous(limits=c(11000,15400),breaks=c(10000,12000,14000),labels=c("10,000","12,000","14,000"))+
  geom_line(aes(x=dist,y=mean_len)) +
  geom_point(aes(x=dist,y=mean_len),size=2.5,color="orange3") +
  geom_errorbar(aes(x=dist,ymin=mean_len-sd_len,ymax=mean_len+sd_len))+
  geom_text(aes(x=dist+2,y=plotting_heights,label=mean_num),size=2.8)+
  #geom_hline(aes(yintercept=11797),lty=3)+
  labs(x="Marker distance cutoff (cM)",y="Map Length (cM)")+
  theme(panel.grid.minor=element_blank())

#next we compare observed versus hypothetical----
setwd("criteria_based_cutoff/")

#First we collect the empricial distribution
cr_progeny <- read_cross2("asp_June9_criteria_3000_2.yaml")
cr_pr <- calc_genoprob(cr_progeny,error_prob = 0.00)
cr_pr_geno <- maxmarg(cr_pr)
cr_pr_xo_loc <- locate_xo(cr_pr_geno,cr_progeny$pmap)
cr_pr_dco_dist <- unlist(c(lapply(c(0,cr_pr_xo_loc$chr1,max(cr_progeny$pmap$chr1)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr2,max(cr_progeny$pmap$chr2)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr3,max(cr_progeny$pmap$chr3)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr4,max(cr_progeny$pmap$chr4)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr5,max(cr_progeny$pmap$chr5)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr6,max(cr_progeny$pmap$chr6)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr7,max(cr_progeny$pmap$chr7)),diff),
                           lapply(c(0,cr_pr_xo_loc$chr8,max(cr_progeny$pmap$chr8)),diff)))
#get lengths of chromosomes
chr_lens <- unlist(lapply(cr_progeny$pmap,max))
#now we simulate a 1000cM map, so 20 crossovers. We first distributed according to lengths
#but we need to only distribute 12, since minimum one crossover necessary for each chromsome
#we want to simulate 1000,3000,5000,7000,9000,11000
#which is 20,60,100,140,180,220 per meiosis
#so 10,30,50,70,90,110 per meiosis
all_DCO_dists <- list()
for (map in seq(from=12,to=162,by=10)){
  total_DCO_dist <- c()
  for (rep in c(1:195)){
    spread_map <- rnorm(1,map,map/10)
    per_chr <- sample(c(1:8),spread_map,replace=T,prob=chr_lens)
    per_chr <- c(per_chr,c(1:8))
    per_chr <- table(per_chr)
    for (chr in c(1:8)){
      #print(paste0("Chromosome: ",chr," has ",per_chr[[chr]]))
      co <- floor(runif(per_chr[[chr]],min=cr_progeny$pmap[[chr]][[1]],max=chr_lens[[chr]]))
      co <- sort(co)
      corr_co <- c()
      if (length(co) > 0){
        for (ind_co in co){
          right_mar <- which(cr_progeny$pmap[[chr]] >= ind_co)[[1]]
          corr_ind_co <- mean(c(cr_progeny$pmap[[chr]][[right_mar]],cr_progeny$pmap[[chr]][[right_mar-1]]))
          corr_co <- c(corr_co,corr_ind_co)
        }
      }
      #print(corr_co)
      corr_co <- c(0,corr_co,chr_lens[[chr]])
      #print(corr_co)
      total_DCO_dist <- c(total_DCO_dist,diff(corr_co))
    }
  }
  all_DCO_dists[[map]] <- total_DCO_dist
}

D_values <- data.frame(D = c(ks.test(all_DCO_dists[[12]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[22]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[32]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[42]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[52]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[62]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[72]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[82]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[92]], cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[102]],cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[112]],cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[122]],cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[132]],cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[142]],cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[152]],cr_pr_dco_dist)$statistic,
                             ks.test(all_DCO_dists[[162]],cr_pr_dco_dist)$statistic),
                       length = seq(from=2000,to=17000,by=1000))

D_values_plot <- ggplot(D_values) + theme_classic() + labs(x="Map Length (cM)",y="D Statistic")+ylim(0,0.6)+xlim(0,17000)+
  geom_line(aes(x=length,y=D)) +
  geom_point(aes(x=length,y=D),size=2.5,col="cornflowerblue")+theme(panel.grid.minor=element_blank())

sim_12kmap <- density(from=3,to=7,n=100,log10(all_DCO_dists[[112]]))
actual <- density(from=3,to=7,n=100,log10(cr_pr_dco_dist))
min <- pmin(sim_12kmap$y,actual$y)
max <- pmax(sim_12kmap$y,actual$y)

D_example <- ggplot() + theme_classic() + labs(y="Density",x="Interevent Distance (bp)") + scale_x_log10(limits=c(1e3,1e7)) + 
  geom_density(aes(x=cr_pr_dco_dist),lwd=1) + 
  geom_density(aes(x=all_DCO_dists[[112]]),lty=3,lwd=1)+
  geom_ribbon(aes(x=10**actual$x,ymin=min,ymax=max),fill="grey70",alpha=0.3)+
  theme(axis.title.x=element_text(),
        axis.text.x=element_text())
D_example
p1 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[12]]))  + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p2 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[22]]))  + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p3 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[42]]))  + xlim(0,2e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p4 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[62]]))  + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p5 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[82]]))  + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p6 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[102]])) + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p7 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[122]])) + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p8 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity") + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + geom_density(aes(x=all_DCO_dists[[142]])) + xlim(0,1e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)
p9 <- ggplot() + theme_bw() + labs(y="\n\n\n\n\n\n\n\ndensity",x="Interevent Distance (bp)") + geom_density(aes(x=all_DCO_dists[[162]])) + xlim(0,2e6) + ylim(0,6e-6) + geom_density(aes(x=cr_pr_dco_dist),lty=3)

p1_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[12]]))
p2_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[22]])) 
p3_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[42]])) 
p4_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[62]]))
p5_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[82]]))
p6_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[102]]))
p7_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[122]]))
p8_l <- ggplot() + theme_bw() + theme(axis.text.x=element_blank(),axis.title.x=element_blank()) + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[142]]))
p9_l <- ggplot() + theme_bw() + labs(x="Interevent Distance (bp)") + scale_x_log10(limits=c(1000,6e6)) + geom_density(aes(x=cr_pr_dco_dist),lty=3) + geom_density(aes(x=all_DCO_dists[[162]]))

setwd("../")
#make supplemental plot
pdf(paste0("Supplemental.distributions-",Sys.Date(),".pdf"),width=8,height=9)
plot_grid(p1,p1_l,p2,p2_l,p3,p3_l,p4,p4_l,p5,p5_l,p6,p6_l,p7,p7_l,p8,p8_l,p9,p9_l,
          label_fontface="plain",ncol=2,hjust=-0.1,rel_widths=c(0.55,0.45),rel_heights=c(rep(1,8),1.3),
          labels=c("A)  2000 cM","","B)  3000 cM","","C)  5000 cM","",
                   "D)  7000 cM","","E)  9000 cM","","F) 11000 cM","",
                   "G) 13000 cM","","H) 15000 cM","","I) 17000 cM",""))
dev.off()

#now we use criteria based cutoff----
progeny <- read_cross2("asp_June9_RFadjust.yaml")
pr <- calc_genoprob(progeny,error_prob=0)
g <- maxmarg(pr)
xo <- locate_xo(g,progeny$pmap)

chrs <- names(xo)
strains <- names(xo$chr1)
#we need to ignore parents, since they have no crossovers
strains <- strains[strains != "NIR"]
strains <- strains[strains != "CNX"]

#cols are num markers, rows are lengths
mar1 <- c(14720,14612,14410,13825,13309,12510,12011,11406,11037)
mar2 <- c(11980,11964,11894,11642,11410,10986,10678,10349,10155)
mar3 <- c(10839,10834,10802,10670,10515,10233,10028, 9773, 9645)
mar4 <- c(10204,10202,10180,10095, 9986, 9766, 9615, 9403, 9301)
mar5 <- c( 9654, 9653, 9641, 9583, 9517, 9382, 9265, 9090, 9007)
mar6 <- c( 9259, 9259, 9255, 9205, 9157, 9055, 8968, 8841, 8769)
cutoff_df <- data.frame("1 marker"  = mar1,"2 markers" = mar2, "3 markers" = mar3,
                        "4 markers" = mar4,"5 markers" = mar5, "6 markers" = mar6,
                        lengths = c(0,500,1000,2000,3000,5000,7000,10000,12000))
cutoff_df <- pivot_longer(cutoff_df,cols=c(1:6))
cutoff_df
cutoff_names <- data.frame(lengths=c(2000,2700,2000,2700,2000,2700),value=c(13825,11442,10670,10095,9583,9205),
                           labels=c("1","2","3","4","5","6"))
cutoff_criteria <- ggplot(cutoff_df)+theme_classic()+
  labs(x="Min. DCO distance (bp)",y="Map Length (cM)")+
  geom_line(aes(x=lengths,y=value,group=name))+
  scale_x_continuous(breaks=c(0,2500,5000,7500,10000))+
  geom_label(data=cutoff_names,aes(x=lengths,y=value,label=labels),label.size=0,label.padding=unit(0.1,"lines"),size=3)+
  geom_hline(aes(yintercept=11966),lty=2)


#now we want to extract the genotypes and get the DCO distribution
pr2 <- calc_genoprob(progeny2,error_prob = 0.00)
pr2_geno <- maxmarg(pr2)
pr2_xo_loc <- locate_xo(pr2_geno,progeny2$pmap)
pr2_dco_dist <- unlist(c(lapply(pr2_xo_loc$chr1,diff),lapply(pr2_xo_loc$chr2,diff),lapply(pr2_xo_loc$chr3,diff),lapply(pr2_xo_loc$chr4,diff),
                      lapply(pr2_xo_loc$chr5,diff),lapply(pr2_xo_loc$chr6,diff),lapply(pr2_xo_loc$chr7,diff),lapply(pr2_xo_loc$chr8,diff)))
pr2_shape <- mean(pr2_dco_dist)*mean(pr2_dco_dist)/var(pr2_dco_dist)
pr2_rate = var(pr2_dco_dist)/mean(pr2_dco_dist)


cr_pr_xo_loc <- locate_xo(cr_pr_geno,cr_progeny$pmap)
cr_pr_dco_dist <- unlist(c(lapply(cr_pr_xo_loc$chr1,diff),lapply(cr_pr_xo_loc$chr2,diff),lapply(cr_pr_xo_loc$chr3,diff),lapply(cr_pr_xo_loc$chr4,diff),
                        lapply(cr_pr_xo_loc$chr5,diff),lapply(cr_pr_xo_loc$chr6,diff),lapply(cr_pr_xo_loc$chr7,diff),lapply(cr_pr_xo_loc$chr8,diff)))
cr_pr_shape <- mean(cr_pr_dco_dist)*mean(cr_pr_dco_dist)/var(cr_pr_dco_dist)
cr_pr_rate = var(cr_pr_dco_dist)/mean(cr_pr_dco_dist)

interevent <- ggplot() + theme_classic() + scale_x_continuous(expand=c(0.01,0),breaks=c(0,0.5e6,1e6,1.5e6,2e6),labels=c(0,0.5,1,1.5,2),limits=c(0,2e6)) + scale_y_continuous(expand=c(0,0.01))+ #scale_x_log10(lim=c(1,1e7))+ 
  labs(x="Distance between crossovers\n(Mb)",y="Density")+
  #geom_segment(aes(x=0.7e6,xend=0.9e6,y=0.04,yend=0.04),color="green")+
  #geom_text(aes(x=0.95e6,y=0.04),label=expression("Raw Data; "~italic("v")~"=0.58"),hjust=0,size=2.5)+
  #geom_segment(aes(x=0.7e6,xend=0.9e6,y=0.03,yend=0.03),color="black")+
  #geom_text(aes(x=0.95e6,y=0.03),label=expression(">=2 Markers,3kbp; "~italic("v")~"=0.74"),hjust=0,size=2.5)+
  geom_density(aes(x=pr2_dco_dist,y=..count..),lwd=1,color="green") +
  geom_density(aes(x=cr_pr_dco_dist,y=..count..)) +
  theme(panel.grid.minor=element_blank())

interevent_inset <- ggplot() + theme_classic() + 
  scale_x_continuous(expand=c(0.01,0),breaks=c(0,25000,50000,75000),labels=c("0kb","25kb","50kb","75kb")) + 
  coord_cartesian(xlim=c(0,100000))+
  scale_y_continuous(expand=c(0,0.01))+
  geom_density(aes(x=pr2_dco_dist,y=..count..),lwd=1,color="green") +
  geom_density(aes(x=cr_pr_dco_dist,y=..count..)) +
  theme(axis.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text=element_text(size=8.5))
interevent <- interevent + annotation_custom(ggplotGrob(interevent_inset),xmin=0.42e6,xmax=2.02e6,ymin=0.04,ymax=0.14)
interevent


#make comparison between different species
dat<-read.delim("simplified_data_recombination_distances.txt" , header=F)
colnames(dat) <- c("chr","Mb","cM","species","group")
tidy_dat <- dat %>% group_by(species,group) %>% summarize(mean_cM = mean(cM),
                                                    sum_cM = sum(cM),
                                                    mean_Mb = mean(Mb),
                                                    sum_Mb = sum(Mb))
#add data from other fumigatus manually
new_dat <- read.delim("supplemental_crosses_distances.txt", header=F)
colnames(new_dat) <- c("chr","Mb","cM","species","group")
tidy_new_dat <- new_dat %>% group_by(species,group) %>% summarize(mean_cM = mean(cM),
                                                          sum_cM = sum(cM),
                                                          mean_Mb = mean(Mb),
                                                          sum_Mb = sum(Mb))

whole_gen <- ggplot(tidy_dat) + theme_classic() + 
  geom_point(aes(x=sum_Mb,y=sum_cM,shape=group),color="cornflowerblue",size=2.5) +
  geom_point(data=data.frame(),aes(x=c(-100,-100,-100),y=c(4336,5778,18019)),col="red")+
  geom_text_repel(aes(x=sum_Mb,y=sum_cM,label=species),size=4,min.segment.length=0.4) +
  labs(x="Total Genome (Mb)",y="Map Length (cM)") +
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=9))
whole_gen
per_chr <- ggplot(tidy_dat) + theme_classic() + 
  geom_point(aes(x=mean_Mb,y=mean_cM,shape=group),color="cornflowerblue",size=2.5) +
  #scale_x_log10() + scale_y_log10() +
  geom_point(data=data.frame(),aes(x=c(-20,-20,-20),y=c(4336/8,5778/8,18019/8)),col="red")+
  geom_text_repel(aes(x=mean_Mb,y=mean_cM,label=species),size=4,min.segment.length=0.6)+
  labs(x="Mean Mb per Chromosome",y="Mean cM per\nChromosome")+
  theme(legend.position="none",
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=9))

#here we subsample the dataset of our 195 individuals, to similate teh results of a smaller population

total <- data.frame(indiv = c(NA),
                    length = c(NA))
# we read in the data of the original cross
# and then make 100 subsets and re-estimate the map
orig_cross <- read_cross2("asp_June9_RFadjust.yaml")
orig_new_map <- est_map(orig_cross,tol=0.001,map_function="haldane",cores=4)
for (i in 1:150){
  set.seed(runif(1,1,1e5))
  x <- floor(runif(1,5,80))
  temp_inds <- sample(row.names(orig_cross$pheno),x)
  temp_cross <- subset(orig_cross,ind=temp_inds)
  temp_gmap_sub <- reduce_markers(orig_new_map,20,cores=3)
  temp_cross_markers2keep <- unlist(lapply(temp_gmap_sub,names))
  temp_cross_sub <- pull_markers(temp_cross,temp_cross_markers2keep)
  new_map <- est_map(temp_cross_sub,tol=0.01,cores=3)
  total <- rbind(total,c(x,sum(unlist(lapply(new_map,max)))))
}
total
#the data from below is summarized here
new_crosses <- data.frame(ind = c(12,15,7),length=c(5778,4436,18013))

new_plus_subset <- ggplot() + geom_jitter(aes(x=total$indiv,y=total$length),width=0.1,height=0,col="grey30",alpha=0.7,size=0.5)+
  theme_classic() + geom_hline(aes(yintercept=11966),lty=2) + ylim(0,35000)+
  geom_smooth(aes(x=total$indiv,y=total$length),level=0)+
  geom_point(aes(x=new_crosses$ind,y=new_crosses$length),shape=c(15,15,16),size=4)+
  labs(x="Individuals in Mapping Population",y="Map Length (cM)")


left_4 <- plot_grid(rarified,cutoff_criteria,
                    interevent,D_values_plot,
                    nrow=2,align="hv",axis="lb",
                    labels=c("A","B","C","D"),label_fontface="plain")
left_1 <- plot_grid(new_plus_subset,
                    nrow=1,align="hv",axis="lb",
                    labels=c("E"),label_fontface="plain")
#left_6 <- plot_grid(rarified,cutoff_criteria,
#                    interevent,D_values_plot,
#                    new_plus_subset,new_plus_subset,
#                    nrow=3,align="hv",axis="lb",
#                    labels=c("A","B","C","D","E","F"),label_fontface="plain")
left_6 <- plot_grid(left_4,left_1,nrow=2,rel_heights = c(2,1))
right_2 <- plot_grid(whole_gen,per_chr,
                     nrow=2,align="hv",axis="lb",
                     labels=c("F","G"),label_fontface = "plain")
whole <- plot_grid(left_6,right_2,ncol=2)
whole
svg(paste0("Figure1_",format(Sys.Date(),"%d-%m-%y"),".svg"),width=10,height=6.64)
whole
dev.off()

#now we want to plot the distance between markers physically versus the genetic distance
#first from the progeny2 object, which has gene conversions, then from the cr_progeny object without
t_start <- est_map(progeny2,error_prob=1e-6,tol=1e-2)
t_end <- est_map(cr_progeny,error_prob=1e-6,tol=1e-2)

start_phy <- unlist(lapply(progeny2$pmap,diff))
start_gen <- unlist(lapply(t_start,diff))

end_phy  <- unlist(lapply(cr_progeny$pmap,diff))
end_gen  <- unlist(lapply(t_end,diff))

start_gen_rf <- (1-exp(-2*start_gen))/2
end_gen_rf <- (1-exp(-2*end_gen))/2

starting <- ggplot() + theme_classic()+
  geom_point(aes(x=start_phy,y=start_gen),size=0.1)+
  scale_x_log10(breaks=c(1,10,100,1000,10000))+
  labs(x="",y="Genetic distance \nbetween markers (cM)")+
  ylim(0,75)

ending   <- ggplot() + theme_classic()+
  geom_point(aes(x=end_phy,y=end_gen),size=0.1)+
  scale_x_log10(breaks=c(1,10,100,1000,10000))+
  labs(x="Physical distance between markers (bp)", y="Genetic Distance \nbetween markers (cM)")+
  ylim(0,75)

plot_grid(starting,ending,nrow=2,labels=c("A","B"),label_fontface = "plain")
svg("Supplemental.distance.compare.svg",width=6,height=4)
plot_grid(starting,ending,nrow=2,labels=c("A","B"),label_fontface = "plain")
dev.off()

stapley_data <- read.csv("/Users/ben/Downloads/stapley_data_manual.csv",header=T,stringsAsFactors = F,
                         colClasses = c("character",rep("numeric",5)))

stapley_whole <- ggplot(stapley_data)+theme_classic()+
  scale_x_log10(limits=c(1,max(stapley_data$genome.size)))+
  labs(x="Genome Size (Mb)",y="Genetic Map Length (cM)")+
  geom_point(aes(x=genome.size,y=maplength),size=0.5)+
  geom_point(data=subset(stapley_data,name=="a_fum"), aes(x=genome.size,y=maplength),shape=17)+
  annotate("text",x=500,y=12000,label=expression(italic("A. fumigatus")))

stapley_chr <- ggplot(stapley_data)+theme_classic()+
  scale_x_log10()+
  labs(x="Mean Mb per Chromosome",y="Mean cM per Chromosome")+
  geom_point(aes(x=chr_phys,y=chr_gen),size=0.5)+
  geom_point(data=subset(stapley_data,name=="a_fum"), aes(x=chr_phys,y=chr_gen),shape=17)+
  annotate("text",x=35,y=1500,label=expression(italic("A. fumigatus")))

svg("supplemental_stapley_data.svg",width=6,height=2.8)
plot_grid(stapley_whole,stapley_chr,ncol=2,labels=c("A","B"),label_fontface = "plain")
dev.off()
