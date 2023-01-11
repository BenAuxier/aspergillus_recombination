library(qtl2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(cowplot)

if(Sys.getenv("USER")=="ben"){setwd("/Users/ben/Library/CloudStorage/OneDrive-WageningenUniversity&Research/aspergillus_recombination/Figure1and2/")}

#now we want to calculate interference across ALL intervals
cr_progeny <- read_cross2("Fig1.lengths/criteria_based_cutoff/asp_June9_criteria_3000_2.yaml")
#we need reasonable distance between markers, so cut to 10cM distance, which should give 1% DCO, or 2 individuals
reduce_map <- reduce_markers(cr_progeny$gmap,min_distance=10)
markers2keep  <- unlist(lapply(reduce_map,names))
working_map <- pull_markers(cr_progeny,markers2keep)
RF1_l <- c()
RF2_l <- c()
RFo_l <- c()
for (chr in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8")){
  markers <- working_map$geno[[chr]]
  for (col in c(3:ncol(markers))){
    print(col)
    temp1 <- markers[,c(col-2,col-1)]
    AA1 <- sum(temp1[,1]==1 & temp1[,2]==1)
    AB1 <- sum(temp1[,1]==1 & temp1[,2]==2)
    BA1 <- sum(temp1[,1]==2 & temp1[,2]==1)
    BB1 <- sum(temp1[,1]==2 & temp1[,2]==2)
    
    temp2 <- markers[,c(col-1,col)]
    AA2 <- sum(temp2[,1]==1 & temp2[,2]==1)
    AB2 <- sum(temp2[,1]==1 & temp2[,2]==2)
    BA2 <- sum(temp2[,1]==2 & temp2[,2]==1)
    BB2 <- sum(temp2[,1]==2 & temp2[,2]==2)
    
    temp_outer <- markers[,c(col-2,col-1,col)]
    AAAo <- sum(temp_outer[,1]==1 & temp_outer[,2]==1 & temp_outer[,3]==1)
    AABo <- sum(temp_outer[,1]==1 & temp_outer[,2]==1 & temp_outer[,3]==2)
    ABAo <- sum(temp_outer[,1]==1 & temp_outer[,2]==2 & temp_outer[,3]==1)
    ABBo <- sum(temp_outer[,1]==1 & temp_outer[,2]==2 & temp_outer[,3]==2)
    BAAo <- sum(temp_outer[,1]==2 & temp_outer[,2]==1 & temp_outer[,3]==1)
    BABo <- sum(temp_outer[,1]==2 & temp_outer[,2]==1 & temp_outer[,3]==2)
    BBAo <- sum(temp_outer[,1]==2 & temp_outer[,2]==2 & temp_outer[,3]==1)
    BBBo <- sum(temp_outer[,1]==2 & temp_outer[,2]==2 & temp_outer[,3]==2)
    
    RF1 <- (AB1+BA1)/sum(c(AA1,AB1,BA1,BB1))
    RF2 <- (AB2+BA2)/sum(c(AA2,AB2,BA2,BB2))
    RFo <- (ABAo+BABo)/sum(c(AAAo,AABo,ABAo,ABBo,BAAo,BABo,BBAo,BBBo))
    RF1_l <- c(RF1_l,RF1)
    RF2_l <- c(RF2_l,RF2)
    RFo_l <- c(RFo_l,RFo)
}
}

coc <- ggplot() + geom_point(aes(x=RF1_l*RF2_l,y=RFo_l),alpha=0.5) + xlim(0,0.15) + ylim(0,0.15) + theme_classic() +
  geom_segment(aes(x=0,y=0,xend=0.15,yend=0.15),lty=2) +
  labs(x=expression("Predicted DCO\n(RF"[left]~" x RF"[right]~")"),y="Observed DCO")


#the chromosomes are 4.66, 4.87,4.03,3.77,3.92,3.86,1.72,1.79
#chr_lengths <- c(4.66,4.87,4.03,3.77,3.92,3.86,1.72,1.79)
chr_lengths <- unlist(lapply(cr_progeny$pmap,max))/1e6
cr_pr <- calc_genoprob(cr_progeny,cr_progeny$gmap,error_prob=0)
cr_g <- maxmarg(cr_pr)
xo_count <- count_xo(cr_g)
#need to drop the parents
xo_count <- xo_count[seq(3,nrow(xo_count)),]
xo_count_norm <- t(t(xo_count)/chr_lengths)

mean(2*unlist(xo_count_norm))
#mean of xo per chr is 3.77, so there are average of 7.55 crossovers per Mb, since offspring only counts half

per_chr <- ggplot() + theme_classic() +
  labs(x="Crossovers\nper Chromosome (per Mb)",y="Count")+
  #geom_histogram(aes(x=2*unlist(xo_count_norm)),bins=150)
  geom_density(aes(x=2*unlist(as.list(xo_count_norm)),y=..count..),color="cornflowerblue",lwd=2.2) +
  geom_density(aes(x=rpois(1568,lambda=7.55),y=..count..),lty=4,lwd=1.2,adjust=2)
per_chr

#here is the comparison of the physical vs. genetic position of chr1
par(mfrow=c(2,4))
plot(cr_progeny$pmap$chr1,cr_progeny$gmap$chr1,main="Chr1",xlab="",ylab="Genetic Position (cM)")
abline(v=4648320*0.55)
abline(a=0,b=0.000397, lty=3)
plot(cr_progeny$pmap$chr2,cr_progeny$gmap$chr2,main="Chr2",xlab="",ylab="")
abline(v=4854928*0.4)
abline(a=0,b=0.000441, lty=3)
plot(cr_progeny$pmap$chr3,cr_progeny$gmap$chr3,main="Chr3",xlab="",ylab="")
abline(v=3994471*0.3)
abline(a=0,b=0.000601, lty=3)
plot(cr_progeny$pmap$chr4,cr_progeny$gmap$chr4,main="Chr4",xlab="",ylab="")
abline(v=3724836*0.32)
abline(a=0,b=0.000483, lty=3)
plot(cr_progeny$pmap$chr5,cr_progeny$gmap$chr5,main="Chr5",xlab="Physical Position (bp)",ylab="Genetic Position (cM)")
abline(v=3899388*0.29)
abline(a=0,b=0.000496, lty=3)
plot(cr_progeny$pmap$chr6,cr_progeny$gmap$chr6,main="Chr6",xlab="Physical Position (bp)",ylab="")
abline(v=3850893*0.3)
abline(a=0,b=0.000508, lty=3)
plot(cr_progeny$pmap$chr7,cr_progeny$gmap$chr7,main="Chr7",xlab="Physical Position (bp)",ylab="")
abline(v=1713723*0.4)
abline(a=0,b=0.000636, lty=3)
plot(cr_progeny$pmap$chr8,cr_progeny$gmap$chr8,main="Chr8",xlab="Physical Position (bp)",ylab="")
abline(v=1780681*0.5)
abline(a=0,b=0.000609, lty=3)
par(mfrow=c(1,1))

#read in 50kb windows from bedtools
windows <- read.table("p20_50kb_windows.bed")
colnames(windows) <- c("chr","p_start","p_end")
windows$g_start <- NA
windows$g_end <- NA

progeny2 <- read_cross2("asp_June9_RFadjust.yaml")
progeny2$gmap <- est_map(progeny2,error_prob=1e-6,tol=0.001)
#now we make a smooth spline fit for the data with GC, and then predict genetic distance for each window position
s_gc_chr1 <- smooth.spline(progeny2$gmap$chr1 ~ progeny2$pmap$chr1,cv=F,spar=1.0)
s_gc_chr2 <- smooth.spline(progeny2$gmap$chr2 ~ progeny2$pmap$chr2,cv=F,spar=1.0)
s_gc_chr3 <- smooth.spline(progeny2$gmap$chr3 ~ progeny2$pmap$chr3,cv=F,spar=1.0)
s_gc_chr4 <- smooth.spline(progeny2$gmap$chr4 ~ progeny2$pmap$chr4,cv=F,spar=1.0)
s_gc_chr5 <- smooth.spline(progeny2$gmap$chr5 ~ progeny2$pmap$chr5,cv=F,spar=1.0)
s_gc_chr6 <- smooth.spline(progeny2$gmap$chr6 ~ progeny2$pmap$chr6,cv=F,spar=1.0)
s_gc_chr7 <- smooth.spline(progeny2$gmap$chr7 ~ progeny2$pmap$chr7,cv=F,spar=1.0)
s_gc_chr8 <- smooth.spline(progeny2$gmap$chr8 ~ progeny2$pmap$chr8,cv=F,spar=1.0)

cr_progeny$gmap <- est_map(cr_progeny,error_prob=1e-6,tol=0.001)
#now we make a smooth spline fit for the data without GC, and then predict genetic distance for each window position
s_chr1 <- smooth.spline(cr_progeny$gmap$chr1 ~ cr_progeny$pmap$chr1,cv=F,spar=1.0)
s_chr2 <- smooth.spline(cr_progeny$gmap$chr2 ~ cr_progeny$pmap$chr2,cv=F,spar=1.0)
s_chr3 <- smooth.spline(cr_progeny$gmap$chr3 ~ cr_progeny$pmap$chr3,cv=F,spar=1.0)
s_chr4 <- smooth.spline(cr_progeny$gmap$chr4 ~ cr_progeny$pmap$chr4,cv=F,spar=1.0)
s_chr5 <- smooth.spline(cr_progeny$gmap$chr5 ~ cr_progeny$pmap$chr5,cv=F,spar=1.0)
s_chr6 <- smooth.spline(cr_progeny$gmap$chr6 ~ cr_progeny$pmap$chr6,cv=F,spar=1.0)
s_chr7 <- smooth.spline(cr_progeny$gmap$chr7 ~ cr_progeny$pmap$chr7,cv=F,spar=1.0)
s_chr8 <- smooth.spline(cr_progeny$gmap$chr8 ~ cr_progeny$pmap$chr8,cv=F,spar=1.0)

for (i in (1:nrow(windows))){
  if (windows$chr[i]=="chr1"){windows$g_start[i]    <- predict(s_chr1,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr1,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr1,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr1,windows$p_end[i])$y}
  if (windows$chr[i]=="chr2"){windows$g_start[i]    <- predict(s_chr2,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr2,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr2,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr2,windows$p_end[i])$y}
  if (windows$chr[i]=="chr3"){windows$g_start[i]    <- predict(s_chr3,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr3,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr3,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr3,windows$p_end[i])$y}
  if (windows$chr[i]=="chr4"){windows$g_start[i]    <- predict(s_chr4,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr4,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr4,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr4,windows$p_end[i])$y}
  if (windows$chr[i]=="chr5"){windows$g_start[i]    <- predict(s_chr5,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr5,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr5,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr5,windows$p_end[i])$y}
  if (windows$chr[i]=="chr6"){windows$g_start[i]    <- predict(s_chr6,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr6,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr6,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr6,windows$p_end[i])$y}
  if (windows$chr[i]=="chr7"){windows$g_start[i]    <- predict(s_chr7,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr7,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr7,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr7,windows$p_end[i])$y}
  if (windows$chr[i]=="chr8"){windows$g_start[i]    <- predict(s_chr8,windows$p_start[i])$y
                              windows$g_end[i]      <- predict(s_chr8,windows$p_end[i])$y
                              windows$g_gc_start[i] <- predict(s_gc_chr8,windows$p_start[i])$y
                              windows$g_gc_end[i]   <- predict(s_gc_chr8,windows$p_end[i])$y}
}
windows$cM <- windows$g_end-windows$g_start
windows$cM[windows$cM < 0] <- 0
windows$gc_cM <- windows$g_gc_end-windows$g_gc_start
windows$gc_cM[windows$gc_cM < 0] <- 0

windows$gc_cont <- windows$gc_cM-windows$cM
write.csv(windows,"windowed.stats.csv")
ggplot() + geom_point(data=windows[windows$chr=="chr7",],aes(x=p_start,y=g_start)) +
  geom_line(aes(x=cr_progeny$pmap$chr7,y=cr_progeny$gmap$chr7))


#now read in the #genes per window
genes_bed <- read.table("p20_50kb_genecounts.bed")


#first we need to make aata frame of centromere positions, I get them from Fedorova 2008, these are the pixel values
centro <- data.frame(chr=names(progeny2$pmap),
                     pos=c(55.91,46.24,31.64,35.32,30.83,32.67,18.98,19.79))
#now to normalize, the chr1 plot is 123.09 wide, which equal 4.663Mb
centro$pos <- centro$pos*(4.6e6/123.09)
#then I need to blast these regions to get the correct centromeres
#2519554.8 1728036.4 1182419.4 1319944.8 1152148.8 1220911.5  709302.1  739572.7
#becomes
centro$pos <- c(1973146,1858544,1158611,1198212,1197926,1332648,768583,849221)

#now we make figure1
chr1_plot <- ggplot() + theme_bw() + scale_x_continuous(expand=c(0,0),limits=c(0,4.9e6))+
  geom_point(aes(x=centro$pos[1],y=0),size=10,color="grey50")+
  geom_rect(aes(xmin=1e4,xmax=centro$pos[1]-1e5,ymin=-50,ymax=50),fill=NA,col="black")+
  geom_rect(aes(xmin=centro$pos[1]+1e5,xmax=chr_lengths[1]*1e6,ymin=-50,ymax=50),fill=NA,col="black")+
  geom_hline(aes(yintercept=mean(windows$cM[windows$chr=="chr1"]/50000*1e6)),lty=3)+
  geom_line(data=windows[windows$chr=="chr1",],aes(x=p_start+25000,y=cM/50000*1e6)) + #scale_y_log10() +
  #geom_line(data=windows[windows$chr=="chr1",],aes(x=p_start+25000,y=-gc_cont/50000*1e6*6))+
  geom_rect(aes(xmin=cr_progeny$pmap$chr1,xmax=cr_progeny$pmap$chr1+2000,ymin=1250,ymax=1300),alpha=0.5) +theme_bw() +
  labs(y="Recomb. (cM/Mb)",x="Position (Mb)")+
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())

chr1_plot

#now we want to plot all the chromosomes


chr_plots <- list()
for (i in seq(1:8)){
  curr_name <- chr_names[i]
  temp_plot <- eval(substitute(ggplot() + theme_bw() + lims(y=c(0,2300))+
    scale_x_continuous(expand=c(0,0),limits=c(0,4.9e6))+
    geom_point(aes(x=centro$pos[[i]],y=0),size=10,color="grey50")+
    geom_rect(aes(xmin=1e4,xmax=2.6e6,ymin=-50,ymax=50),fill=NA,col="black")+
    geom_rect(aes(xmin=2.8e6,xmax=4.8e6,ymin=-50,ymax=50),fill=NA,col="black")+
    geom_hline(aes(yintercept=mean(windows$cM[windows$chr==curr_name]/50000*1e6)),lty=3)+
    geom_line(data=windows[windows$chr==curr_name,],aes(x=p_start+25000,y=cM/50000*1e6)) + #scale_y_log10() +
    #geom_line(data=windows[windows$chr=="chr1",],aes(x=p_start+25000,y=-gc_cont/50000*1e6*6))+
    geom_rect(aes(xmin=cr_progeny$pmap[[i]],xmax=cr_progeny$pmap[[i]]+2000,ymin=2000,ymax=2100),alpha=0.5) +theme_bw() +
    labs(y="Recomb. (cM/Mb)",x="Position (Mb)")+
    theme(panel.grid.minor.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
    ,list(i = i)))
  print(temp_plot)
  chr_plots[[i]] <- temp_plot
}

svg(paste0("Supplemental.Chr.Plots",Sys.Date(),".svg"),width=8,height=12)
plot_grid(chr_plots[[1]],chr_plots[[2]],chr_plots[[3]],chr_plots[[4]],
          chr_plots[[5]],chr_plots[[6]],chr_plots[[7]],chr_plots[[8]],
          ncol=1,labels=chr_names,label_fontface = "plain")
dev.off()


#data is found in xo_count
#we add three more crosses, to compare
#UK1:  61  58  98 105  85  87  61  86  59  56  74  89 
#UK2: 42 64 38 53 46 14 10 39 61 32 54 48 14 43 40
#NL1: 109  94  52 109  72 117  49 
UK1_xo_count <- c(61,58,98,105,85,87,61,86,59,56,74,89)
UK2_xo_count <- c(42,64,38,53,46,14,10,39,61,32,54,48,14,43,40)
NL1_xo_count <- c(109,94,52,109,72,117,49)

#mean is 105.3
rpois(10000,105.3)
xo_count <- as.data.frame(xo_count)
xo_count$sum <- rowSums(xo_count)
sim_xo_dist <- rpois(196,105.3)
hist_xo <- ggplot() + theme_classic()+ labs(y="Density",x="CO Count") +
  geom_segment(aes(x=NL1_xo_count*2,xend=NL1_xo_count*2,y=0.032,yend=0.033),alpha=0.6,linewidth=0.3)+
  geom_segment(aes(x=UK1_xo_count*2,xend=UK1_xo_count*2,y=0.030,yend=0.031),alpha=0.6,linewidth=0.3)+
  geom_segment(aes(x=UK2_xo_count*2,xend=UK2_xo_count*2,y=0.028,yend=0.029),alpha=0.6,linewidth=0.3)+
  geom_text(aes(x=c(275,275,275),y=c(0.0325,0.0305,0.0285),label=c("NL1","UK1","UK2")),vjust="middle",size=2)+
  geom_histogram(aes(x=xo_count$sum*2,y=..density..),bins=40,fill="orange2",col="black") +
  geom_line(aes(x=seq(90,300,1),y=dpois(seq(90,300,1),lambda=210.6)),lty=4,lwd=0.7)
hist_xo

summary(lm(xo_count$chr1~xo_count$sum-xo_count$chr1))
xo_corr <- ggplot(data=xo_count) + theme_classic() + labs(x="# CO on Chr. 1",y="Total # CO\nexcept Chr. 1") +
  geom_point(aes(x=2*chr1,y=2*(sum-chr1)),size=2,alpha=0.6,col="cornflowerblue") +
  stat_smooth(aes(x=2*chr1,y=2*(sum-chr1)),se=0,method="lm",col="black",lty=2)+
  annotate("text",x=40,y=100,label=bquote(R^2~"=0.33"),size=4)+
  xlim(0,55)+ylim(0,300)
xo_corr

#now lets plot some genotypes:
t <- as.data.frame(t(cr_progeny$geno$chr1))[c(1:1269),]
t <- cbind(start = cr_progeny$pmap$chr1[1:length(cr_progeny$pmap$chr1)-1],
           width = diff(cr_progeny$pmap$chr1),
           t)
#131 has the most crossovers, 50 and 173 has the fewest
#170,187,89 are all average, near 103

t <- t[,colnames(t) %in% c("start","width","NIR","CNX","78","171","177","180","217","251")]
colnames(t)[colnames(t)=="NIR"] <- "Afir974"
colnames(t)[colnames(t)=="CNX"] <- "Afir964"
t_long <- t %>% pivot_longer(!c("start","width"),names_to="strain",values_to="geno")
t_long$strain <- factor(t_long$strain,levels=c("Afir974","Afir964","78","171","177","180","217","251"))
t_long$geno <- as.character(t_long$geno)
colors <- c("2" = "cornflowerblue","1"="orange2","0"="grey30")
geno <- ggplot(data=t_long) + theme_void() + 
  geom_rect(aes(xmin=start,xmax=start+width,fill=geno,col=geno,ymin=0,ymax=1)) + facet_grid(rows=vars(strain))+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  geom_segment(aes(x=1e4,xend=max(t$start+t$width),y=0,yend=0)) + geom_segment(aes(x=1e4,xend=max(t$start+t$width),y=1,yend=1))+
  scale_y_continuous(expand=c(0,0)) + scale_x_continuous(expand=c(0,0),limits=c(-3.15e5,max(t$start+t$width)+0.95e5))+
  theme(legend.position="none",
        panel.spacing=unit(0.001,"in"))
geno

top <- plot_grid(geno,chr1_plot,nrow=2,label_fontface="plain",labels=c("A","B"),rel_heights=c(0.6,1))
bottom <- plot_grid(coc,per_chr,hist_xo,xo_corr,nrow=1,labels=c("C","D","E","F"),align="hv")


svg(paste0("Fig2-",Sys.Date(),".svg"),width=10,height=6)
plot_grid(top,bottom,nrow=2,labels="",rel_heights=c(1,0.5))
dev.off()
