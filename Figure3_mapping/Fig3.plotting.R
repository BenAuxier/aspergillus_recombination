library(qtl2)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gggenes)
#setwd("C:\\Users\\auxie001\\OneDrive - WageningenUR\\aspergillus_recombination\\qtl")
setwd("/Users/ben/OneDrive\ -\ Wageningen\ University\ &\ Research/aspergillus_recombination/qtl/")
progeny <- read_cross2("asp_June9_DCOadjust_RF.yaml")
pr <- calc_genoprob(progeny,error_prob=0)

growth_50_pheno <- progeny$pheno[,6]
growth_50 <- scan1(pr,growth_50_pheno)
#now run a permutation test on the markers and phenotypes
growth_50_perm <- scan1perm(pr,growth_50_pheno,n_perm=1000,cores=4)
#get the 95% of the permutations
sort(growth_50_perm)[950]
#Permutation based value is 3.93
find_peaks(growth_50,progeny$pmap,threshold=3.93,drop=1.5)

genes <- read.table("V108_20.gff",header=F,skip=11,stringsAsFactors = F)
genes <- genes[genes$V3=="gene",]
genes <- genes[,c(1,4,5,7,9)]
colnames(genes) <- c("chr","start","stop","strand","Name")
genes$start <- genes$start
genes$stop <- genes$stop
genes$Name <- gsub("ID=","",genes$Name)

g50_df <- data.frame(LOD=growth_50,pos=unlist(progeny$pmap)/1e6,markers=names(unlist(progeny$pmap)))
head(g50_df)
colnames(g50_df)[1]<- "LOD"
g50_df <- separate(g50_df,markers,c("chr","marker"),".m")

p1 <- ggplot(data=g50_df)+
  geom_hline(aes(yintercept=3.93),lty=2,lwd=0.5) + 
  geom_rect(data = data.frame(xmin=0.39,xmax=0.85,ymin=3,ymax=12.5,xint=1,chr="chr6"), aes(xmin = xmin,xmax=xmax,ymin=ymin,ymax=ymax),lty = 3,fill=NA,col="black",lwd=0.4)+
  geom_point(aes(x=pos,y=LOD,col=chr),size=1) +
  facet_grid(~chr,scales="free",space="free",switch="both") + 
  theme_classic() + 
  #geom_rect(data=data.frame(),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes = FALSE)+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing = unit(0,"cm"))
p1
p2 <- ggplot(data=subset(g50_df,chr=="chr6"))+
  scale_y_continuous(breaks=c(0,4,8,12))+
  coord_cartesian(xlim=c(0.655,0.675))+
  geom_hline(aes(yintercept=3.93),lty=2)+
  labs(x="Position (Mb)",y="LOD")+
  annotate("text",x=0.662,y=11.9,label="Non-synonymous",size=3)+
  annotate("text",x=0.664,y=11.3,label="Synonymous",size=3)+
  geom_point(aes(x=pos,y=LOD),col="black",fill="cornflowerblue",pch=21,size=3,alpha=0.7) + 
  theme_classic() + 
  theme(legend.position = "none")
p2
#so the region is between 655kb and 675kb on Chr6 from p20_June5, which maps to 666464-688573 from p21_June5

genes <- read.table("Afum_p20_June5.acr.resistance.gff3")
genes <- genes[genes$V3 == "gene",]
genes$V7[genes$V7 == "+"] <- 1
genes$V7[genes$V7 == "-"] <- -1
genes$V7 <- as.numeric(genes$V7)


p2 <- p2 + geom_gene_arrow(data=genes,
                     aes(xmin=V4/1e6,xmax=V5/1e6,y=2,fill=V9,forward=V7),
                     arrow_body_height = grid::unit(8,"mm"),
                     arrowhead_height=grid::unit(8,"mm"),
                     arrowhead_width=grid::unit(6,"mm"))+
  geom_gene_label(data=genes,
                  aes(xmin=V4/1e6,xmax=V5/1e6,y=2,fill=V9,forward=V7,label=V9),min.size=2.5)
p2

top <- plot_grid("none",p2,rel_widths=c(0.3,0.7),labels="AUTO",label_fontface = "plain")
whole <- plot_grid(top,p1,labels=c("","D"),label_fontface="plain",nrow=2)
whole
svg(paste0("acr.resist.",Sys.Date(),".svg"),height=6,width=10)
whole
dev.off()
