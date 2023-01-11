library(qtl2)
library(tidyverse)
library(cowplot)

setwd("/Users/ben/Library/CloudStorage/OneDrive-WageningenUniversity&Research/aspergillus_recombination/new_crosses")
if(Sys.getenv("USER")=="ben"){setwd("/Users/ben/OneDrive\ -\ Wageningen\ University\ &\ Research/aspergillus_recombination/new_crosses")}

#first we want to make some random simulations of less individuals, around 10
set.seed(13335)

#original cross----
orig_cross <- read_cross2("asp_June9_RFadjust.yaml")
orig_new_map <- est_map(orig_cross,tol=0.001,map_function="haldane",cores=4)
orig_cross$gmap <- orig_new_map
sum(unlist(lapply(orig_new_map,max)))

gmap_sub <- reduce_markers(orig_cross$gmap,20)
orig_cross_markers2keep <- unlist(lapply(gmap_sub,names))
orig_cross_sub <- pull_markers(orig_cross,orig_cross_markers2keep)
orig_cross_sub_map <- est_map(orig_cross_sub,tol=0.01,cores=3)
sum(unlist(lapply(orig_cross_sub_map,max)))

pr <- calc_genoprob(orig_cross,error_prob=1e-06)
g <- maxmarg(pr)
ggplot() + geom_histogram(aes(x=rowSums(count_xo(g))))
orig_xo <- locate_xo(g,orig_cross$pmap)

orig_cross_dco_dist <- c(unlist(lapply(c(0,orig_xo$chr1,max(orig_cross$pmap$chr1)),diff)),
                         unlist(lapply(c(0,orig_xo$chr2,max(orig_cross$pmap$chr2)),diff)),
                         unlist(lapply(c(0,orig_xo$chr3,max(orig_cross$pmap$chr3)),diff)),
                         unlist(lapply(c(0,orig_xo$chr4,max(orig_cross$pmap$chr4)),diff)),
                         unlist(lapply(c(0,orig_xo$chr5,max(orig_cross$pmap$chr5)),diff)),
                         unlist(lapply(c(0,orig_xo$chr6,max(orig_cross$pmap$chr6)),diff)),
                         unlist(lapply(c(0,orig_xo$chr7,max(orig_cross$pmap$chr7)),diff)),
                         unlist(lapply(c(0,orig_xo$chr8,max(orig_cross$pmap$chr8)),diff)))

orig_cross_dco_plot <- ggplot() + theme_classic() + geom_density(aes(x=orig_cross_dco_dist)) + scale_x_log10(limits=c(9e1,1e7)) +
  geom_vline(aes(xintercept=c(1e5,1e6)),lty=2) + xlab("Distance between CO in UK data (bp)")

#first cross ----
uk_cross1 <- read_cross2("uk_cross1.yaml")

uk_cross1_new_map <- est_map(uk_cross1, error_prob = 1e-6,tol=0.0001,map_function="haldane")
#uk_cross1$gmap <- uk_cross1_new_map

sum(unlist(lapply(uk_cross1$gmap,max)))

pr <- calc_genoprob(uk_cross1,error_prob=0)
g <- maxmarg(pr)
rowSums(count_xo(g))
uk_cross1_xo <- locate_xo(g,uk_cross1$pmap)

uk_cross1_dco_dist <- c(unlist(lapply(c(0,uk_cross1_xo$chr1,max(uk_cross1$pmap$chr1)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr2,max(uk_cross1$pmap$chr2)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr3,max(uk_cross1$pmap$chr3)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr4,max(uk_cross1$pmap$chr4)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr5,max(uk_cross1$pmap$chr5)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr6,max(uk_cross1$pmap$chr6)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr7,max(uk_cross1$pmap$chr7)),diff)),
                     unlist(lapply(c(0,uk_cross1_xo$chr8,max(uk_cross1$pmap$chr8)),diff)))

sum(diff(uk_cross1$pmap$chr1) > 10000)


uk_cross1_dco_plot <- ggplot() + theme_classic() + geom_histogram(aes(x=uk_cross1_dco_dist)) + scale_x_log10(limits=c(9e1,1e7)) +
  geom_vline(aes(xintercept=c(1e5,1e6)),lty=2) + xlab("Distance between CO in UK data (bp)")

uk_cross1_gmap_sub <- reduce_markers(uk_cross1$gmap,20)
uk_cross1_markers2keep <- unlist(lapply(uk_cross1_gmap_sub,names))
#drop mitochondria, since when only one marker, it leads to errors. No recombination in mitochondira, so produces an error
uk_cross1_markers2keep <- uk_cross1_markers2keep[!grepl("mito",names(uk_cross1_markers2keep))]
uk_cross1_sub <- pull_markers(uk_cross1,uk_cross1_markers2keep)
uk_cross1_sub_map <- est_map(uk_cross1_sub)
sum(unlist(lapply(uk_cross1_sub_map,max)))



#second cross----
uk_cross2 <- read_cross2("uk_cross2.yaml")

uk_cross2_new_map <- est_map(uk_cross2, tol=0.0001,map_function="haldane")

sum(unlist(lapply(uk_cross2_new_map,max)))

pr <- calc_genoprob(uk_cross2,error_prob=0.1)
g <- maxmarg(pr)
ph <- guess_phase(uk_cross2,g)
plot_onegeno(ph,uk_cross2$pmap,ind="4")
rowSums(count_xo(g))
uk_cross2_xo <- locate_xo(g,uk_cross2$pmap)

uk_cross2_dco_dist <- c(unlist(lapply(c(0,uk_cross2_xo$chr1,max(uk_cross2$pmap$chr1)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr2,max(uk_cross2$pmap$chr2)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr3,max(uk_cross2$pmap$chr3)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr4,max(uk_cross2$pmap$chr4)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr5,max(uk_cross2$pmap$chr5)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr6,max(uk_cross2$pmap$chr6)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr7,max(uk_cross2$pmap$chr7)),diff)),
                     unlist(lapply(c(0,uk_cross2_xo$chr8,max(uk_cross2$pmap$chr8)),diff)))

sum(diff(uk_cross2$pmap$chr1) > 50000)

uk_cross2_dco_plot <- ggplot() + theme_classic() + geom_density(aes(x=uk_cross2_dco_dist)) + scale_x_log10(limits=c(1e1,1e7)) +
  geom_vline(aes(xintercept=c(1e5,1e6)),lty=2) + xlab("Distance between CO in UK data (bp)")

uk_cross2_gmap_sub <- reduce_markers(uk_cross2$gmap,20)
uk_cross2_markers2keep <- unlist(lapply(uk_cross2_gmap_sub,names))
#need to purge mitochondria same as for uk_cross1
uk_cross2_markers2keep <- uk_cross2_markers2keep[!grepl("mito",names(uk_cross2_markers2keep))]
uk_cross2_sub <- pull_markers(uk_cross2,uk_cross2_markers2keep)
uk_cross2_sub_map <- est_map(uk_cross2_sub,tol=0.01,cores=3)
sum(unlist(lapply(uk_cross2_sub_map,max)))

#new crosses from baseclear dataset----
C1 <- read_cross2("C1.yaml")

C1_new_map <- est_map(C1,error_prob=0.001, tol=0.0001,map_function="haldane")

pr_C1 <- calc_genoprob(C1,error_prob=1e-4)
g_C1 <- maxmarg(pr_C1)
rowSums(count_xo(g_C1))
xo_C1 <- locate_xo(g_C1,C1$pmap)

C1_dco_dist <- c(unlist(lapply(c(0,xo_C1$chr1,max(C1$pmap$chr1)),diff)),
                 unlist(lapply(c(0,xo_C1$chr2,max(C1$pamp$chr2)),diff)),
                 unlist(lapply(c(0,xo_C1$chr3,max(C1$pamp$chr3)),diff)),
                 unlist(lapply(c(0,xo_C1$chr4,max(C1$pamp$chr4)),diff)),
                 unlist(lapply(c(0,xo_C1$chr5,max(C1$pamp$chr5)),diff)),
                 unlist(lapply(c(0,xo_C1$chr6,max(C1$pamp$chr6)),diff)),
                 unlist(lapply(c(0,xo_C1$chr7,max(C1$pamp$chr7)),diff)),
                 unlist(lapply(c(0,xo_C1$chr8,max(C1$pamp$chr8)),diff)))


C1_dco_plot <- ggplot() + theme_classic() + geom_density(aes(x=C1_dco_dist)) + scale_x_log10(limits=c(1e1,1e7)) +
  geom_vline(aes(xintercept=c(1e5,1e6)),lty=2) + xlab("Distance between CO in UK data (bp)")

sum(unlist(lapply(C1_new_map,max)))

C1_gmap_sub <- reduce_markers(C1$gmap,20)
C1_markers2keep <- unlist(lapply(C1_gmap_sub,names))
#drop mitochondrial markers
C1_markers2keep <- C1_markers2keep[!grepl("mito",names(C1_markers2keep))]
C1_sub <- pull_markers(C1,C1_markers2keep)
C1_sub_map <- est_map(C1_sub,tol=0.01,cores=3)
sum(unlist(lapply(C1_sub_map,max)))


#other new cross from baseclear dataset----
C2 <- read_uk_cross2("C2.yaml")

C2_new_map <- est_map(C2,error_prob=0.001, tol=0.0001,map_function="haldane")
C2$gmap <- C2_new_map
pr_C2 <- calc_genoprob(C2,error_prob=1e-4)
g_C2 <- maxmarg(pr_C2)
phased_C2 <- guess_phase(C2,g_C2)

plot_onegeno(phased_C2,C2$gmap,"2")
xo_C2 <- locate_xo(g_C2,C2$pmap)

C2_dco_dist <- c(unlist(lapply(xo_C2$chr1,diff)),
                 unlist(lapply(xo_C2$chr2,diff)),
                 unlist(lapply(xo_C2$chr3,diff)),
                 unlist(lapply(xo_C2$chr4,diff)),
                 unlist(lapply(xo_C2$chr5,diff)),
                 unlist(lapply(xo_C2$chr6,diff)),
                 unlist(lapply(xo_C2$chr7,diff)),
                 unlist(lapply(xo_C2$chr8,diff)))

ggplot() + theme_classic() + geom_density(aes(x=C2_dco_dist)) + scale_x_log10() +
  geom_vline(aes(xintercept=c(1e5,1e6)),lty=2) + xlab("Distance between CO in UK data (bp)")

sum(unlist(lapply(C2_new_map,max)))

C2_gmap_sub <- reduce_markers(C2$gmap,20)
C2_markers2keep <- unlist(lapply(gmap_sub,names))
C2_sub <- pull_markers(C2,markers2keep)
sub_map <- est_map(C2_sub)
sum(unlist(lapply(sub_map,max)))

plot_grid(orig_cross_dco_plot,uk_cross1_dco_plot,uk_cross2_dco_plot,C1_dco_plot,nrow=4)


unlist(lapply(orig_new_map,max))
unlist(lapply(uk_cross1_new_map,max))
unlist(lapply(uk_cross2_new_map,max))
unlist(lapply(C1_new_map,max))
