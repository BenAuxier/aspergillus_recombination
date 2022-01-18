This set of scripts is used for the analysis of LD decay in Figure 1.

The first script is 02.downloader.sh, which downloads the SRA data from Baber et al., 20201 https://doi.org/10.1038/s41564-021-00993-xf
This download 188 sets of paired end data from BioProject PRJNA697844 

Then, the 03.aligner.sh aligns these reads to the reference Af293 genome, and uses Freebayes for variant calling, which is filtered using a minimum number of missing individuals (AN) as well as a maximum number of heterozygotes (the awk code in the last line).

The file produced, barber.data.vcf is produced.

The 04.LD.analysis.sh script then analyzes this script using plink to calculate LD between markers, additionaly, this script subsets 3 random sets of 5000 variants postions and removes locations information to allow for calculation of genome-wide LD between chromosomes

Finally, the 05.summarize.R script summarizes the large LD plink output files, produces rsq_means.tab and subsampled variants. These are then used for plotting LD decay.

Due to size restrictions, only the first chromosome of our final filtered VCF of our mapping population can be found here, as well as the first 10,000 lines of the population level VCF. The full VCF files are available from ben.auxier@wur.nl currently, and will be put on a repository such as Zenodo following publication.

NOTE: this pipeline is largely inspired by the extremely well documented pileine of Lofgren et al. found at https://github.com/MycoPunk/Afum_PopPan
