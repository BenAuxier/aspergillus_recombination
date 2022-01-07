These scripts were used to assembly the genomes of AfIR964 and AfIR974, the two parental strains

First the upstream_assembly.sh file was run with the Nanopore and Illumina data, then the contigs were scaffolded with the custom_merge.sh script, followed by final racon/pilon polishing.

Following read mapping and variant calling, variants were filtered and converted to a table for input to R/qtl2 using the vcf 2 table.R script
