### script describing minimap / miniasm and canu / racon pipelines, including pilon polishing, 
### and finally mapping and variant calling, joost.vandenheuvel@wur.nl

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###
###
###  minimap / miniasm /pilon assembly
###
###

# use filtlong to remove reads shorted that 1kb

~/programs/Filtlong/bin/filtlong --min_length 1000 raw.fastq > raw.1kb.fastq

# make map

~/programs/minimap2/minimap2 -x ava-ont -t8  raw.1kb.fastq  raw.1kb.fastq | gzip -1 > reads.1kb.gz

# assembly

~/programs/miniasm/miniasm -f raw.1kb.fastq reads.1kb.gz > reads.1kb.gfa

# make fasta

awk '/^S/{print ">"$2"\n"$3}' reads.1kb.gfa | fold > reads.1kb.fa

# then three rounds of pilon using short reads 

mv reads.1kb.fa fastas/new.fasta

for round in pilon1 pilon2 pilon3
do
mkdir $round
mv fastas/new.changes fastas/$round.changesmadebefore.txt
mv fastas/new.fasta fastas/$round.fasta
~/programs/bwa-0.7.15/bwa index fastas/$round.fasta
~/programs/bwa-0.7.15/bwa mem -t 24 fastas/$round.fasta  ../p0.R1.fq.gz ../p0.R2.fq.gz   | samtools view -Sbq 20  > $round/raw.bam
samtools sort $round/CNX.bam -@ 24 -o $round/sort.bam
samtools index $round/sort.bam
java -jar /mnt/LTR_userdata/heuve063/programs/pilon-1.23.jar --genome fastas/$round.fasta --frags $round/sort.bam --outdir fastas --output new --changes
done

###
###
### end minimap / miniasm / pilon
###
###

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###
###
### canu / racon / pilon
###
###

canu -trim-assemble -p Asp.trim.assem -d Asp.trim.assem genomeSize=30m errorRate=0.05 -nanopore-raw canu/raw.1kb.fastq

mkdir nanopore_mapping

~/programs/bwa-0.7.15/bwa index Asp.trim.assem.contigs.fasta
~/programs/bwa-0.7.15/bwa mem -t 25 -x ont2d Asp.trim.assem.contigs.fasta Asp.trim.assem.trimmedReads.fasta.gz > nanopore_mapping/mapping.sam

mkdir racon

# racon v1.3.1

racon -m 8 -x -6 -g -8 -w 500 -t 25 Asp.trim.assem.trimmedReads.fasta.gz nanopore_mapping/mapping.sam Asp.trim.assem.contigs.fasta > racon/racon.fasta

# pilon v1.23

bwa index racon/racon.fasta

bwa mem -t 16 racon/racon.fasta  ../p0.R1.fq.gz ../p0.R2.fq.gz   | samtools view -Sbq 20  > racon/p0.bam
samtools sort racon/p0.bam -@ 28 -o racon/p0.sort.bam
samtools index racon/p0.sort.bam


java -jar ~/programs/pilon-1.23.jar --genome racon/racon.fasta --frags racon/p0.sort.bam --outdir pilon --output pilon1 --changes

mkdir fastas
cp pilon/pilon1.fasta fastas/new.fasta

for round in pilon2 pilon3 pilon4 pilon5 pilon6 pilon7 pilon8
do
mkdir $round
mv fastas/new.changes fastas/$round.changesmadebefore.txt
mv fastas/new.fasta fastas/$round.fasta
~/programs/bwa-0.7.15/bwa index fastas/$round.fasta
~/programs/bwa-0.7.15/bwa mem -t 24 fastas/$round.fasta  ../p0.R1.fq.gz ../p0.R2.fq.gz   | samtools view -Sbq 20  > $round/p0.bam
samtools sort $round/p0.bam -@ 24 -o $round/p0.sort.bam
samtools index $round/p0.sort.bam
~/programs/pilon-1.23.jar --genome fastas/$round.fasta --frags $round/p0.sort.bam --outdir fastas --output new --changes
done

###
###
### end canu / racon / pilon
###
###

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###
###
### custom.merge.sh used to merge miniasm and canu assemblies
###
###


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 



###
###
###  annotation of the assembled genome
###
###

# masking genome using repeatmasker

~/programs/RepeatMasker/repeatmasker -species "Aspergillus fumigatus" -pa 10 -gff -xsmall  -dir p0.final.pilon p0.final.pilon.fasta


# annotation using augustus. start with the af293 genes, blat and using hints from blatting

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_rna.fna.gz
gunzip GCF_000002655.1_ASM265v1_rna.fna.gz

/programs/miniconda3/pkgs/blat-36-0/bin/blat -t=dna -q=rna -minIdentity=90 -minMatch=1 -minScore=10 p0.masked.fasta GCF_000002655.1_ASM265v1_rna.fna  af293.90.psl

cat af293.90.psl | sort -n -k 16,16 | sort -s -k 14,14 > sorted.af293.90.psl 
~/programs/augustus.2.5.5/scripts/blat2hints.pl --in=sorted.af293.90.psl  --out=NIR.final.hints.gff
export AUGUSTUS_CONFIG_PATH=/programs/augustus.2.5.5/config
~/programs/augustus.2.5.5/bin/augustus --species=aspergillus_fumigatus --gff3=on --extrinsicCfgFile=/programs/augustus.2.5.5/config/extrinsic/extrinsic.E.cfg --hintsfile=p0.final.hints.gff p0.masked.fasta > p0.final.annotated.gff

###
###
### end annotation assembled genome
###
###

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

###
###
### mapping script and variant calling
###
###

## numlist is used, which contain all the sample numbers

while read p 
do
var="D${p}"
~/programs/bwa-0.7.15/bwa mem -t 16 genome/genome.fna rawdat/$var/*.fq.gz | samtools view -Sbq 20 > bams/$p.bam 
samtools sort -@ 16 -o bams/$p.sort.bam bams/$p.bam
rm bams/$p.bam
gatk-4.1.6.0/gatk AddOrReplaceReadGroups -I bams/$p.sort.bam -O bams/$p.sort.rg.bam --RGID $p --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM $p
java -jar ~/bin/picard-2.8.2.jar  MarkDuplicates I=bams/$p.sort.rg.bam O=bams/$p.sort.rg.nodup.bam M=bams/$p.marked_dup_metrics.txt 
rm bams/$p.sort.rg.bam
done < numlist 

freebayes -f genome/genome.fna -L bam.list >vcfs/final.vcf



java -Xmx4g -jar /programs/snpEff/snpEff.jar v108_20.Asp.fum /vcfs/final.vcf > /vcfs/final.annotated.vcf 

###
###
### end mapping and variant calling


