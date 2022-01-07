#bwa-mem2 index Afum_Af293.fna

for i in fastq/*1.fastq.gz
do echo $i
echo ""
SAMPLE=${i/_1.fastq.gz/}
echo $SAMPLE
bwa-mem2 mem -t 24 -R "@RG\tID:"$SAMPLE"\tSM:"$SAMPLE Afum_Af293.fna $SAMPLE\_1.fastq.gz $SAMPLE\_2.fastq.gz | samtools view -@ 8 -bh | samtools sort -@ 12 - > ${SAMPLE/fastq/bams}.bam
done

ls bams/*.bam > bamlist.txt

cd bams
for i in *.bam; do samtools index $i; done
cd ..

#old single threaded version
#freebayes -C 15 -Z -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.vcf

#now using parallel
/mnt/LTR_userdata/auxie001/programs/freebayes/scripts/freebayes-parallel <(/mnt/LTR_userdata/auxie001/programs/freebayes/scripts/fasta_generate_regions.py Afum_Af293.fna 200000) 18 -f Afum_Af293.fna -L bamlist.txt | vcffilter -f "QUAL > 1 & SAR > 0 & SAF > 0" > barber.raw.vcf
vcffilter -f "QUAL > 100 & AN > 290" barber.raw.vcf |  awk '{if (gsub(/0\/1/,/0\/1/) < 5) print $0}' | /mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfbiallelic | /mnt/LTR_userdata/auxie001/programs/vcflib/bin/vcfallelicprimitives > barber.data.vcf
