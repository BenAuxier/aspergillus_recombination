#collect first random sample
grep "#" barber.data.vcf > barber.sample1.vcf
grep -v "#" barber.data.vcf | shuf -n 5000 | awk 'OFS = "\t" {$1 = "ChrX"}{$2 = NR}{print $0}' >> barber.sample1.vcf
tail -n 5 barber.sample1.vcf | cut -b 1-75
/mnt/LTR_userdata/auxie001/programs/plink --vcf barber.sample1.vcf --const-fid --allow-extra-chr --vcf-idspace-to _ --vcf-half-call "h" --keep-allele-order --make-bed --out barber.sample1_LD
/mnt/LTR_userdata/auxie001/programs/plink --bfile "barber.sample1_LD" --r2 gz --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 25 --threads 8 --out "barber.sample1_LD_out"
zcat barber.sample1_LD_out.ld.gz | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > barber.sample1_LD_out.ld.summary


#now collect second random sample
grep "#" barber.data.vcf > barber.sample2.vcf
grep -v "#" barber.data.vcf | shuf -n 5000 | awk 'OFS = "\t" {$1 = "ChrX"}{$2 = NR}{print $0}' >> barber.sample2.vcf
/mnt/LTR_userdata/auxie001/programs/plink --vcf barber.sample2.vcf --const-fid --allow-extra-chr --vcf-idspace-to _ --vcf-half-call "h" --keep-allele-order --make-bed --out barber.sample2_LD
/mnt/LTR_userdata/auxie001/programs/plink --bfile "barber.sample2_LD" --r2 gz --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 25 --threads 8 --out "barber.sample2_LD_out"
zcat barber.sample2_LD_out.ld.gz | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > barber.sample2_LD_out.ld.summary


#get third random sample
grep "#" barber.data.vcf > barber.sample3.vcf
grep -v "#" barber.data.vcf | shuf -n 5000 | awk 'OFS = "\t" {$1 = "ChrX"}{$2 = NR}{print $0}' >> barber.sample3.vcf
/mnt/LTR_userdata/auxie001/programs/plink --vcf barber.sample3.vcf --const-fid --allow-extra-chr --vcf-idspace-to _ --vcf-half-call "h" --keep-allele-order --make-bed --out barber.sample3_LD
/mnt/LTR_userdata/auxie001/programs/plink --bfile "barber.sample3_LD" --r2 gz --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 25 --threads 8 --out "barber.sample3_LD_out"
zcat barber.sample3_LD_out.ld.gz | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > barber.sample3_LD_out.ld.summary


#now the full analysis
/mnt/LTR_userdata/auxie001/programs/plink --vcf barber.data.vcf --const-fid --allow-extra-chr --vcf-idspace-to _ --vcf-half-call "h" --keep-allele-order --make-bed --out barber.data_LD
/mnt/LTR_userdata/auxie001/programs/plink --bfile "barber.data_LD" --r2 gz --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 25 --threads 8 --out "barber.data_LD_out"
zcat barber.data_LD_out.ld.gz | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > barber.data_LD_out.ld.summary


