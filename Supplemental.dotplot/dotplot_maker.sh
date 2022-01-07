minimap2 -x asm5 -t 4 -X ../genomes/Afum_p20_June5.fna ../genomes/Afum_p21_June5.fna > p20_p21.paf
#minimap2 -x asm5 -t 4 -X ../genomes/Afum_p20_June5.fna ../genomes/Afum_Af293.fna > p20_Af293.paf
#minimap2 -x asm5 -t 4 -X ../genomes/Afum_p21_June5.fna ../genomes/Afum_Af293.fna > p21_Af293.paf


./pafCoordsDotPlotly_updated.R -i p20_p21.paf -o p20_p21 -m 5000 -q 500000 -l
./pafCoordsDotPlotly_updated.R -i p20_Af293.paf -o p20_Af293 -m 5000 -q 500000 -l
./pafCoordsDotPlotly_updated.R -i p21_Af293.paf -o p21_Af293 -m 5000 -q 500000 -l
