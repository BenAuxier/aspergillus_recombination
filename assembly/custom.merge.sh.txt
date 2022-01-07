### custom.merge.sh, joost.vandenheuvel@wur.nl
### this is more a lab journal but merge two miniasm / minimap genomes with a canu / racon genome, after which they are clean with racon / pilon again. 

###################################

### after assembly and seven pilon round, we come to looking at the genome

#### first most chromosomes are complete, we need to first ditch a part of contig 144
### in local ubuntu blast_P02af293
blastn -db AF293.fasta -query P0.fasta -word_size 500 -outfmt 6 > blast2af293.txt 

### chromosome1 is intact and is tig001
### chromosome2 is in two pieces that starts with tig065 and ends with tig 027. 
###  by comparing tig065 and 027 with asp fum af293 chr 2, they almost 'touch' ,but not completely
### it is likely that the linkage map will show the linkage

### chromosome3 is entirely composed of tig004

### the start of chromosome 4 poses a bit of an issue, where tig 523, 524 and 525 have overlapping regions (seemingly)
### this piece of an blast shows that most of 524 overlaps with the same region as 525
tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007197.1	99.592	44885	129	18	412354	457211	385140	429997	0.0	81824
tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007197.1	97.129	11562	162	114	59565	71048	704816	716285	0.0	19355
tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007197.1	97.910	11628	136	38	41988	53581	704816	716370	0.0	20028

### these regions are extremely similar and map onto the 28s ribosomal protein
### therefore we decide to colapse these sequences and remove 524 from the assembly
### the end of 523 ends with the 18s ribosomal subunit....
### check this.....

### usin the ITS 28s and 18 s cluster for blast

blastn -db P0.fasta -query AB305101.1.fasta -outfmt 6 > rib2P0.txt  

## we indeed get 

#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.867  752     0       1       1       752     38040   37290   0.0     1382          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.734  752     1       1       1       752     14693   13943   0.0     1376          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.468  752     0       4       1       752     45831   45084   0.0     1363          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.468  752     0       4       1       752     53614   52867   0.0     1363          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.336  753     0       4       1       752     30262   29514   0.0     1358          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.335  752     0       5       1       752     22467   21721   0.0     1356          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.202  752     1       5       1       752     69186   68440   0.0     1351          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.936  752     0       8       1       752     61404   60661   0.0     1338          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.803  752     4       5       1       752     6926    6180    0.0     1334          
#AB305101.1      tig00000092_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.547  757     3       6       1       752     76988   76235   0.0     1330          
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.070  753     1       6       1       752     33204   32457   0.0     1347          
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.069  752     0       6       1       752     40975   40231   0.0     1343          
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   95.958  767     9       17      1       752     25444   24685   0.0     1225          
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   94.682  771     15      23      1       752     9777    9014    0.0     1173          
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   94.785  767     17      15      1       752     17633   16875   0.0     1173          
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   92.196  756     20      27      2       752     2050    1329    0.0     1033          
#AB305101.1      tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.805  753     5       4       1       752     490838  490089  0.0     1338          
#AB305101.1      tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.808  755     0       7       1       752     483068  482320  0.0     1336          
#AB305101.1      tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.670  752     1       8       1       752     475302  474560  0.0     1325          
#AB305101.1      tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon   97.368  760     5       13      1       752     498608  497856  0.0     1279          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.404  752     2       8       1       752     42893   42152   0.0     1314          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.143  754     3       10      1       752     50661   49917   0.0     1304          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.011  754     5       9       1       752     19614   18869   0.0     1301          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   97.878  754     3       9       1       752     58550   57808   0.0     1291          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   97.875  753     2       10      1       752     35142   34403   0.0     1290          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   97.742  753     6       8       1       752     27390   26648   0.0     1286          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   92.529  763     28      23      1       752     4028    3284    0.0     1066          
#AB305101.1      tig00000524_pilon_pilon_pilon_pilon_pilon_pilon_pilon   90.513  780     32      34      1       752     11845   11080   0.0     992    

####### so either we stick this together and 

### when we blast this to af293 we get 

#AB305101.1      NC_007197.1     100.000 752     0       0       1       752     442592  441841  0.0     1389 

### so only one is present

### so we use the first one in the 523 contig, and end there
#AB305101.1      tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.670  752     1       8       1       752     475302  474560   
### and then start the next contig (525) at the last position of the last blast hit  
#AB305101.1      tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.069  752     0       6       1       752     40975   40231            

### so we end 523 at 475302 and start 525 at 40975

### between 525 and 037 there is about 8 KB left, we keep these separate

#tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007197.1	99.844	31495	39	5	845768	877252	1560538	1592032	0.0	57880
#tig00000037_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007197.1	99.792	12472	7	3	1	12453	1601214	1613685	0.0	22870

#### chromosome 5 exists in two pieces, which overlap

#tig00000071_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007198.1	99.590	14388	18	21	1553261	1567627	1576474	1590841	0.0	26205
#tig00000032_pilon_pilon_pilon_pilon_pilon_pilon_pilon	NC_007198.1	99.697	75823	132	57	2295079	2370843	1652256	1576474	0.0	1,39E+08
### to see what is happening there we put out two pieces

samtools faidx P0.fasta tig00000071_pilon_pilon_pilon_pilon_pilon_pilon_pilon:1553261-1567627 > overlapping.chr5.txt
samtools faidx P0.fasta tig00000032_pilon_pilon_pilon_pilon_pilon_pilon_pilon:2356843-2370843 > overlapping.rc.chr5.txt

### these two region, as inidicated by the blast output are completely overlapping

### so we keep only till 1553261 of 71 and then 
### keep till 2370843 of 32

### the other chromosomes are fine!

#tig00000001_pilon_pilon_pilon_pilon_pilon_pilon_pilon    -> chr1 - correct order
#tig00000065_pilon_pilon_pilon_pilon_pilon_pilon_pilon    -> chr2a - reverser complement
#tig00000027_pilon_pilon_pilon_pilon_pilon_pilon_pilon    -> chr2b ->reverser complement
#tig00000004_pilon_pilon_pilon_pilon_pilon_pilon_pilon    -> chr3    reverse complement
#tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon (only first 475302 bases [1-475302]) -> chr 4a   correct order
#tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon start at 40975  [40975-  877252] -> chr 4b correct order
#tig00000037_pilon_pilon_pilon_pilon_pilon_pilon_pilon   -> chr4c correct order
#tig00000071_pilon_pilon_pilon_pilon_pilon_pilon_pilon   -> chr 5a [1-1553261]  correct order
#tig00000032_pilon_pilon_pilon_pilon_pilon_pilon_pilon   -> chr 5b [1-2370843] (reverse complement)
#tig00000007_pilon_pilon_pilon_pilon_pilon_pilon_pilon   -> chr6 (reverse complement)
#tig00000060_pilon_pilon_pilon_pilon_pilon_pilon_pilon  -> chr 7 correct order
#tig00000056_pilon_pilon_pilon_pilon_pilon_pilon_pilon  -> chr 8 correct order

mkdir partial_chromosomes

samtools faidx P0.fasta tig00000001_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr1.fasta
samtools faidx P0.fasta tig00000065_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr2a.rc.fasta
samtools faidx P0.fasta tig00000027_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr2b.rc.fasta
samtools faidx P0.fasta tig00000004_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr3.rc.fasta
samtools faidx P0.fasta tig00000523_pilon_pilon_pilon_pilon_pilon_pilon_pilon:1-475302 > partial_chromosomes/P0.chr4a.fasta
samtools faidx P0.fasta tig00000525_pilon_pilon_pilon_pilon_pilon_pilon_pilon:40975-877252 > partial_chromosomes/P0.chr4b.fasta
samtools faidx P0.fasta tig00000037_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr4c.fasta 
samtools faidx P0.fasta tig00000071_pilon_pilon_pilon_pilon_pilon_pilon_pilon:1-1553261 > partial_chromosomes/P0.chr5a.fasta 
samtools faidx P0.fasta tig00000032_pilon_pilon_pilon_pilon_pilon_pilon_pilon:1-2370843 > partial_chromosomes/P0.chr5b.rc.fasta
samtools faidx P0.fasta tig00000007_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr6.rc.fasta                                                    
samtools faidx P0.fasta tig00000060_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr7.fasta  
samtools faidx P0.fasta tig00000056_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P0.chr8.fasta   

#### lastly, the mitohcondrion is also present in several repeats, so lets have a look again

blastn -db P0.fasta -query mito.fasta -outfmt 6 > mito2p0.txt   
### which give us (note that the mitochondrion is 30696 long)
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.722  17638   8       29      1       17626   66819   84427   0.0     32262         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.615  17681   14      34      1       17646   34976   52637   0.0     32225         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.592  17644   11      45      1       17626   3242    20842   0.0     32127         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.508  17677   15      49      1       17646   98586   116221  0.0     32095         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.397  13425   44      24      17286   30696   21574   34975   0.0     24308         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.102  9244    35      13      17286   26516   53342   62550   0.0     16567         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.587  7988    5       18      22713   30696   90622   98585   0.0     14543         
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   98.998  5491    35      15      17286   22773   85159   90632   0.0     9817          
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.860  4274    1       5       26425   30696   62548   66818   0.0     7854          
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   92.996  5126    147     154     20100   25083   120617  125672  0.0     7282          
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.753  3244    1       7       27457   30696   1       3241    0.0     5939          
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   93.995  3747    88      84      17286   20940   116928  120629  0.0     5546          
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   77.912  498     91      17      17659   18149   20944   21429   6.02e-77        292   
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   77.912  498     91      17      17659   18149   84529   85014   6.02e-77        292   
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   77.711  498     90      19      17659   18149   52713   53196   1.01e-74        285   
#NC_017016.1     tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon   77.645  501     88      21      17659   18149   116297  116783  3.62e-74        283


samtools faidx P0.fasta tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon:3242-20842 > partial_chromosomes/P0.mita.fasta  
samtools faidx P0.fasta tig00000144_pilon_pilon_pilon_pilon_pilon_pilon_pilon:21574-34975 > partial_chromosomes/P0.mitb.fasta  

########################## so now in partial_chromosomes, we have to reverse complement a bunch
## on local machine.. 

### for all rc fasta files, we make them reverse complement

cd partial_chromosomes

~/programs/seqtk/seqtk seq -r P0.chr2a.rc.fasta > P0.chr2a.fasta
~/programs/seqtk/seqtk seq -r P0.chr2b.rc.fasta > P0.chr2b.fasta
~/programs/seqtk/seqtk seq -r P0.chr3.rc.fasta > P0.chr3.fasta
~/programs/seqtk/seqtk seq -r P0.chr5b.rc.fasta > P0.chr5b.fasta
~/programs/seqtk/seqtk seq -r P0.chr6.rc.fasta > P0.chr6.fasta

mkdir renamed_chromosomes

for name in chr1 chr2a chr2b chr3 chr4a chr4b chr4c chr5a chr5b chr6 chr7 chr8 mita mitb
do
foo=">${name}"
echo $foo > renamed.$name.fasta
awk 'NR>1 {print $0}' P0.$name.fasta >> renamed.$name.fasta
done

cd ..
mv partial_chromosomes/renamed.* renamed_chromsomes/

#############

### to see how many and potentially which nucleotides are in between the holes of the scaffolds, we can blast the borders of the holes, to 
### the minimap / miniasm assembly

### first we concatenate the genome

cat *.fasta > P0.renamed.total.fasta

samtools faidx P0.renamed.total.fasta

#chr1    4663364 6       60      61                                                                                      
#chr2a   1922009 4741100 1922009 1922010                                                                                 
#chr2b   2942269 6663117 2942269 2942270                                                                                 
#chr3    4026132 9605393 4026132 4026133                                                                                 
#chr4a   475302  13631533        60      61                                                                              
#chr4b   836278  14114764        60      61                                                                              
#chr4c   2416870 14964987        60      61                                                                              
#chr5a   1553261 17422146        60      61                                                                              
#chr5b   2370843 19001302        2370843 2370844                                                                         
#chr6    3857747 21372152        3857747 3857748                                                                         
#chr7    1724474 25229906        60      61                                                                              
#chr8    1792139 26983128        60      61                                                                              
#mita    17601   28805142        60      61                                                                              
#mitb    13402   28823043        60      61                                                                             

samtools faidx P0.renamed.total.fasta chr2a:1872009-1922009 > blast.chr2.min.fasta
samtools faidx P0.renamed.total.fasta chr2b:1-50000 >> blast.chr2.min.fasta  
                                              
samtools faidx P0.renamed.total.fasta chr4a:425302-475302 > blast.chr4.min.fasta
samtools faidx P0.renamed.total.fasta chr4b:1-50000 >> blast.chr4.min.fasta
samtools faidx P0.renamed.total.fasta chr4b:786278-836278 >> blast.chr4.min.fasta
samtools faidx P0.renamed.total.fasta chr4c:1-50000 >> blast.chr4.min.fasta    

samtools faidx P0.renamed.total.fasta chr5a:1503261-1553261 > blast.chr5.min.fasta
samtools faidx P0.renamed.total.fasta chr5b:1-50000 >> blast.chr5.min.fasta

blastn -db P0.minimap.fna -query blast.chr2.min.fasta -outfmt 6 -word_size 1000 > chr2.2.min.txt
### this has as a result
#chr2a:1872009-1922009   utg000007l_pilon_pilon_pilon    99.359  50057   158     120     1       50001   1873549 1923498 0.0     90509                                                                                                           
#chr2b:1-50000   utg000007l_pilon_pilon_pilon    100.000 50000   0       0       1       50000   1927301 1977300 0.0     92333   

### which indicates that indeed we miss a part that is utg000007l_pilon_pilon_pilon:1923498-1927301

blastn -db P0.minimap.fna -query blast.chr4.min.fasta -outfmt 6 -word_size 1000 > chr4.2.min.txt

#chr4a:425302-475302     utg000010l_pilon_pilon_pilon    99.797  39905   43      30      1       39889   91770   51888   0.0     73207                                                                                                           
#chr4a:425302-475302     utg000010l_pilon_pilon_pilon    95.196  4080    81      78      41993   46048   49924   45936   0.0     6342                                                                                                            
#chr4b:1-50000   utg000008l_pilon_pilon_pilon    97.998  34161   273     283     1       33989   23645   57566   0.0     58914                                                                                                                   
#chr4b:1-50000   utg000008l_pilon_pilon_pilon    97.132  16072   203     203     33989   50000   53329   69202   0.0     26888                                                                                                                   
#chr4b:1-50000   utg000010l_pilon_pilon_pilon    96.268  11226   173     201     1       11117   11089   1       0.0     18183                                                                                                                   
#chr4b:786278-836278     utg000008l_pilon_pilon_pilon    99.042  21082   113     62      28943   50001   832984  853999  0.0     37730                                                                                                           
#chr4b:786278-836278     utg000008l_pilon_pilon_pilon    98.836  15117   86      70      12355   27444   816529  831582  0.0     26858                                                                                                           
#chr4b:786278-836278     utg000008l_pilon_pilon_pilon    99.976  12259   0       1       1       12259   803957  816212  0.0     22618                                                                                                           
#chr4c:1-50000   utg000008l_pilon_pilon_pilon    98.925  29012   113     117     19069   48030   882028  910890  0.0     51662                                                                                                                   
#chr4c:1-50000   utg000008l_pilon_pilon_pilon    99.671  18860   34      17      1       18849   863181  882023  0.0     34459                                                                                                                   
#chr4c:1-50000   utg000008l_pilon_pilon_pilon    98.958  1631    7       5       48370   50000   911185  912805  0.0     2909        

### again this nicely matches, we taken 50000- 46048 = 3952 nucleotides lower than 45936 from 10L. 
### this is 45396 - 3952 = 41444 and starting from 11089. 

### between b and c   853999 - 863181 of 8L needs to be added

###########

blastn -db P0.minimap.fna -query blast.chr5.min.fasta -outfmt 6 -word_size 1000 > chr5.2.min.txt

#chr5a:1503261-1553261   utg000004l_pilon_pilon_pilon    99.794  50020   48      34      1       50001   1498378 1548361 0.0     91748                                                                                                           
#chr5b:1-50000   utg000004l_pilon_pilon_pilon    98.973  50139   212     225     1       50000   1548361 1598335 0.0     89456                                                                                                               

### as expected, not bases need to be added......!!!!!


samtools faidx P0.minimap.fna utg000007l_pilon_pilon_pilon:1923498-1927301 > P0.chr2a2.fasta                
samtools faidx P0.minimap.fna utg000010l_pilon_pilon_pilon:11089-41444 > P0.chr4a2.rc.fasta
samtools faidx P0.minimap.fna utg000008l_pilon_pilon_pilon:853999-863181 > P0.chr4b2.fasta

~/programs/seqtk/seqtk seq -r P0.chr4a2.rc.fasta > P0.chr4a2.fasta

for name in chr2a2 chr4a2 chr4b2
do
foo=">${name}"
echo $foo > renamed.$name.fasta
awk 'NR>1 {print $0}' P0.$name.fasta >> renamed.$name.fasta
done

### finalising
################################

### need to fix some of the files

sed 's/in.fasta/renamed.chr1.fasta/g' fix_fasta.py > fix.fasta.chr1.py
sed 's/out.fasta/fixed.chr1.fasta/g' fix.fasta.chr1.py > fix.fasta.chr1.2.py 

sed 's/in.fasta/renamed.chr2a.fasta/g' fix_fasta.py > fix.fasta.chr2a.py
sed 's/out.fasta/fixed.chr2a.fasta/g' fix.fasta.chr2a.py > fix.fasta.chr2a.2.py 

sed 's/in.fasta/renamed.chr2a2.fasta/g' fix_fasta.py > fix.fasta.chr2a2.py
sed 's/out.fasta/fixed.chr2a2.fasta/g' fix.fasta.chr2a2.py > fix.fasta.chr2a2.2.py 

sed 's/in.fasta/renamed.chr2b.fasta/g' fix_fasta.py > fix.fasta.chr2b.py
sed 's/out.fasta/fixed.chr2b.fasta/g' fix.fasta.chr2b.py > fix.fasta.chr2b.2.py 

sed 's/in.fasta/renamed.chr3.fasta/g' fix_fasta.py > fix.fasta.chr3.py
sed 's/out.fasta/fixed.chr3.fasta/g' fix.fasta.chr3.py > fix.fasta.chr3.2.py 

sed 's/in.fasta/renamed.chr4a.fasta/g' fix_fasta.py > fix.fasta.chr4a.py
sed 's/out.fasta/fixed.chr4a.fasta/g' fix.fasta.chr4a.py > fix.fasta.chr4a.2.py 

sed 's/in.fasta/renamed.chr4a2.fasta/g' fix_fasta.py > fix.fasta.chr4a2.py
sed 's/out.fasta/fixed.chr4a2.fasta/g' fix.fasta.chr4a2.py > fix.fasta.chr4a2.2.py 

sed 's/in.fasta/renamed.chr4b.fasta/g' fix_fasta.py > fix.fasta.chr4b.py
sed 's/out.fasta/fixed.chr4b.fasta/g' fix.fasta.chr4b.py > fix.fasta.chr4b.2.py 

sed 's/in.fasta/renamed.chr4b2.fasta/g' fix_fasta.py > fix.fasta.chr4b2.py
sed 's/out.fasta/fixed.chr4b2.fasta/g' fix.fasta.chr4b2.py > fix.fasta.chr4b2.2.py 

sed 's/in.fasta/renamed.chr4c.fasta/g' fix_fasta.py > fix.fasta.chr4c.py
sed 's/out.fasta/fixed.chr4c.fasta/g' fix.fasta.chr4c.py > fix.fasta.chr4c.2.py 

sed 's/in.fasta/renamed.chr5a.fasta/g' fix_fasta.py > fix.fasta.chr5a.py
sed 's/out.fasta/fixed.chr5a.fasta/g' fix.fasta.chr5a.py > fix.fasta.chr5a.2.py 

sed 's/in.fasta/renamed.chr5b.fasta/g' fix_fasta.py > fix.fasta.chr5b.py
sed 's/out.fasta/fixed.chr5b.fasta/g' fix.fasta.chr5b.py > fix.fasta.chr5b.2.py 

sed 's/in.fasta/renamed.chr6.fasta/g' fix_fasta.py > fix.fasta.chr6.py
sed 's/out.fasta/fixed.chr6.fasta/g' fix.fasta.chr6.py > fix.fasta.chr6.2.py 

sed 's/in.fasta/renamed.chr7.fasta/g' fix_fasta.py > fix.fasta.chr7.py
sed 's/out.fasta/fixed.chr7.fasta/g' fix.fasta.chr7.py > fix.fasta.chr7.2.py 

sed 's/in.fasta/renamed.chr8.fasta/g' fix_fasta.py > fix.fasta.chr8.py
sed 's/out.fasta/fixed.chr8.fasta/g' fix.fasta.chr8.py > fix.fasta.chr8.2.py 

sed 's/in.fasta/renamed.mita.fasta/g' fix_fasta.py > fix.fasta.mita.py
sed 's/out.fasta/fixed.mita.fasta/g' fix.fasta.mita.py > fix.fasta.mita.2.py 

sed 's/in.fasta/renamed.mitb.fasta/g' fix_fasta.py > fix.fasta.mitb.py
sed 's/out.fasta/fixed.mitb.fasta/g' fix.fasta.mitb.py > fix.fasta.mitb.2.py 

for pop in chr1 chr2a chr2a2 chr2b chr3 chr4a chr4a2 chr4b chr4b2 chr4c chr5a chr5b chr6 chr7 chr8 mita mitb
do
python fix.fasta.$pop.2.py
rm fix.fasta.$pop.py
rm fix.fasta.$pop.2.py
done
################

#### now we need to give each file of a chromosome the same name and then concat the sequence

mv fixed.chr1.fasta final.chr1.fasta
mv fixed.chr3.fasta final.chr3.fasta
mv fixed.chr6.fasta final.chr6.fasta
mv fixed.chr7.fasta final.chr7.fasta
mv fixed.chr8.fasta final.chr8.fasta

for name in chr2a chr2a2 chr2b
do
echo ">chr2" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr2a.fasta sub.chr2a2.fasta sub.chr2b.fasta > final.chr2.fasta

for name in chr4a chr4a2 chr4b chr4b2 chr4c
do
echo ">chr4" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr4a.fasta sub.chr4a2.fasta sub.chr4b.fasta sub.chr4b2.fasta sub.chr4c.fasta > final.chr4.fasta

for name in chr5a chr5b 
do
echo ">chr5" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr5a.fasta sub.chr5b.fasta > final.chr5.fasta

for name in mita mitb 
do
echo ">mitochondrion" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.mita.fasta sub.mitb.fasta > final.mit.fasta

cat final.chr1.fasta final.chr2.fasta final.chr3.fasta final.chr4.fasta final.chr5.fasta final.chr6.fasta final.chr7.fasta final.chr8.fasta final.mit.fasta > P0.final.total.fasta

cd /Asp.fumigatus/fasta_cleanup/20canu

mkdir final_cleaning

~/programs/bwa-0.7.15/bwa index P0.final.total.fasta
~/programs/bwa-0.7.15/bwa mem -t 25 -x ont2d P0.final.total.fasta Asp.trim.assem/Asp.trim.assem.trimmedReads.fasta.gz > final_cleaning/mapping.sam

cd final_cleaning

racon -m 8 -x -6 -g -8 -w 500 -t 25 ../Asp.trim.assem/Asp.trim.assem.trimmedReads.fasta.gz mapping.sam ../P0.final.total.fasta > racon.fasta

bwa index racon.fasta
bwa mem -t 16 racon.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > P0.bam
samtools sort P0.bam -@ 28 -o P0.sort.bam
samtools index P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome racon.fasta --frags P0.sort.bam --outdir pilon1 --output pilon1 --changes

bwa index pilon1/pilon1.fasta
bwa mem -t 16 pilon1/pilon1.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon1/P0.bam
samtools sort pilon1/P0.bam -@ 28 -o pilon1/P0.sort.bam
samtools index pilon1/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon1/pilon1.fasta --frags pilon1/P0.sort.bam --outdir pilon2 --output pilon2 --changes

bwa index pilon2/pilon2.fasta
bwa mem -t 16 pilon2/pilon2.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon2/P0.bam
samtools sort pilon2/P0.bam -@ 28 -o pilon2/P0.sort.bam
samtools index pilon2/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon2/pilon2.fasta --frags pilon2/P0.sort.bam --outdir pilon3 --output pilon3 --changes

bwa index pilon3/pilon3.fasta
bwa mem -t 16 pilon3/pilon3.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon3/P0.bam
samtools sort pilon3/P0.bam -@ 28 -o pilon3/P0.sort.bam
samtools index pilon3/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon3/pilon3.fasta --frags pilon3/P0.sort.bam --outdir pilon4 --output pilon4 --changes


bwa index pilon4/pilon4.fasta
bwa mem -t 16 pilon4/pilon4.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon4/P0.bam
samtools sort pilon4/P0.bam -@ 28 -o pilon4/P0.sort.bam
samtools index pilon4/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon4/pilon4.fasta --frags pilon4/P0.sort.bam --outdir pilon5 --output pilon5 --changes

bwa index pilon5/pilon5.fasta
bwa mem -t 16 pilon5/pilon5.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon5/P0.bam
samtools sort pilon5/P0.bam -@ 28 -o pilon5/P0.sort.bam
samtools index pilon5/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon5/pilon5.fasta --frags pilon5/P0.sort.bam --outdir pilon6 --output pilon6 --changes

bwa index pilon6/pilon6.fasta
bwa mem -t 16 pilon6/pilon6.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon6/P0.bam
samtools sort pilon6/P0.bam -@ 28 -o pilon6/P0.sort.bam
samtools index pilon6/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon6/pilon6.fasta --frags pilon6/P0.sort.bam --outdir pilon7 --output pilon7 --changes

bwa index pilon7/pilon7.fasta
bwa mem -t 16 pilon7/pilon7.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon7/P0.bam
samtools sort pilon7/P0.bam -@ 28 -o pilon7/P0.sort.bam
samtools index pilon7/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon7/pilon7.fasta --frags pilon7/P0.sort.bam --outdir pilon8 --output pilon8 --changes

bwa index pilon8/pilon8.fasta
bwa mem -t 16 pilon8/pilon8.fasta  ../P0.R1.fq.gz ../P0.R2.fq.gz   | samtools view -Sbq 20  > pilon8/P0.bam
samtools sort pilon8/P0.bam -@ 28 -o pilon8/P0.sort.bam
samtools index pilon8/P0.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon8/pilon8.fasta --frags pilon8/P0.sort.bam --outdir pilon9 --output pilon9 --changes

##############################
###############################
###########################

cd /programs/BUSCO-docker
cp /Asp.fumigatus/fasta_cleanup/20canu/final_cleaning/pilon5/pilon5.fasta pilon5.final_cleaning.fasta

sudo docker run -it --rm -v $(pwd):/home/working -w /home/working chrishah/busco-docker run_BUSCO.py \
--in pilon5.final_cleaning.fasta  --out pilon5.final_cleaning.busco -l ./eurotiomycetes_odb9 --mode genome --species aspergillus_fumigatus 



# BUSCO version is: 3.1.0                                                                                                                                                            # The lineage dataset is: eurotiomycetes_odb9 (Creation date: 2016-02-13, number of species: 25, number of BUSCOs: 4046)                                                             
# To reproduce this run: python /usr/bin/run_BUSCO.py -i pilon5.final_cleaning.fasta -o pilon5.final_cleaning.busco -l ./eurotiomycetes_odb9/ -m genome -c 1 -sp aspergillus_fumigatus                                                                                                                                                                                    
#                                                                                                                                                                                    
# Summarized benchmarking in BUSCO notation for file pilon5.final_cleaning.fasta                                        
# BUSCO was run in mode: genome                                                                                                                                                                                                                                                                                                                                                   
C:98.9%[S:98.8%,D:0.1%],F:0.5%,M:0.6%,n:4046                                                                                                                                                                                                                                                                                                                              
4000    Complete BUSCOs (C)                                                                                                                                                          
3996    Complete and single-copy BUSCOs (S)                                                                                                                                          
4       Complete and duplicated BUSCOs (D)                                                                                                                                           
21      Fragmented BUSCOs (F)                                                                                                                                                        
25      Missing BUSCOs (M)                                                                                                                                                           
4046    Total BUSCO groups searched                                                                                                                                          
~                                                     


### pilon5 seems like the right genome## now we copy it to the 20canu directory and start making a masked version

/Asp.fumigatus/fasta_cleanup/20canu$ cp final_cleaning/pilon5/pilon5.fasta P0.final.pilon.fasta    


/programs/RepeatMasker/repeatmasker -species "Aspergillus fumigatus" -pa 10 -gff -xsmall  -dir P0.final.pilon P0.final.pilon.fasta

vim names.txt
#>chr1
#>chr2
#>chr3
#>chr4
#>chr5
#>chr6
#>chr7
#>chr8
#>mitochondrion

cd /Asp.fumigatus/fasta_cleanup/20canu/P0.final.pilon

perl /programs/rename_fasta_headers.pl names.txt P0.final.pilon.fasta.masked > P0.masked.fasta

export AUGUSTUS_CONFIG_PATH=/programs/augustus.2.5.5/config
/programs/augustus.2.5.5/bin/augustus --species=aspergillus_fumigatus P0.masked.fasta > P0.wout.hints.txt

head -10 P0.wout.hints.txt > P0.final.wout.gff
awk '$2=="AUGUSTUS" {print $0}' P0.wout.hints.txt >> P0.final.wout.gff 

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_rna.fna.gz
gunzip GCF_000002655.1_ASM265v1_rna.fna.gz

/programs/miniconda3/pkgs/blat-36-0/bin/blat -t=dna -q=rna -minIdentity=90 -minMatch=1 -minScore=10 P0.masked.fasta GCF_000002655.1_ASM265v1_rna.fna  af293.90.psl

cat af293.90.psl | sort -n -k 16,16 | sort -s -k 14,14 > sorted.af293.90.psl 
/programs/augustus.2.5.5/scripts/blat2hints.pl --in=sorted.af293.90.psl  --out=P0.final.hints.gff
export AUGUSTUS_CONFIG_PATH=/programs/augustus.2.5.5/config
/programs/augustus.2.5.5/bin/augustus --species=aspergillus_fumigatus --gff3=on --extrinsicCfgFile=/programs/augustus.2.5.5/config/extrinsic/extrinsic.E.cfg --hintsfile=P0.final.hints.gff P0.masked.fasta > P0.final.annotated.gff

head -10 V108_21.final.annotated.gff > V108_21.final.annotated.clean.gff
awk '$2=="AUGUSTUS" {print $0}' V108_21.final.annotated.gff >> V108_21.final.annotated.clean.gff


### snpeff


cd /programs/snpEff/data/genomes
cp /Asp.fumigatus/fasta_cleanup/20canu/P0.final.pilon/P0.masked.fasta v108_20.Asp.fum.fa



cd /programs/snpEff/data
mkdir v108_20.Asp.fum
cd v108_20.Asp.fum
cp /Asp.fumigatus/fasta_cleanup/20canu/P0.final.pilon/V108_20.final.annotated.gff genes.gff
~/bin/nobackup/gffread/gffread genes.gff -T -o genes.gtf


cd /programs/snpEff
### change configure file at 

vim snpEff.config

## add the following

#aspergillus v108_20                                                                                                    
v108_20.Asp.fum.genome : apsergillus fumigatus v108 20 de novo assembly                                                                                                    

## standard esc :wq

### then we build the annotation

cd /programs/snpEff
java -jar snpEff.jar build -gff3 -v v108_20.Asp.fum

### now 'v108_20.Asp.fum' can be added as a species for which annotation is done

java -Xmx4g -jar /programs/snpEff/snpEff.jar v108_20.Asp.fum /Asp.fum.off/vcfs/final.vcf > /Asp.fum.off/vcfs/final.annotated.vcf 



#####################     ####  #   #  #   #
####################     #      ##  #   # #
###################      #      # # #    #
##################       #      #  ##   # #
#####################     ####  #   #  #   #


###################################

cd /Asp.fumigatus/fasta_cleanup/21canu/Asp.trim.assem/

mkdir finalising

cp fastas/pilon8.fasta finalising/P1.fasta

cd /Asp.fumigatus/fasta_cleanup/21canu/Asp.trim.assem/fastas/

### after assembly and seven pilon round, we come to looking at the genome

#### first most chromosomes are complete, we need to first ditch a part of contig 144
### in local ubuntu blast_P02af293

### this will also be done locally

## ingredients are af293 (GCF_000002655.1_ASM265v1_genomic.fna), the minion assembly after two pilon rounds  21.1.fna and the canu assembly after 7 rounds of pilon
###(P1.fasta) and mito.fasta (whcih is the mitochondrion)

blastn -db af293.fna -query P1.fasta -word_size 500 -outfmt 6 > P12af293.txt

blastn -db 21.1.fna -query P1.fasta -word_size 500 -outfmt 6 > min2P1.txt 

samtools faidx 21.1.fna
samtools faidx af293.fna
samtools faidx P1.fasta

### af293.fna.fai
#NC_007194.1     4918979 85      80      81            aka chr1                                                                                                                           
#NC_007195.1     4844472 4980637 80      81            aka chr2                                                                                                                     
#NC_007196.1     4079167 9885750 80      81            aka chr3                                                                                                             
#NC_007197.1     3923705 14015992        80      81    aka chr4                                                                                                            
#NC_007198.1     3948441 17988829        80      81    aka chr5                                                                                                            
#NC_007199.1     3778736 21986711        80      81    aka chr6                                                                                                                    
#NC_007200.1     2058334 25812767        80      81    aka chr7                                                                                                             
#NC_007201.1     1833124 27896916        80      81    aka chr8


### 21.1.fna.fai
#utg000001l_pilon_pilon  4567384 24      80      81                                                                                                                        
#utg000002l_pilon_pilon  3720278 4624525 80      81                                                                                                                        
#utg000003l_pilon_pilon  4651496 8391331 80      81                                                                                                                        
#utg000004l_pilon_pilon  1751612 13100995        80      81                                                                                                                
#utg000005l_pilon_pilon  1677177 14874527        80      81                                                                                                                
#utg000006l_pilon_pilon  702959  16572693        80      81                                                                                                                
#utg000007l_pilon_pilon  1728523 17284463        80      81                                                                                                                
#utg000008l_pilon_pilon  2148029 19034617        80      81                                                                                                                
#utg000009l_pilon_pilon  2336490 21209521        80      81                                                                                                                
#utg000010l_pilon_pilon  1622022 23575242        80      81                                                                                                                
#utg000011l_pilon_pilon  745732  25217564        80      81                                                                                                                
#utg000012l_pilon_pilon  26745   25972642        80      81                                                                                                                
#utg000013l_pilon_pilon  30609   25999746        80      81                                                                                                                
#utg000014l_pilon_pilon  31330   26030762        80      81                                                                                                                
#utg000015l_pilon_pilon  717415  26062508        80      81                                                                                                                
#utg000016l_pilon_pilon  884262  26788915        80      81                                                                                                                
#utg000017l_pilon_pilon  957728  27684255        80      81                                                                                                                
#utg000018l_pilon_pilon  31052   28653979        80      81                                                                                                                
#utg000019l_pilon_pilon  33038   28685444        80      81                                                                                                                
#utg000020l_pilon_pilon  35515   28718919        80      81                                                                                                                
#utg000021l_pilon_pilon  31072   28754902        80      81                                                                                                                
#utg000022l_pilon_pilon  29672   28786387        80      81                                                                                                                
#utg000023l_pilon_pilon  34566   28816454        80      81                                                                                                                
#utg000024l_pilon_pilon  50838   28851477        80      81                                                                                                                
#utg000025l_pilon_pilon  32236   28902975        80      81                                                                                                                
#utg000026l_pilon_pilon  31548   28935638        80      81                                                                                                                
#utg000027l_pilon_pilon  29090   28967605        80      81                                                                                                                
#utg000028c_pilon_pilon  15472   28997083        80      81                                                                                                                
#utg000029l_pilon_pilon  295360  29012773        80      81                                                                                                                
#utg000030l_pilon_pilon  44463   29311849        80      81                                                                                                                
#utg000031l_pilon_pilon  33613   29356892        80      81                                                                                                                
#utg000032l_pilon_pilon  31450   29390950        80      81                                                                                                                
#utg000033l_pilon_pilon  39429   29422818        80      81   

######## P1.fasta.fai
#tig00000001_pilon_pilon_pilon_pilon_pilon_pilon_pilon   4005027 55      80      81                                                                                        
#tig00000006_pilon_pilon_pilon_pilon_pilon_pilon_pilon   3794862 4055200 80      81                                                                                        
#tig00000020_pilon_pilon_pilon_pilon_pilon_pilon_pilon   46664   7897553 80      81                                                                                        
#tig00000045_pilon_pilon_pilon_pilon_pilon_pilon_pilon   2122923 7944856 80      81                                                                                        
#tig00000050_pilon_pilon_pilon_pilon_pilon_pilon_pilon   2016544 10094371        80      81                                                                                
#tig00000058_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1743926 12136177        80      81                                                                                
#tig00000070_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1429807 13901958        80      81                                                                                
#tig00000072_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1434291 15349693        80      81                                                                                
#tig00000073_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1510027 16801968        80      81                                                                                
#tig00000088_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1755998 18330926        80      81                                                                                
#tig00000096_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1303808 20108929        80      81                                                                                
#tig00000113_pilon_pilon_pilon_pilon_pilon_pilon_pilon   11641   21429090        80      81                                                                                
#tig00000117_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1253769 21440932        80      81                                                                                
#tig00000121_pilon_pilon_pilon_pilon_pilon_pilon_pilon   1343789 22710429        80      81                                                                                
#tig00000136_pilon_pilon_pilon_pilon_pilon_pilon_pilon   943784  24071071        80      81                                                                                
#tig00000163_pilon_pilon_pilon_pilon_pilon_pilon_pilon   983248  25026708        80      81                                                                                
#tig00000178_pilon_pilon_pilon_pilon_pilon_pilon_pilon   659668  26022302        80      81                                                                                
#tig00000195_pilon_pilon_pilon_pilon_pilon_pilon_pilon   700958  26690271        80      81                                                                                
#tig00000208_pilon_pilon_pilon_pilon_pilon_pilon_pilon   510470  27400046        80      81                                                                                
#tig00000215_pilon_pilon_pilon_pilon_pilon_pilon_pilon   615440  27916952        80      81                                                                                
#tig00000263_pilon_pilon_pilon_pilon_pilon_pilon_pilon   352068  28540140        80      81                                                                                
#tig00000299_pilon_pilon_pilon_pilon_pilon_pilon_pilon   47543   28896664        80      81                                                                                
#tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   107905  28944857        80      81     

################
## first we map the tigs to the af293 and recall the chromosome names

awk '$4>50000 {print $0}' P12af293.txt > fil.P1.af293.txt

### chr1 (all in correct order)
#tig178     NC_007194.1 395539-889830
#tig050     NC_007194.1 956080-2887017
#tig163     NC_007194.1 2939439- 3872731
#tig136     NC_007194.1 3974323-4807735

### chr2
#tig072     NC_007195.1  221676- 1370995 (rev comp)
#tig045     NC_007195.1  1431930- 3501464
#tig096     NC_007195.1   3514625 - 4805440

### chr3
#tig001    NC_007196.1 83227 - 3986759

### chr4
#tig208  NC_007197.1  124441 - 299092 
#tig263  NC_007197.1  754303 - 974252 (rev comp)
#tig073  NC_007197.1  1015604 - 2540893
#tig070  NC_007197.1  2640343 - 3909119

### chr 5 
#tig117    NC_007198.1 307212 - 1197693 (rev comp)
#tig195    NC_007198.1 1474801 - 1951765
#tig215    NC_007198.1 2340156 - 2538633  (rev comp)
#tig121    NC_007198.1 2650717 - 3821307

#### chr 6
#tig006   NC_007199.1 109639 - 3708038

### chr 7
#tig088 NC_007200.1 130852 - 1620423

### chr8 
#tig058 71249 - 1630154

#### now make partial chromosome files
### chr1
samtools faidx P1.fasta tig00000178_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr1a.fasta
samtools faidx P1.fasta tig00000050_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr1b.fasta
samtools faidx P1.fasta tig00000163_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr1c.fasta
samtools faidx P1.fasta tig00000136_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr1d.fasta

###chr2

samtools faidx P1.fasta tig00000072_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr2a.rc.fasta
samtools faidx P1.fasta tig00000045_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr2b.fasta
samtools faidx P1.fasta tig00000096_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr2c.fasta

### chr3
samtools faidx P1.fasta tig00000001_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr3.fasta

### chr4
samtools faidx P1.fasta tig00000208_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr4a.fasta
samtools faidx P1.fasta tig00000263_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr4b.rc.fasta
samtools faidx P1.fasta tig00000073_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr4c.fasta
samtools faidx P1.fasta tig00000070_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr4d.fasta 
### chr5
samtools faidx P1.fasta tig00000117_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr5a.rc.fasta
samtools faidx P1.fasta tig00000195_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr5b.fasta
samtools faidx P1.fasta tig00000215_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr5c.rc.fasta
samtools faidx P1.fasta tig00000121_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr5d.fasta  
## chr6,7,8
samtools faidx P1.fasta tig00000006_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr6.fasta
samtools faidx P1.fasta tig00000088_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr7.fasta
samtools faidx P1.fasta tig00000058_pilon_pilon_pilon_pilon_pilon_pilon_pilon > partial_chromosomes/P1.chr8.fasta

cp ../blast_P02af293/mito.fasta ../blastP12/mito.fasta  

blastn -db P1.fasta -query mito.fasta -outfmt 6 > mito2P0.txt   

#NC_017016.1     tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.683  17639   14      31      1       17626   68255   50646   0.0     32225             
#NC_017016.1     tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.665  17634   14      33      1       17626   36517   18921   0.0     32195             
#NC_017016.1     tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.596  15837   8       40      1805    17626   98165   82370   0.0     28840             
#NC_017016.1     tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.538  13416   37      14      17286   30696   49913   36518   0.0     24408             
#NC_017016.1     tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.396  13420   38      30      17286   30696   81641   68256   0.0     24293             
#NC_017016.1     tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon   99.262  13418   40      39      17286   30696   18195   4830    0.0     24175     

samtools faidx P1.fasta tig00000317_pilon_pilon_pilon_pilon_pilon_pilon_pilon:18921-49913 > partial_chromosomes/mito.rc.fasta

################################################################################
###############################################################################


cd partial_chromosomes

~/programs/seqtk/seqtk seq -r P1.chr2a.rc.fasta > P1.chr2a.fasta
~/programs/seqtk/seqtk seq -r P1.chr4b.rc.fasta > P1.chr4b.fasta
~/programs/seqtk/seqtk seq -r P1.chr5a.rc.fasta > P1.chr5a.fasta
~/programs/seqtk/seqtk seq -r P1.chr5c.rc.fasta > P1.chr5c.fasta
~/programs/seqtk/seqtk seq -r mito.rc.fasta > P1.mito.fasta



for name in chr1a chr1b chr1c chr1d chr2a chr2b chr2c chr3 chr4a chr4b chr4c chr4d chr5a chr5b chr5c chr5d chr6 chr7 chr8 mito
do
foo=">${name}"
echo $foo > renamed_chromosomes/renamed.$name.fasta
awk 'NR>1 {print $0}' partial_chromosomes/P1.$name.fasta >> renamed_chromosomes/renamed.$name.fasta
done


#############

### to see how many and potentially which nucleotides are in between the holes of the scaffolds, we can blast the borders of the holes, to 
### the minimap / miniasm assembly

### first we concatenate the genome

cat *.fasta > P1.renamed.total.fasta

############################################

##################### now we blast this against the minimap / miniasm genome

blastn -db P1.minimap.fna -query P1.renamed.total.fasta -outfmt 6 -word_size 1000 > renamedP1.2.min.txt

blastn -db P1.minimap.fna -query P1.renamed.total.fasta -outfmt 6 -word_size 100 > renamedP1.2.100.min.txt

## samtools faidx P1.renamed.total.fasta
#chr1a   659668  7       60      61                                                                                      
#chr1b   2016544 670677  60      61                                                                                      
#chr1c   983248  2720838 60      61                                                                                      
#chr1d   943784  3720481 60      61                                                                                      
#chr2a   1434291 4680002 1434291 1434292                                                                                 
#chr2b   2122923 6114301 60      61                                                                                      
#chr2c   1303808 8272614 60      61                                                                                      
#chr3    4005027 9598159 60      61                                                                                      
#chr4a   510470  13669944        60      61                                                                              
#chr4b   352068  14188929        352068  352069                                                                          
#chr4c   1510027 14541005        60      61                                                                              
#chr4d   1429807 16076207        60      61                                                                              
#chr5a   1253769 17529852        1253769 1253770                                                                         
#chr5b   700958  18783629        60      61                                                                              
#chr5c   615440  19496277        615440  615441                                                                          
#chr5d   1343789 20111725        60      61                                                                              
#chr6    3794862 21477917        60      61                                                                              
#chr7    1755998 25336033        60      61                                                                              
#chr8    1743926 27121304        60      61                                                                              
#mito    30993   28894302        30993   30994    

## samtools P1.minimap.fna
#utg000001l_pilon_pilon  4567384 24      80      81                                                                      
#utg000002l_pilon_pilon  3720278 4624525 80      81                                                                      
#utg000003l_pilon_pilon  4651496 8391331 80      81                                                                      
#utg000004l_pilon_pilon  1751612 13100995        80      81                                                              
#utg000005l_pilon_pilon  1677177 14874527        80      81                                                              
#utg000006l_pilon_pilon  702959  16572693        80      81                                                              
#utg000007l_pilon_pilon  1728523 17284463        80      81                                                              
#utg000008l_pilon_pilon  2148029 19034617        80      81                                                              
#utg000009l_pilon_pilon  2336490 21209521        80      81                                                              
#utg000010l_pilon_pilon  1622022 23575242        80      81                                                              
#utg000011l_pilon_pilon  745732  25217564        80      81                                                              
#utg000012l_pilon_pilon  26745   25972642        80      81                                                              
#utg000013l_pilon_pilon  30609   25999746        80      81                                                              
#utg000014l_pilon_pilon  31330   26030762        80      81                                                              
#utg000015l_pilon_pilon  717415  26062508        80      81                                                              
#utg000016l_pilon_pilon  884262  26788915        80      81                                                              
#utg000017l_pilon_pilon  957728  27684255        80      81                                                              
#utg000018l_pilon_pilon  31052   28653979        80      81                                                              
#utg000019l_pilon_pilon  33038   28685444        80      81                                                              
#utg000020l_pilon_pilon  35515   28718919        80      81                                                              
#utg000021l_pilon_pilon  31072   28754902        80      81                                                              
#utg000022l_pilon_pilon  29672   28786387        80      81                                                              
#utg000023l_pilon_pilon  34566   28816454        80      81                                                              
#utg000024l_pilon_pilon  50838   28851477        80      81                                                              
#utg000025l_pilon_pilon  32236   28902975        80      81                                                              
#utg000026l_pilon_pilon  31548   28935638        80      81                                                              
#utg000027l_pilon_pilon  29090   28967605        80      81                                                              
#utg000028c_pilon_pilon  15472   28997083        80      81                                                              
#utg000029l_pilon_pilon  295360  29012773        80      81                                                              
#utg000030l_pilon_pilon  44463   29311849        80      81                                                              
#utg000031l_pilon_pilon  33613   29356892        80      81                                                              
#utg000032l_pilon_pilon  31450   29390950        80      81                                                              
#utg000033l_pilon_pilon  39429   29422818        80      81   

### now lets look at chr 1

#chr1a   659668  7       60      61                                                                                      
#chr1b   2016544 670677  60      61                                                                                      
#chr1c   983248  2720838 60      61                                                                                      
#chr1d   943784  3720481 60      61    


##  3l 939965 - 969746 is missing
#chr1d	utg000003l_pilon_pilon	98.417	56969	390	354	1	56793	939965	883333	0.0	99733
#chr1c	utg000003l_pilon_pilon	97.449	150747	1717	1577	832151	982194	1119066	969746	0.0	2,55E+08
## 3l  1948849 - 1959391 is missing
#chr1c	utg000003l_pilon_pilon	99.430	16133	38	37	1	16105	1948849	1932743	0.0	29233
#chr1b	utg000003l_pilon_pilon	97.291	146211	1700	1655	1871081	2016544	2104087	1959391	0.0	2,46E+08
#   3l 3968539 - 3993534
#chr1b	utg000003l_pilon_pilon	99.060	165502	679	642	1	165185	3968539	3803597	0.0	2,96E+08
#chr1a	utg000003l_pilon_pilon	99.545	188842	359	357	471048	659668	4182095	3993534	0.0	3,44E+08

samtools faidx P1.renamed.total.fasta chr1a > filtered_chromosomes/chr1a.fasta
samtools faidx P1.renamed.total.fasta chr1b > filtered_chromosomes/chr1b.fasta
samtools faidx P1.renamed.total.fasta chr1c > filtered_chromosomes/chr1c.fasta
samtools faidx P1.renamed.total.fasta chr1d > filtered_chromosomes/chr1d.fasta

samtools faidx P1.minimap.fna utg000003l_pilon_pilon:3968539-3993534 > filtered_chromosomes/P1.chr1a2.rc.fasta
samtools faidx P1.minimap.fna utg000003l_pilon_pilon:1948849-1959391 > filtered_chromosomes/P1.chr1b2.rc.fasta
samtools faidx P1.minimap.fna utg000003l_pilon_pilon:939965-968692 > filtered_chromosomes/P1.chr1c2.rc.fasta

#################

#chr2a   1434291 4680002 1434291 1434292                                                                                 
#chr2b   2122923 6114301 60      61                                                                                      
#chr2c   1303808 8272614 60      61  


#chr2c	utg000001l_pilon_pilon	98.566	205153	1240	1236	1	204550	1297166	1093113	0.0	3,61E+08
#chr2b	utg000001l_pilon_pilon	97.129	31558	414	362	2087075	2118484	1341517	1310304	0.0	52804

#chr2b	utg000001l_pilon_pilon	99.376	48719	115	120	1	48633	3420677	3372062	0.0	8,81E+04
#chr2a	utg000001l_pilon_pilon	99.714	47578	67	50	1386732	1434291	3478621	3431095	0.0	87042

samtools faidx P1.renamed.total.fasta chr2a > filtered_chromosomes/chr2a.fasta
samtools faidx P1.renamed.total.fasta chr2b > filtered_chromosomes/chr2b.fasta
samtools faidx P1.renamed.total.fasta chr2c > filtered_chromosomes/chr2c.fasta


###1310304 - (2122923-2118484) = 1305865

samtools faidx P1.minimap.fna utg000001l_pilon_pilon:1297166-1305865 > filtered_chromosomes/P1.chr2b2.rc.fasta
samtools faidx P1.minimap.fna utg000001l_pilon_pilon:3420677-3431095 > filtered_chromosomes/P1.chr2a2.rc.fasta


########################################

samtools faidx P1.renamed.total.fasta chr3 > filtered_chromosomes/chr3.fasta

########################################

#chr4a   510470  13669944        60      61                                                                              
#chr4b   352068  14188929        352068  352069                                                                          
#chr4c   1510027 14541005        60      61                                                                              
#chr4d   1429807 16076207        60      61  

### at this chromsomes, chr4a and chr4b overlap, whcih is again due to the ribosomal cluster. again, we should blast the ribosomal cluster to the genome
#chr4a	utg000002l_pilon_pilon	97.274	181316	2112	2063	306003	486305	304983	484480	0.0	3,05E+08
#chr4b	utg000002l_pilon_pilon	94.733	51790	1181	1185	31008	82145	458333	509227	0.0	79096

#chr4b	utg000002l_pilon_pilon	99.114	31944	119	114	320223	352068	746948	778825	0.0	5,73E+04
#chr4c	utg000002l_pilon_pilon	99.772	51800	33	46	1458254	1510027	838719	786979	0.0	94924

#chr4c	utg000002l_pilon_pilon	98.640	35362	212	183	1	35277	2296366	2261189	0.0	62388
#chr4d	utg000002l_pilon_pilon	99.562	149252	277	258	1	149080	2299994	2449040	0.0	2,72E+08

makeblastdb -in P1.renamed.total.fasta -dbtype 'nucl' 
blastn -db P1.renamed.total.fasta -query AB305101.1.fasta -outfmt 6 > blast.rib.2.P1.renamed.txt

#AB305101.1      chr4a   100.000 752     0       0       1       752     471722  470971  0.0     1389                    
#AB305101.1      chr4a   99.601  752     2       1       1       752     463946  463196  0.0     1371                    
#AB305101.1      chr4a   99.468  752     0       4       1       752     487288  486541  0.0     1363                    
#AB305101.1      chr4a   99.202  752     0       6       1       752     479507  478762  0.0     1351                    
#AB305101.1      chr4a   98.404  752     2       8       1       752     495043  494302  0.0     1314                    
#AB305101.1      chr4a   97.358  757     6       13      1       752     502770  502023  0.0     1275                    
#AB305101.1      chr4b   98.803  752     3       6       1       752     42834   42089   0.0     1334                    
#AB305101.1      chr4b   98.672  753     3       7       1       752     19538   18792   0.0     1328                    
#AB305101.1      chr4b   98.670  752     4       6       1       752     35068   34323   0.0     1328                    
#AB305101.1      chr4b   98.537  752     3       8       1       752     27301   26558   0.0     1321                    
#AB305101.1      chr4b   98.005  752     3       10      1       752     50582   49843   0.0     1295                    
#AB305101.1      chr4b   92.328  769     27      26      1       752     3891    3138    0.0     1064                    
#AB305101.1      chr4b   97.561  410     4       6       1       407     11759   11353   0.0     697                     
#AB305101.1      chr4b   98.701  385     3       1       368     752     11352   10970   0.0     682      #

### due to the ribosomal cluster, we cut a piece from chr4a, this is till the end of the firt 463946

samtools faidx P1.renamed.total.fasta chr4a:1-463946 > filtered_chromosomes/chr4a.fasta
### and b is then started from the end of the last
samtools faidx P1.renamed.total.fasta chr4b:50582-352068 > filtered_chromosomes/chr4b.fasta

samtools faidx P1.minimap.fna utg000002l_pilon_pilon:778825-786979 > filtered_chromosomes/chr4b2.fasta

### strangly 4c is in reverse complement
samtools faidx P1.renamed.total.fasta chr4c > filtered_chromosomes/P1.chr4c.rc.fasta

samtools faidx P1.minimap.fna utg000002l_pilon_pilon:2296366-2299994 > filtered_chromosomes/chr4c2.fasta

samtools faidx P1.renamed.total.fasta chr4d > filtered_chromosomes/chr4d.fasta

################################

#chr5a   1253769 17529852        1253769 1253770                                                                         
#chr5b   700958  18783629        60      61                                                                              
#chr5c   615440  19496277        615440  615441                                                                          
#chr5d   1343789 20111725        60      61
    
## the canu assembly of chr5 consists of 4 portions, that of minasm and minimap of 3, with overlapping breakpoints

#chr5a	utg000009l_pilon_pilon	99.403	7039	9	19	1246751	1253769	1243147	1250172	0.0	1,27E+04
#chr5b	utg000009l_pilon_pilon	97.928	47497	440	381	1	47309	1257146	1304286	0.0	8,18E+04

#chr5b	utg000009l_pilon_pilon	97.931	69123	651	583	632109	700958	1886789	1955405	0.0	1,19E+08
#chr5c	utg000009l_pilon_pilon	98.649	153837	882	867	6539	160029	1982629	2135615	0.0	2,72E+08

#chr5c	utg000009l_pilon_pilon	99.164	11959	39	46	349728	361666	2324573	2336490	0.0	2,15E+04

#chr5d	utg000016l_pilon_pilon	99.317	52297	129	152	573875	626031	52209	1	0.0	94387
#chr5d	utg000016l_pilon_pilon	99.382	13583	19	20	560184	573730	65877	52324	0.0	24557
#chr5d	utg000016l_pilon_pilon	98.813	208316	1065	1007	352248	560097	273310	65937	0.0	3,70E+08
#chr5d	utg000016l_pilon_pilon	93.375	19743	633	514	331969	351431	293486	274139	0.0	2,86E+04
#chr5d	utg000016l_pilon_pilon	98.747	61140	330	319	270698	331679	354627	293766	0.0	1,08E+08
#chr5d	utg000016l_pilon_pilon	99.916	19132	1	5	251441	270568	373855	354735	0.0	3,52E+04
#chr5d	utg000016l_pilon_pilon	98.667	149330	853	790	102042	250987	522838	374263	0.0	2,64E+08
#chr5d	utg000016l_pilon_pilon	96.528	102226	1677	1394	1	101669	624050	523140	0.0	1,67E+08
#chr5c	utg000016l_pilon_pilon	98.792	190847	1000	906	424963	615440	826426	636517	0.0	3,38E+08
#chr5c	utg000016l_pilon_pilon	99.826	6323	1	2	418267	424579	833084	826762	0.0	11607
#chr5c	utg000016l_pilon_pilon	99.948	44096	7	14	366249	410335	884262	840174	0.0	8,13E+04

#chr5d	utg000006l_pilon_pilon	99.680	195882	241	252	1148049	1343787	195654	15	0.0	3,58E+08
#chr5d	utg000006l_pilon_pilon	99.070	86517	336	322	1061551	1147905	281972	195763	0.0	1,55E+08
#chr5d	utg000006l_pilon_pilon	99.082	62513	226	222	998264	1060619	345106	282785	0.0	1,12E+08
#chr5d	utg000006l_pilon_pilon	94.595	13525	358	279	980048	993429	363015	349721	0.0	2,06E+04
#chr5d	utg000006l_pilon_pilon	97.643	85702	944	783	894319	979729	448191	363275	0.0	1,46E+08
#chr5d	utg000006l_pilon_pilon	99.739	11510	9	17	882534	894031	459935	448435	0.0	2,11E+04
#chr5d	utg000006l_pilon_pilon	98.395	106070	726	705	776621	882391	565457	460065	0.0	1,86E+08
#chr5d	utg000006l_pilon_pilon	98.016	138464	1253	1138	638504	776508	702959	565531	0.0	2,39E+08

samtools faidx P1.renamed.total.fasta chr5a > filtered_chromosomes/chr5a.fasta
samtools faidx P1.minimap.fna utg000009l_pilon_pilon:1250172-1257146 > filtered_chromosomes/chr5a2.fasta
samtools faidx P1.renamed.total.fasta chr5b > filtered_chromosomes/chr5b.fasta
### 1982629 - 6539 =  1976090
samtools faidx P1.minimap.fna utg000009l_pilon_pilon:1955405-1976090 > filtered_chromosomes/chr5b2.fasta
samtools faidx P1.renamed.total.fasta chr5c > filtered_chromosomes/chr5c.fasta
samtools faidx P1.minimap.fna utg000016l_pilon_pilon:624050-636517 > filtered_chromosomes/P1.chr5c2.rc.fasta
samtools faidx P1.renamed.total.fasta chr5d > filtered_chromosomes/chr5d.fasta

samtools faidx P1.renamed.total.fasta chr6 > filtered_chromosomes/chr6.fasta
samtools faidx P1.renamed.total.fasta chr7 > filtered_chromosomes/chr7.fasta
samtools faidx P1.renamed.total.fasta chr8 > filtered_chromosomes/chr8.fasta
samtools faidx P1.renamed.total.fasta mito > filtered_chromosomes/mito.fasta

############### some of the files in filtered_chromosomes are reverse complement so we reverse complement these

cd filtered_chromosomes

for chr in chr1a2 chr1b2 chr1c2 chr2a2 chr2b2 chr4c chr5c2
do
~/programs/seqtk/seqtk seq -r P1.$chr.rc.fasta > $chr.fasta
rm P1.$chr.rc.fasta
done

### first list all the files that are needed

cd filtered_chromosomes

sed 's/in.fasta/chr1a.fasta/g' fix_fasta.py > fix.fasta.chr1a.py
sed 's/out.fasta/fixed.chr1a.fasta/g' fix.fasta.chr1a.py > fix.fasta.chr1a.2.py 

sed 's/in.fasta/chr1a2.fasta/g' fix_fasta.py > fix.fasta.chr1a2.py
sed 's/out.fasta/fixed.chr1a2.fasta/g' fix.fasta.chr1a2.py > fix.fasta.chr1a2.2.py 

sed 's/in.fasta/chr1b.fasta/g' fix_fasta.py > fix.fasta.chr1b.py
sed 's/out.fasta/fixed.chr1b.fasta/g' fix.fasta.chr1b.py > fix.fasta.chr1b.2.py 

sed 's/in.fasta/chr1b2.fasta/g' fix_fasta.py > fix.fasta.chr1b2.py
sed 's/out.fasta/fixed.chr1b2.fasta/g' fix.fasta.chr1b2.py > fix.fasta.chr1b2.2.py 

sed 's/in.fasta/chr1c.fasta/g' fix_fasta.py > fix.fasta.chr1c.py
sed 's/out.fasta/fixed.chr1c.fasta/g' fix.fasta.chr1c.py > fix.fasta.chr1c.2.py 

sed 's/in.fasta/chr1c2.fasta/g' fix_fasta.py > fix.fasta.chr1c2.py
sed 's/out.fasta/fixed.chr1c2.fasta/g' fix.fasta.chr1c2.py > fix.fasta.chr1c2.2.py 


sed 's/in.fasta/chr1d.fasta/g' fix_fasta.py > fix.fasta.chr1d.py
sed 's/out.fasta/fixed.chr1d.fasta/g' fix.fasta.chr1d.py > fix.fasta.chr1d.2.py 

sed 's/in.fasta/chr2a.fasta/g' fix_fasta.py > fix.fasta.chr2a.py
sed 's/out.fasta/fixed.chr2a.fasta/g' fix.fasta.chr2a.py > fix.fasta.chr2a.2.py 

sed 's/in.fasta/chr2a2.fasta/g' fix_fasta.py > fix.fasta.chr2a2.py
sed 's/out.fasta/fixed.chr2a2.fasta/g' fix.fasta.chr2a2.py > fix.fasta.chr2a2.2.py 

sed 's/in.fasta/chr2b.fasta/g' fix_fasta.py > fix.fasta.chr2b.py
sed 's/out.fasta/fixed.chr2b.fasta/g' fix.fasta.chr2b.py > fix.fasta.chr2b.2.py 

sed 's/in.fasta/chr2b2.fasta/g' fix_fasta.py > fix.fasta.chr2b2.py
sed 's/out.fasta/fixed.chr2b2.fasta/g' fix.fasta.chr2b2.py > fix.fasta.chr2b2.2.py 

sed 's/in.fasta/chr2c.fasta/g' fix_fasta.py > fix.fasta.chr2c.py
sed 's/out.fasta/fixed.chr2c.fasta/g' fix.fasta.chr2c.py > fix.fasta.chr2c.2.py 

sed 's/in.fasta/chr3.fasta/g' fix_fasta.py > fix.fasta.chr3.py
sed 's/out.fasta/fixed.chr3.fasta/g' fix.fasta.chr3.py > fix.fasta.chr3.2.py 


sed 's/in.fasta/chr4a.fasta/g' fix_fasta.py > fix.fasta.chr4a.py
sed 's/out.fasta/fixed.chr4a.fasta/g' fix.fasta.chr4a.py > fix.fasta.chr4a.2.py 

sed 's/in.fasta/chr4b.fasta/g' fix_fasta.py > fix.fasta.chr4b.py
sed 's/out.fasta/fixed.chr4b.fasta/g' fix.fasta.chr4b.py > fix.fasta.chr4b.2.py 

sed 's/in.fasta/chr4b2.fasta/g' fix_fasta.py > fix.fasta.chr4b2.py
sed 's/out.fasta/fixed.chr4b2.fasta/g' fix.fasta.chr4b2.py > fix.fasta.chr4b2.2.py 

sed 's/in.fasta/chr4c.fasta/g' fix_fasta.py > fix.fasta.chr4c.py
sed 's/out.fasta/fixed.chr4c.fasta/g' fix.fasta.chr4c.py > fix.fasta.chr4c.2.py 

sed 's/in.fasta/chr4c2.fasta/g' fix_fasta.py > fix.fasta.chr4c2.py
sed 's/out.fasta/fixed.chr4c2.fasta/g' fix.fasta.chr4c2.py > fix.fasta.chr4c2.2.py 

sed 's/in.fasta/chr4d.fasta/g' fix_fasta.py > fix.fasta.chr4d.py
sed 's/out.fasta/fixed.chr4d.fasta/g' fix.fasta.chr4d.py > fix.fasta.chr4d.2.py 

sed 's/in.fasta/chr5a.fasta/g' fix_fasta.py > fix.fasta.chr5a.py
sed 's/out.fasta/fixed.chr5a.fasta/g' fix.fasta.chr5a.py > fix.fasta.chr5a.2.py 

sed 's/in.fasta/chr5a2.fasta/g' fix_fasta.py > fix.fasta.chr5a2.py
sed 's/out.fasta/fixed.chr5a2.fasta/g' fix.fasta.chr5a2.py > fix.fasta.chr5a2.2.py 

sed 's/in.fasta/chr5b.fasta/g' fix_fasta.py > fix.fasta.chr5b.py
sed 's/out.fasta/fixed.chr5b.fasta/g' fix.fasta.chr5b.py > fix.fasta.chr5b.2.py 


sed 's/in.fasta/chr5b2.fasta/g' fix_fasta.py > fix.fasta.chr5b2.py
sed 's/out.fasta/fixed.chr5b2.fasta/g' fix.fasta.chr5b2.py > fix.fasta.chr5b2.2.py 

sed 's/in.fasta/chr5c.fasta/g' fix_fasta.py > fix.fasta.chr5c.py
sed 's/out.fasta/fixed.chr5c.fasta/g' fix.fasta.chr5c.py > fix.fasta.chr5c.2.py 

sed 's/in.fasta/chr5c2.fasta/g' fix_fasta.py > fix.fasta.chr5c2.py
sed 's/out.fasta/fixed.chr5c2.fasta/g' fix.fasta.chr5c2.py > fix.fasta.chr5c2.2.py 


sed 's/in.fasta/chr5d.fasta/g' fix_fasta.py > fix.fasta.chr5d.py
sed 's/out.fasta/fixed.chr5d.fasta/g' fix.fasta.chr5d.py > fix.fasta.chr5d.2.py 

sed 's/in.fasta/chr6.fasta/g' fix_fasta.py > fix.fasta.chr6.py
sed 's/out.fasta/fixed.chr6.fasta/g' fix.fasta.chr6.py > fix.fasta.chr6.2.py 

sed 's/in.fasta/chr7.fasta/g' fix_fasta.py > fix.fasta.chr7.py
sed 's/out.fasta/fixed.chr7.fasta/g' fix.fasta.chr7.py > fix.fasta.chr7.2.py 

sed 's/in.fasta/chr8.fasta/g' fix_fasta.py > fix.fasta.chr8.py
sed 's/out.fasta/fixed.chr8.fasta/g' fix.fasta.chr8.py > fix.fasta.chr8.2.py 

sed 's/in.fasta/mito.fasta/g' fix_fasta.py > fix.fasta.mito.py
sed 's/out.fasta/fixed.mito.fasta/g' fix.fasta.mito.py > fix.fasta.mito.2.py 


for pop in chr1a chr1a2 chr1b chr1b2 chr1c chr1c2 chr1d chr2a chr2a2 chr2b chr2b2 chr2c chr3 chr4a chr4b chr4b2 chr4c chr4c2 chr4d chr5a chr5a2 chr5b chr5b2 chr5c chr5c2 chr5d chr6 chr7 chr8 mito
do
python fix.fasta.$pop.2.py
rm fix.fasta.$pop.py
rm fix.fasta.$pop.2.py
done



##############################################


#### now we need to give each file of a chromosome the same name and then concat the sequence

cp fixed.chr3.fasta ../final_fastas/final.chr3.fasta
cp fixed.chr6.fasta ../final_fastas/final.chr6.fasta
cp fixed.chr7.fasta ../final_fastas/final.chr7.fasta
cp fixed.chr8.fasta ../final_fastas/final.chr8.fasta
cp fixed.mito.fasta ../final_fastas/final.mito.fasta



for name in chr1a chr1a2 chr1b chr1b2 chr1c chr1c2 chr1d 
do
echo ">chr1" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr1a.fasta sub.chr1a2.fasta sub.chr1b.fasta sub.chr1b2.fasta sub.chr1c.fasta sub.chr1c2.fasta sub.chr1d.fasta > ../final_fastas/final.chr1.fasta
rm sub.*


for name in chr2a chr2a2 chr2b chr2b2 chr2c
do
echo ">chr2" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr2a.fasta sub.chr2a2.fasta sub.chr2b.fasta sub.chr2b2.fasta sub.chr2c.fasta > ../final_fastas/final.chr2.fasta

for name in chr4a chr4b chr4b2 chr4c chr4c2 chr4d 
do
echo ">chr4" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr4a.fasta sub.chr4b.fasta sub.chr4b2.fasta sub.chr4c.fasta sub.chr4c2.fasta sub.chr4d.fasta > ../final_fastas/final.chr4.fasta

for name in chr5a chr5a2 chr5b chr5b2 chr5c chr5c2 chr5d
do
echo ">chr5" > sub.$name.fasta
awk 'NR>1 {print $0}' fixed.$name.fasta >> sub.$name.fasta
done
seqkit concat sub.chr5a.fasta sub.chr5a2.fasta sub.chr5b.fasta sub.chr5b2.fasta sub.chr5c.fasta sub.chr5c2.fasta sub.chr5d.fasta > ../final_fastas/final.chr5.fasta

cd ../final_fastas

cat final.chr1.fasta final.chr2.fasta final.chr3.fasta final.chr4.fasta final.chr5.fasta final.chr6.fasta final.chr7.fasta final.chr8.fasta final.mito.fasta > P1.final.total.fasta

######################################################################################################

## final P1 sizes
#chr1    4667511 6       60      61                                                 
#chr2    4880141 4745315 60      61                                                 
#chr3    4005027 9706798 60      61                                                 
#chr4    3717051 13778582        60      61                                         
#chr5    3954085 17557590        60      61                                         
#chr6    3794862 21577583        60      61                                         
#chr7    1755998 25435699        60      61                                         
#chr8    1743926 27220970        60      61                                         
#mito    30993   28993968        60      61       

## final P0 size

#chr1    4663364 6       60      61                                                 
#chr2    4868082 4741099 60      61                                                 
#chr3    4026132 9690322 60      61                                                 
#chr4    3767989 13783563        60      61                                         
#chr5    3924104 17614358        60      61                                         
#chr6    3857747 21603870        60      
#chr7    1724474 25525919        60      61                                         
#chr8    1792139 27279141        60      61                                         
#mitochondrion   31003   29101164        60      61  
                             

########################
#######################
###################### final cleaning


cd 21canu

mkdir final_cleaning

/programs/bwa-0.7.15/bwa index P1.final.total.fasta
/programs/bwa-0.7.15/bwa mem -t 25 -x ont2d P1.final.total.fasta Asp.trim.assem/Asp.trim.assem.trimmedReads.fasta.gz > final_cleaning/mapping.sam

cd final_cleaning

racon -m 8 -x -6 -g -8 -w 500 -t 25 ../Asp.trim.assem/Asp.trim.assem.trimmedReads.fasta.gz mapping.sam ../P1.final.total.fasta > racon.fasta

/programs/bwa-0.7.15/bwa index racon.fasta
/programs/bwa-0.7.15/bwa mem -t 16 racon.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > P1.bam
samtools sort P1.bam -@ 28 -o P1.sort.bam
samtools index P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome racon.fasta --frags P1.sort.bam --outdir pilon1 --output pilon1 --changes


/programs/bwa-0.7.15/bwa index pilon1/pilon1.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon1/pilon1.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon1/P1.bam
samtools sort pilon1/P1.bam -@ 28 -o pilon1/P1.sort.bam
samtools index pilon1/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon1/pilon1.fasta --frags pilon1/P1.sort.bam --outdir pilon2 --output pilon2 --changes

/programs/bwa-0.7.15/bwa index pilon2/pilon2.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon2/pilon2.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon2/P1.bam
samtools sort pilon2/P1.bam -@ 28 -o pilon2/P1.sort.bam
samtools index pilon2/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon2/pilon2.fasta --frags pilon2/P1.sort.bam --outdir pilon3 --output pilon3 --changes

/programs/bwa-0.7.15/bwa index pilon3/pilon3.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon3/pilon3.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon3/P1.bam
samtools sort pilon3/P1.bam -@ 28 -o pilon3/P1.sort.bam
samtools index pilon3/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon3/pilon3.fasta --frags pilon3/P1.sort.bam --outdir pilon4 --output pilon4 --changes


/programs/bwa-0.7.15/bwa index pilon4/pilon4.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon4/pilon4.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon4/P1.bam
samtools sort pilon4/P1.bam -@ 28 -o pilon4/P1.sort.bam
samtools index pilon4/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon4/pilon4.fasta --frags pilon4/P1.sort.bam --outdir pilon5 --output pilon5 --changes

/programs/bwa-0.7.15/bwa index pilon5/pilon5.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon5/pilon5.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon5/P1.bam
samtools sort pilon5/P1.bam -@ 28 -o pilon5/P1.sort.bam
samtools index pilon5/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon5/pilon5.fasta --frags pilon5/P1.sort.bam --outdir pilon6 --output pilon6 --changes

/programs/bwa-0.7.15/bwa index pilon6/pilon6.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon6/pilon6.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon6/P1.bam
samtools sort pilon6/P1.bam -@ 28 -o pilon6/P1.sort.bam
samtools index pilon6/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon6/pilon6.fasta --frags pilon6/P1.sort.bam --outdir pilon7 --output pilon7 --changes

/programs/bwa-0.7.15/bwa index pilon7/pilon7.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon7/pilon7.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon7/P1.bam
samtools sort pilon7/P1.bam -@ 28 -o pilon7/P1.sort.bam
samtools index pilon7/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon7/pilon7.fasta --frags pilon7/P1.sort.bam --outdir pilon8 --output pilon8 --changes

/programs/bwa-0.7.15/bwa index pilon8/pilon8.fasta
/programs/bwa-0.7.15/bwa mem -t 16 pilon8/pilon8.fasta  ../P1.R1.fq.gz ../P1.R2.fq.gz   | samtools view -Sbq 20  > pilon8/P1.bam
samtools sort pilon8/P1.bam -@ 28 -o pilon8/P1.sort.bam
samtools index pilon8/P1.sort.bam
java -jar /programs/pilon-1.23.jar --genome pilon8/pilon8.fasta --frags pilon8/P1.sort.bam --outdir pilon9 --output pilon9 --changes



cd /programs/BUSCO-docker
cp /Asp.fumigatus/fasta_cleanup/21canu/final_cleaning/pilon5/pilon5.fasta pilon5.final_cleaning.P1.fasta

sudo docker run -it --rm -v $(pwd):/home/working -w /home/working chrishah/busco-docker run_BUSCO.py \
--in pilon5.final_cleaning.P1.fasta  --out pilon5.final_cleaning.P1.busco -l ./eurotiomycetes_odb9 --mode genome --species aspergillus_fumigatus 



# BUSCO version is: 3.1.0                                                                                                                           
# The lineage dataset is: eurotiomycetes_odb9 (Creation date: 2016-02-13, number of species: 25, number of BUSCOs: 4046)                            
# To reproduce this run: python /usr/bin/run_BUSCO.py -i pilon5.final_cleaning.P1.fasta -o pilon5.final_cleaning.P1.busco -l ./eurotiomycetes_odb9/ -m genome -c 1 -sp aspergillus_fumigatus                                                                                                          
#                                                                                                                                                   
# Summarized benchmarking in BUSCO notation for file pilon5.final_cleaning.P1.fasta                                                                
# BUSCO was run in mode: genome                                                                                                                                                                                                                                                                                 
C:99.3%[S:99.2%,D:0.1%],F:0.2%,M:0.5%,n:4046                                                                                                                                                                                                                                                            
4017    Complete BUSCOs (C)                                                                                                                         
4013    Complete and single-copy BUSCOs (S)                                                                                                         
4       Complete and duplicated BUSCOs (D)                                                                                                          
10      Fragmented BUSCOs (F)                                                                                                                       
19      Missing BUSCOs (M)                                                                                                                          
4046    Total BUSCO groups searched             



### pilon5 seems like the right genome## now we copy it to the 20canu directory and start making a masked version

cd /Asp.fumigatus/fasta_cleanup/21canu/

cp /Asp.fumigatus/fasta_cleanup/21canu/final_cleaning/pilon5/pilon5.fasta P1.final.pilon.fasta    

/programs/RepeatMasker/RepeatMasker -species "Aspergillus fumigatus" -pa 10 -gff -xsmall  -dir P1.final.pilon P1.final.pilon.fasta

#vim names.txt
#>chr1
#>chr2
#>chr3
#>chr4
#>chr5
#>chr6
#>chr7
#>chr8
#>mitochondrion

cd /Asp.fumigatus/fasta_cleanup/21canu/P1.final.pilon

perl /programs/rename_fasta_headers.pl names.txt P1.masked.fasta> P1.masked.2.fasta


~/programs/seqtk/seqtk seq -r  P0.1.rc.fasta > P0.1.fasta

#############################################3

samtools faidx 

cat P0.1.fixed.fasta P0.2.fasta P0.3.fasta P0.4.fasta P0.5.fasta P0.6.fasta P0.7.fasta P0.8.fasta mito.P0.fasta > P0.masked.fasta

