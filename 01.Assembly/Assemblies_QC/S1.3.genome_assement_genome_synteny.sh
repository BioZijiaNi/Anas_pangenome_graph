#!/bin/bash
#SBATCH --job-name=S1.3.synteny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.3.genome_assessment_synteny_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.3.genome_assessment_synteny_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate genome_evaluate_synteny
### Softawre
# plotsr == 1.1.1
# syri == 1.7.0
# minimap2 == 2.28-r1209
# /pdata1/home/nizijia/download/NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl == 1.41
### Subscript
#
#######################################################################################

#######################################################################################
##### Setting Params #####
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S1.3.genome_synteny
mkdir -p $CacheDir
ResultDir=$WorkDir/Result_genome_evaluation/genome_synteny
mkdir -p $ResultDir

Threads=16
SampleList=$WorkDir/data/genome_list_plotsr
genome2=/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies_whole/W02.fa
genome3=/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies_whole/W04.fa
genome1=/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies_whole/Ma.fa
genome4=/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies_whole/SB.fa
sample2=W02
sample3=W04
sample1=Ma
sample4=SB

AFDir=$ResultDir/adjusted_fasta
mkdir -p $AFDir

######################## Running script ########################
################################################################
##### 01 Chromosome-level assembly stats (Chr01-Chr30 ChrZ ChrW)
################################################################
##### 1.1 Extract Chr1-Chr30 ChrZ ChrW AND Fix Chrom number
### W02
samtools faidx \
$genome2 \
Chr{1..30} ChrZ ChrW \
| seqkit replace -p 'Chr(\d)$' -r 'Chr0${1}' > $CacheDir/${sample2}_cache.fa
### W04
samtools faidx \
$genome3 \
W04_Chr{1..30} W04_ChrZ W04_ChrW \
| seqkit replace -p '^W04_' -r "" \
| seqkit replace -p 'Chr(\d)$' -r 'Chr0${1}' > $CacheDir/${sample3}_cache.fa
### Ma
samtools faidx $genome1 Chr{1..30} Chr32_Z Chr33_W \
| seqkit replace -p '^(Chr32_|Chr33_)' -r "Chr" \
| seqkit replace -p 'Chr(\d)$' -r 'Chr0${1}' > $CacheDir/${sample1}_chr.fa
### SB
samtools faidx $genome4 Chr{1..30} Chr32_Z Chr33_W \
| seqkit replace -p '^(Chr32_|Chr33_)' -r "Chr" \
| seqkit replace -p 'Chr(\d)$' -r 'Chr0${1}' > $CacheDir/${sample4}_cache.fa

##### 1.2 Reorder and reverse (Ma as reference)
AFDir=$ResultDir/adjusted_fasta
mkdir -p $AFDir
FnaReorder.sh --ref $CacheDir/${sample1}_chr.fa --query $CacheDir/${sample2}_cache.fa -o $AFDir/${sample2} -t $Threads
FnaReorder.sh --ref $CacheDir/${sample1}_chr.fa --query $CacheDir/${sample3}_cache.fa -o $AFDir/${sample3} -t $Threads
FnaReorder.sh --ref $CacheDir/${sample1}_chr.fa --query $CacheDir/${sample4}_cache.fa -o $AFDir/${sample4} -t $Threads
cat $CacheDir/${sample1}_chr.fa > $AFDir/${sample1}.adjust.fa

##### 1.3 Create list
echo "$sample1 $sample2 $sample3 $sample4" | \
awk -F' ' '{for (i=1; i<NF; i++) print $i "\t" $(i+1)}' > $AFDir/list

##### 1.4 Run syri
cat $AFDir/list | while read line;
do
### Vars
arr=($line)
s1=${arr[0]}
s2=${arr[1]}

### Align genomes
# minimap2 align
minimap2 -ax asm5 -t $Threads --eqx $AFDir/${s1}.adjust.fa $AFDir/${s2}.adjust.fa \
| samtools sort --threads $Threads -O BAM - > $CacheDir/${s1}_${s2}_chr.bam
# index
samtools index -@ $Threads $CacheDir/${s1}_${s2}_chr.bam

### Running syri for finding structural rearrangements between A and B
rm $CacheDir/*${s1}_${s2}syri*
syri \
-c $CacheDir/${s1}_${s2}_chr.bam \
-r $AFDir/${s1}.adjust.fa \
-q $AFDir/${s2}.adjust.fa \
-F B \
--dir $AFDir \
--prefix ${s1}_${s2}
done

##### 1.5 Running plotsr
# -R Merge adjacent homologous regions
# -s minimum size of a SR to be plotted (default: 10000)
# --chrord $WorkDir/data/chrom_order_Chr1_30 \
plotsr \
--sr $AFDir/${sample1}_${sample2}syri.out \
--sr $AFDir/${sample2}_${sample3}syri.out \
--sr $AFDir/${sample3}_${sample4}syri.out \
--genomes $SampleList \
-R \
-s 100 \
-S 0.7 \
-H 12 -W 9 \
-o $ResultDir/Ma_W02_W04_SB_plotsr.pdf

##### 1.6 Running NGenomeSyn
### pairwise
cat $AFDir/list | while read line;
do
### Vars
arr=($line)
s1=${arr[0]}
s2=${arr[1]}

perl /pdata1/home/nizijia/download/NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl \
-InGenomeA $AFDir/${s1}.adjust.fa \
-InGenomeB $AFDir/${s2}.adjust.fa \
-MinLenA 1000 \
-MinLenB 1000 \
-MinAlnLen 100000 \
-MappingBin minimap2 \
-BinDir minimap2 \
-OutPrefix $CacheDir/${s1}_vs_${s2} \
-NumThreads $Threads
done

##### 1.7 Running four genome synteny
### Fix Config File (*.conf) using paramters below：
# vim /pdata1/home/nizijia/project/01.duck_graph_pangenome/Result_genome_evaluation/genome_synteny/Ma_W02_W04_SB.conf
# Copy AND paste the paramters in The Config File
<<EOF
SetParaFor = global
GenomeInfoFile1=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/Ma_vs_W02.A.chr.len
GenomeInfoFile2=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/Ma_vs_W02.B.chr.len
GenomeInfoFile3=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W04_vs_SB.A.chr.len
GenomeInfoFile4=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W04_vs_SB.B.chr.len
LinkFileRef1VsRef2=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/Ma_vs_W02.link
LinkFileRef2VsRef3=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W02_vs_W04.link
LinkFileRef3VsRef4=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W04_vs_SB.link
SetParaFor = Genome1
NormalizedScale=1
GenomeName=Ma
SetParaFor = Genome2
GenomeName=W02
NormalizedScale=1
SetParaFor = Genome3
NormalizedScale=1
GenomeName=W04
SetParaFor = Genome4
GenomeName=SB
NormalizedScale=1
EOF

### Running (Copy from *.mapping.sh)
/pdata1/home/nizijia/download/NGenomeSyn-1.41/bin/NGenomeSyn \
-InConf /pdata1/home/nizijia/project/01.duck_graph_pangenome/Result_genome_evaluation/genome_synteny/Ma_W02_W04_SB.conf \
-OutPut /pdata1/home/nizijia/project/01.duck_graph_pangenome/Result_genome_evaluation/genome_synteny/Ma_W02_W04_SB.pdf
EOF

#######################################################################
##### 02 Chromosome-level assembly stats (small chromosome Chr31-)
#######################################################################
<<EOF
##### 2.0 Remove cache
rm $CacheDir/*
##### 2.1 Extract Chr1-Chr30 ChrZ ChrW AND Fix Chrom number
### W02
samtools faidx \
$genome2 \
Chr{31..42} --continue > $CacheDir/${sample2}_cache.fa
### W04
samtools faidx \
$genome3 \
W04_Chr{31..42} --continue \
| seqkit replace -p '^W04_' -r "" > $CacheDir/${sample3}_cache.fa
### Ma
samtools faidx $genome1 Chr{31..42} --continue \
| seqkit replace -p 'Chr(\d)$' -r 'Chr0${1}' > $CacheDir/${sample1}_chr.fa
### SB
samtools faidx $genome4 Chr{31..42} --continue \
| seqkit replace -p 'Chr(\d)$' -r 'Chr0${1}' > $CacheDir/${sample4}_cache.fa

##### 2.2 Reorder and reverse (以Ma为参考)
#AFDir=$ResultDir/adjusted_fasta
#mkdir -p $AFDir
FnaReorder.sh --ref $CacheDir/${sample1}_chr.fa --query $CacheDir/${sample2}_cache.fa -o $AFDir/${sample2} -t $Threads
FnaReorder.sh --ref $CacheDir/${sample1}_chr.fa --query $CacheDir/${sample3}_cache.fa -o $AFDir/${sample3} -t $Threads
FnaReorder.sh --ref $CacheDir/${sample1}_chr.fa --query $CacheDir/${sample4}_cache.fa -o $AFDir/${sample4} -t $Threads
cat $CacheDir/${sample1}_chr.fa > $AFDir/${sample1}.adjust.fa

##### 2.3 Create list
echo "$sample1 $sample2 $sample3 $sample4" | \
awk -F' ' '{for (i=1; i<NF; i++) print $i "\t" $(i+1)}' > $AFDir/list
EOF

### 2.4 Running NGenomeSyn
### pairwise
cat $AFDir/list | while read line;
do
### Vars
arr=($line)
s1=${arr[0]}
s2=${arr[1]}

perl /pdata1/home/nizijia/download/NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl \
-InGenomeA $AFDir/${s1}.adjust.fa \
-InGenomeB $AFDir/${s2}.adjust.fa \
-MinLenA 1000 \
-MinLenB 1000 \
-MinAlnLen 100000 \
-MappingBin minimap2 \
-BinDir minimap2 \
-OutPrefix $CacheDir/${s1}_vs_${s2} \
-NumThreads $Threads

done

<<EOF
### 2.5 Fix Config File (*.conf) using paramters below：
# vim /pdata1/home/nizijia/project/01.duck_graph_pangenome/Result_genome_evaluation/genome_synteny/Ma_W02_W04_SB.small.conf
# Copy AND paste the paramters in The Config File
SetParaFor = global
GenomeInfoFile1=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/Ma_vs_W02.A.chr.len
GenomeInfoFile2=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/Ma_vs_W02.B.chr.len
GenomeInfoFile3=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W04_vs_SB.A.chr.len
GenomeInfoFile4=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W04_vs_SB.B.chr.len
LinkFileRef1VsRef2=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/Ma_vs_W02.link
LinkFileRef2VsRef3=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W02_vs_W04.link
LinkFileRef3VsRef4=/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/S1.3.genome_synteny/W04_vs_SB.link
SetParaFor = Genome1
NormalizedScale=1
GenomeName=Ma
ChrNameShow=1
ChrNameShiftY=10
SetParaFor = Genome2
GenomeName=W02
NormalizedScale=1
SetParaFor = Genome3
NormalizedScale=1
GenomeName=W04
SetParaFor = Genome4
GenomeName=SB
NormalizedScale=1
ChrNameShow=1
ChrNameShiftY=45
ChrNameShiftX=1
EOF

### 2.6 Running (Copy from *.mapping.sh)
/pdata1/home/nizijia/download/NGenomeSyn-1.41/bin/NGenomeSyn \
-InConf /pdata1/home/nizijia/project/01.duck_graph_pangenome/Result_genome_evaluation/genome_synteny/Ma_W02_W04_SB.small.conf \
-OutPut /pdata1/home/nizijia/project/01.duck_graph_pangenome/Result_genome_evaluation/genome_synteny/Ma_W02_W04_SB.small.pdf
