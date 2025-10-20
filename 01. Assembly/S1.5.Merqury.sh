##############################################
##### ENVs
# mamba activate genome_evaluation
##### Software
# merqury
#############################################

ResultDir=/home2/nizijia/projects/1.B.1.graph_genome_evaluation/Result1.3.genome_quality_control

########################################
##### 01. K select
########################################
best_k.sh 1225994370

########################################
##### 02. Merqury ONT (Drop)
########################################
kmer=21
TargetList=/home3/nizijia/repository/1.B.4.W02_W04_data/SampleList_W04_ONT.txt

cat $TargetList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
TGS=${arr[1]}
fa=${arr[2]}

### Running
meryl k=$kmer count output ONT_${sample}.mery $TGS
merqury.sh ONT_${sample}.mery $fa ONT_${sample}_result
done

########################################
##### 03. Merqury NGS
########################################
kmer=21
TargetList=/home3/nizijia/repository/1.B.4.W02_W04_data/SampleList_W02_W04_NGS.txt

cat $TargetList | tail -n +2 |while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}
fa=${arr[3]}
SampleDir=$ResultDir/NGS_${sample}
mkdir -p $SampleDir

### Build db
meryl k=${kmer} count output meryl1 $fq1
meryl k=${kmer} count output meryl2 $fq2
meryl union-sum output read.meryl meryl1 meryl2

merqury.sh read.meryl $fa NGS_${sample}_result

mv *NGS_${sample}_result* $SampleDir/
done
