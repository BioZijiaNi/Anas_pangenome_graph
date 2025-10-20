#!/bin/bash
#SBATCH --job-name=Script2.2
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --time=70:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.2_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.2_err.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.1/venv-cactus-v2.9.1/bin/activate
########################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S2.real_reads_mapping_rate
mkdir -p $CacheDir
ResultDir=$CacheDir/minimap2_mapping
mkdir -p $ResultDir
GraphDir=$WorkDir/graph/MC_graph_250314/
mkdir -p $GraphDir
Threads=40
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.clean.txt"
GraphName="duck_MC_graph_16.d2"
OutName="Linear_minimap2"

##################################################################
##################### YOUR RUNNING SCRIPTS #######################
##################################################################
##### 1 Extract BJ sequence
# cannot multi-threads
vg paths -F -Q 'BJ#0#' -x ${GraphDir}/${GraphName}.gbz > ${GraphDir}/BJ.paths.fa

##### 2 Running minimap2
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo "##############"
echo "> $sample"

### NGS reads mapping linear genome
minimap2 \
-ax sr -t $Threads \
$GraphDir/BJ.paths.fa \
$CacheDir/${sample}.r1.fq.gz \
$CacheDir/${sample}.r2.fq.gz \
|samtools view -@ $Threads -bh > ${ResultDir}/${GraphName}-${sample}.minimap2.bam

### Stats
# SAM
samtools stats \
--threads $Threads \
${ResultDir}/${GraphName}-${sample}.minimap2.bam > $ResultDir/${GraphName}-${sample}.minimap2.bam.stats
# Remove SAM
rm ${ResultDir}/${GraphName}-${sample}.minimap2.bam
#rm $CacheDir/${sample}.r1.fq.gz
#rm $CacheDir/${sample}.r2.fq.gz
done
