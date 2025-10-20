#!/bin/bash
#SBATCH --job-name=S2.1.Graph_call_variants
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.1_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.1_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.1/venv-cactus-v2.9.1/bin/activate
### Ref
# https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/fly
##############################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S2.real_reads_mapping_rate
mkdir -p $CacheDir
ResultDir=$CacheDir/vg_giraffe_mapping
mkdir -p $ResultDir
GraphDir=$WorkDir/graph/MC_graph_250314/
mkdir -p $GraphDir
GraphName="duck_MC_graph_16.d2"
Threads=40
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.clean.txt"
########################################################################################

##################### YOUR RUNNING SCRIPTS #######################
##################################################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo "####################"
echo "> $sample gam stats:"
##### 1. Extract 1M reads
seqkit sample -s 123 -n 1000000 -j $Threads -o $CacheDir/${sample}.r1.fq.gz $fq1
seqkit sample -s 123 -n 1000000 -j $Threads -o $CacheDir/${sample}.r2.fq.gz $fq2
seqtk sample -s 123 $fq1 1000000 |gzip > $CacheDir/${sample}.r1.fq.gz
seqtk sample -s 123 $fq2 1000000 |gzip > $CacheDir/${sample}.r2.fq.gz

##### 2. Align reads to graph
# https://github.com/vgteam/vg/issues/4368
# --parameter-preset fast
echo "vg giraffe "
TimeOutMSG 10m vg giraffe \
-Z ${GraphDir}/${GraphName}.gbz \
-m ${GraphDir}/${GraphName}.min \
-d ${GraphDir}/${GraphName}.dist \
-f $CacheDir/${sample}.r1.fq.gz \
-f $CacheDir/${sample}.r2.fq.gz \
-N $sample \
-o gam \
-t $Threads --progress \
--parameter-preset fast > ${ResultDir}/${GraphName}-${sample}.gam

### GAM stats
vg stats \
--threads ${Threads} \
--alignments ${ResultDir}/${GraphName}-${sample}.gam \
${GraphDir}/${GraphName}.gbz > $ResultDir/${GraphName}-${sample}.gam.stats

### Remove cache
rm ${ResultDir}/${GraphName}-${sample}.gam

done
