#!/bin/bash
#SBATCH --job-name=S2.3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=70:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.3_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.3_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate base
### Software
# multiqc
### Ref
# https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/fly
#######################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S2.real_reads_mapping_rate
mkdir -p $CacheDir
GraphDir=$WorkDir/graph/MC_graph_250314/
mkdir -p $GraphDir

Threads=24
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.clean.txt"
GraphName="duck_MC_graph_16"
OutName="Linear_minimap2"

##################################################################
##################### RUNNING SCRIPTS ############################
##################################################################
##### Minimap2 multiqc
OutName="Linear_minimap2"
ResultDir=$CacheDir/minimap2_mapping
mkdir -p $ResultDir
multiqc \
$ResultDir --force --outdir $ResultDir --filename $OutName

##### VG giraffe multiqc
OutName="d2_Graph_giraffe"
ResultDir=$CacheDir/vg_giraffe_mapping
mkdir -p $ResultDir
multiqc $ResultDir/*d2* --force --outdir $ResultDir --filename $OutName
