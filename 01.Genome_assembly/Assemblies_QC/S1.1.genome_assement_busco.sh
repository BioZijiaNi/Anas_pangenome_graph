#!/bin/bash
#SBATCH --job-name=S1.1.busco
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.1.genome_assessment_busco_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.1.genome_assessment_busco_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate genome_evaluate
### Softawre
# busco == v5.8.2
# augustus == 3.4.0 
(Augustus used in BUSCO is 3.5.0, which encountered errors, so it was reverted to 3.4.0)
#######################################################################################

#######################################################################################
##### Setting Params #####
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/A.1.1.cache
mkdir -p $CacheDir
ResultDir=$WorkDir/Result_genome_evaluation/
mkdir -p $ResultDir
BuscodbDir=$WorkDir/data/busco_db

Threads=30
SampleList=$WorkDir/data/genome_list
########################################################################################

################## Running script ##################
####################################################
##### 01 Genome Busco
####################################################
cat $SampleList | while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fa=${arr[1]}

### Run Busco
busco \
-i $fa \
-m genome \
--offline -l $BuscodbDir/aves_odb10 \
--augustus \
--cpu $Threads \
--force \
--out_path $ResultDir --out busco_blast_${sample} 

done

