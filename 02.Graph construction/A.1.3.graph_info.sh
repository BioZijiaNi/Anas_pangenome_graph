#!/bin/bash
#SBATCH --job-name=A.1.3
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --time=48:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.3_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.3_error.txt

################## Enviroments setting #######################################################
################## 运行前提前激活环境 ########################################################
### ENV
# mamba activate MC_graph
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.1/venv-cactus-v2.9.1/bin/activate
##############################################################################################

################## Setting input and output directory ###############
########## Params ##########
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
ResultDir=$WorkDir/graph/MC_graph_250314_graph_info
mkdir -p $ResultDir
CacheDir=$WorkDir/Cache/A1.2.graph.growth
mkdir -p $CacheDir
GraphName=duck_MC_graph_16
GraphDir=/pdata1/home/nizijia/project/01.duck_graph_pangenome/graph/MC_graph_250314/
Threads=40

################## Running script ###################
##### Total graph info
# output: A.1.3_out.txt
for Graph in ${GraphName}.sv.gfa.gz ${GraphName}.d1.gfa.gz ${GraphName}.d2.gfa.gz ${GraphName}.gfa.gz ${GraphName}.full.gfa.gz
do
gzip -dc $GraphDir/$Graph > $CacheDir/graph.gfa
echo "> ${Graph}"
vg stats --size -l $CacheDir/graph.gfa --threads $Threads
rm $CacheDir/graph.gfa
done
