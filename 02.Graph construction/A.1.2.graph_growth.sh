#!/bin/bash
#SBATCH --job-name=A.1.2
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --time=48:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.2_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.2_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
### Software
# panacus v0.3.3
#######################################################################################

################## Setting input and output directory ###############
########## Params ##########
SampleList="/pdata1/home/nizijia/repository/SampleList_anas_level_assemblies.txt"
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
ResultDir=$WorkDir/graph/MC_graph_250314_graph_growth
mkdir -p $ResultDir
CacheDir=$WorkDir/Cache/A1.2.graph.growth
mkdir -p $CacheDir
GraphName=duck_MC_graph_16
Graph=/pdata1/home/nizijia/project/01.duck_graph_pangenome/graph/MC_graph_250314/duck_MC_graph_16.gfa.gz
Threads=40

################## Running script ###################
##### Preparing
### Assemblies order
cat $SampleList | awk '{print $1}'|tail -n +2 |head -n 16 > $CacheDir/${GraphName}.order.txt

###################################
##### 01 Ordered growth (Exculde BJ)
###################################
for i in bp
do
# Exclude Reference path
#zcat $Graph |grep -e '^W' | cut -f2-6 | awk '{ print $1 "#" $2 "#" $3 ":" $4 "-" $5 }' > $CacheDir/${GraphName}.paths.txt
grep -ie "BJ" $CacheDir/${GraphName}.paths.txt > $CacheDir/${GraphName}.paths.haplotype.txt

### Panacus
# -e: 排除BJ的nodes
panacus ordered-histgrowth \
--count $i -t ${Threads} -O $CacheDir/${GraphName}.order.txt \
--coverage 1,2,3,5,10,16 \
--quorum 0 \
-S -e $CacheDir/${GraphName}.paths.haplotype.txt \
$Graph > $ResultDir/${GraphName}.ordered-histgrowth.${i}.tsv

panacus-visualize \
$ResultDir/${GraphName}.ordered-histgrowth.${i}.tsv > $ResultDir/${GraphName}.ordered-histgrowth.${i}.pdf
done


########################################
##### 02 Graph growth curve (exclude BJ)
########################################
for i in bp
do

grep -ve "BJ" $CacheDir/${GraphName}.paths.txt > $CacheDir/${GraphName}.VBJ.paths.haplotype.txt
panacus histgrowth \
-t $Threads \
-l 1,2,1,1,1 \
-q 0,0,1,0.5,0.1 \
--count $i \
-S -a \
-s $CacheDir/${GraphName}.VBJ.paths.haplotype.txt $Graph > $ResultDir/${GraphName}.histgrowth.${i}.tsv

panacus-visualize \
$ResultDir/${GraphName}.histgrowth.${i}.tsv > $ResultDir/${GraphName}.histgrowth.${i}.pdf
done
EOF

##########################################
##### 03 Graph growth curve (entire graph)
##########################################
for i in bp
do
panacus histgrowth \
-t $Threads \
-l 1,1,2,2 \
-q 0,0.9,0,0.9 \
--count $i \
-S -a \
$Graph > $ResultDir/${GraphName}.whole.histgrowth.${i}.tsv

panacus-visualize \
--estimate_growth_params \
$ResultDir/${GraphName}.whole.histgrowth.${i}.tsv > $ResultDir/${GraphName}.whole.histgrowth.${i}.pdf
done
