#!/bin/bash
#SBATCH --job-name=Script1.2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=70:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1.2_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1.2_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.1/venv-cactus-v2.9.1/bin/activate
### Software
# cactus v2.9.1
# vg v1.60
# bcftools
# vcflib
### Ref
# https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/fly

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline
mkdir -p $CacheDir
ResultDir=$CacheDir/Result_vg_pipeline
mkdir -p $ResultDir
RSDir=$CacheDir/Result_linear_small_variants
mkdir -p $RSDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.clean.txt"
GraphName="duck_MC_graph_16.d2"
Ref="/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies/Chr.BJ.fa"
Threads=30

##################################################################
##################### YOUR RUNNING SCRIPTS #######################
##################################################################
####################################
##### 01 SVs hard filtering
####################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

##### Sample-level filter (only keep SVs)
### Filter
# Reheader, PASS AND Qual>=30 AND GQ>2, Split & Attach variant's length
bcftools view \
-i 'FILTER="PASS" & QUAL>=30 & GQ/(FORMAT/DP) > 2' \
${ResultDir}/${GraphName}-${sample}.vcf.gz \
--threads ${Threads} \
| sed 's/^##contig=<ID=BJ#0#/##contig=<ID=/g' \
| bcftools norm -m -any --threads $Threads \
| vcflength > ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30.vcf.gz
bcftools index ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30.vcf.gz --threads ${Threads} -f -t
### Select SVs
bcftools view \
${ResultDir}/${GraphName}-${sample}_split_PASS_Q30.vcf.gz \
-i 'length>=50 | length<=-50' -Oz --threads $Threads > ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30_SVs.vcf.gz
bcftools index ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30_SVs.vcf.gz --threads ${Threads} -f -t
done
