#!/bin/bash
#SBATCH --job-name=Script1.3.fix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=70:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1.3_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1.3_error.txt

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
##### 01 small variants
####################################
cat $SampleList | while read line;
do
##### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

##### 1 hardfiltering & rename samplename & only keep snps
# QUAL >= 30
# DP <= 100
# GD(GQ/DP) > 2
# SAF/(SAF+SAR) ≥ 0.2
### Select snps
echo "unknown   $sample" > $CacheDir/samplename_old_new
bcftools reheader ${RSDir}/${GraphName}-${sample}.small.variants.vcf.gz --sample $CacheDir/samplename_old_new \
| bcftools norm -m -any --threads $Threads \
| bcftools view --threads $Threads \
-i 'QUAL >= 30 & FORMAT/DP < 100 & GQ/FORMAT/DP > 2 & SAF/(SAF+SAR) >= 0.2' \
--types snps \
-Oz -o $RSDir/${sample}_hardfiltering.vcf.gz
### Create index
bcftools index $RSDir/${sample}_hardfiltering.vcf.gz --threads $Threads -f -t
done
