#!/bin/bash
#SBATCH --job-name=S3.3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.3_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.3_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate plink_env
##############################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script3.cache/SNPs_VS_SVs
mkdir -p $CacheDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir
Threads=4
SVVCF=$WorkDir/Cache/Script3.cache/Graph_SVs/merged_population_sample_population.vcf
SNPVCF=$WorkDir/Cache/Script3.cache/Graph_surject_small_variants/merge_population_snps.filtered.vcf.gz
Ref="/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies/Chr.BJ.fa"
NumChrom=$WorkDir/data/A.1.2.old.new.ChromName

###################### Running ##########################
#########################################################
##### 01 SNPs and SVs correlation
##########################################################
##### 1.1 1Mb bins
Newfai=$CacheDir/new.fai
paste $NumChrom ${Ref}.fai | awk 'BEGIN { OFS="\t" } {print $2,$4,$5,$6,$7}' > $Newfai
bedtools makewindows -g $Newfai -w 1000000 > $CacheDir/1mb_windows.bed

##### 1.2 Sorting VCF
### SVs
bcftools view -h $SVVCF > $CacheDir/header
bcftools view $SVVCF | bedtools sort -i - -faidx ${Newfai} |cat $CacheDir/header - > $CacheDir/SVs.sort.vcf
### SNPs
bcftools view -h $SNPVCF > $CacheDir/header
bcftools view $SNPVCF | bedtools sort -i - -faidx ${Newfai} |cat $CacheDir/header - > $CacheDir/SNPs.sort.vcf
##### SVs and SNPs
bedtools coverage -a $CacheDir/SVs.sort.vcf -b $CacheDir/SNPs.sort.vcf | awk '{print $1,$2,$3,$(NF-3),$(NF-2),$(NF-1),$NF}' > $CacheDir/SVs_SNPs_overlap.bed

##### 1.3 coverage （SNPs AND SVs）
bedtools coverage -a $CacheDir/1mb_windows.bed -b $CacheDir/SNPs.sort.vcf -counts > $CacheDir/snp_counts.bed
bedtools coverage -a $CacheDir/1mb_windows.bed -b $CacheDir/SVs.sort.vcf -counts > $CacheDir/sv_counts.bed
