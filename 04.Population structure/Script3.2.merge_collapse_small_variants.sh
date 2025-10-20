#!/bin/bash
#SBATCH --job-name=S3.2.merge_VCF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.2_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.2_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate plink_env
#######################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script3.cache/Graph_surject_small_variants
mkdir -p $CacheDir
ResultDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_vg_pipeline
mkdir -p $ResultDir
RSDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_linear_small_variants
mkdir -p $RSDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir
NumChrom=$WorkDir/data/A.1.2.old.new.ChromName
SampleList="/pdata1/home/nizijia/repository/1.B.6.ncbi_duck_WGS/script/Table1.duck_WGS_SRR_code.tsv"
Threads=24
VCFSuffix="hardfiltering.vcf.gz"

################## Running #####################
################################################
##### 01 Merge population-level SNP VCF
################################################
### Get SRR code
cat $SampleList | awk -F'\t' '{print $1}' > $CacheDir/target_sample.txt
### Population-level merge - bcftools merge
bcftools merge \
--threads $Threads -m snps `ls $RSDir/*_${VCFSuffix}|grep -f $CacheDir/target_sample.txt` \
| bgzip --threads $Threads --force > $CacheDir/merge_population_snps.vcf.gz
tabix $CacheDir/merge_population_snps.vcf.gz --force
### Remove Cache file
rm $CacheDir/target_sample.txt

#################################################
##### 02 Filter population-level vcf
#################################################
# plink
# Modify chromosome names to be numerical, and filter.
bcftools view $CacheDir/merge_population_snps.vcf.gz --threads $Threads \
| bcftools annotate --rename-chrs $NumChrom --threads $Threads \
| bcftools view -i "F_MISSING <= 0.20" --threads $Threads \
| bcftools view -i "MAF >= 0.01" --threads $Threads -Oz > $CacheDir/merge_population_snps.filtered.vcf.gz
bcftools index --threads $Threads $CacheDir/merge_population_snps.filtered.vcf.gz --force

##############################################
##### 03 Population sturcture
##############################################
### 2.1 vcf to bed & add ID
plink \
--vcf $CacheDir/merge_population_snps.filtered.vcf.gz \
--make-bed \
--set-missing-var-ids @:#[ChrBJ]\$1,\$2 \
--chr-set 41 --allow-extra-chr \
--out $CacheDir/allsnp

### 2.2 plink filter
# For structure analysis
# 独立性检验
plink \
--bfile $CacheDir/allsnp \
--indep-pairwise 50 5 0.5 \
--make-bed \
--chr-set 41 --allow-extra-chr \
--out $CacheDir/allsnp_pruned

plink \
--bfile $CacheDir/allsnp \
--extract $CacheDir/allsnp_pruned.prune.in \
--make-bed \
--chr-set 41 --allow-extra-chr \
--out $CacheDir/allsnp_filter

### 2.3 PCA
plink \
--bfile $CacheDir/allsnp_filter \
--chr-set 41 --allow-extra-chr \
--pca \
--freq \
--out $CacheDir/pca_filter_snps

### 2.4 Admixture
# using indepentent snps
mkdir -p $CacheDir/admixture_folder
# Run
for K in `seq 2 8`
do
admixture --cv $CacheDir/allsnp_filter.bed $K -j${Threads} |tee $CacheDir/admixture_folder/log${K}.out
done
mv $WorkDir/script/allsnp* $CacheDir/admixture_folder/
#cat $CacheDir/allsnp.fam |awk '{print $1}' > $CacheDir/admixture_folder/ind2pop.txt

for i in {2..8}
do
cat $CacheDir/allsnp_filter.fam |awk '{print $1}'|paste - $CacheDir/admixture_folder/allsnp_filter.${i}.Q -d ' ' > $CacheDir/admixture_folder/allsnp_filter.${i}.Q.Sample
done
