#!/bin/bash
#SBATCH --job-name=S3.1.merge_VCF
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.1_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.1_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate plink_env
#######################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script3.cache/Graph_SVs
mkdir -p $CacheDir
ResultDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_vg_pipeline
mkdir -p $ResultDir
RSDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_linear_small_variants
mkdir -p $RSDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir
NumChrom=$WorkDir/data/A.1.2.old.new.ChromName
SampleList="/pdata1/home/nizijia/repository/1.B.6.ncbi_duck_WGS/script/Table1.duck_WGS_SRR_code.tsv"
InfoDir=$WorkDir/data/InfoDir
Threads=30

#################### Running script #######################
###########################################################
##### 01 merge population-level SVs for sample
###########################################################
##### 1 Population-level merge
group="population"
cat $SampleList | awk -F'\t' '{print $1}' > $CacheDir/target_sample.txt
# Population-level merge - bcftools merge
bcftools merge \
--threads $Threads -m none \
`ls $ResultDir/*_split_PASS_Q30_SVs.vcf.gz|grep -f $CacheDir/target_sample.txt` \
| bgzip --threads $Threads --force > $CacheDir/merge_${group}.vcf.gz
tabix $CacheDir/merge_${group}.vcf.gz --force
# Remove redundant SVs - Truvari
/pdata1/home/nizijia/download/bin/micromamba run -n MC_graph \
truvari collapse \
-i $CacheDir/merge_${group}.vcf.gz \
-o $CacheDir/truvari_merge_${group}.vcf \
-c $CacheDir/truvari_collapsed_${group}.vcf
# bgzip and tabix
bgzip $CacheDir/truvari_merge_${group}.vcf --threads $Threads --force
tabix $CacheDir/truvari_merge_${group}.vcf.gz --threads $Threads --force
# Remove Cache file
rm $CacheDir/target_sample.txt
#rm $CacheDir/merge_${group}.vcf*
rm $CacheDir/truvari_collapsed_${group}.vcf*

##### 2 Only SVs, Split, sort
bcftools norm -m- $CacheDir/truvari_merge_${group}.vcf.gz --threads $Threads \
| bcftools sort --max-mem 2G -Oz > $CacheDir/ngs_graph_${group}.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.vcf.gz --force
# Remove cache file
#rm $CacheDir/truvari_merge_${group}.vcf*

##### 3 Filter
bcftools +fill-tags $CacheDir/ngs_graph_${group}.vcf.gz \
| bcftools view -i "F_MISSING <= 0.8" \
| bcftools view -i "MAF >= 0.01" -Oz > $CacheDir/merged_population.vcf.gz
bcftools index --threads $Threads $CacheDir/merged_population.vcf.gz --force --tbi

### 4 Chrom to numeric
bcftools view $CacheDir/merged_population.vcf.gz \
| bcftools annotate --rename-chrs $NumChrom --threads $Threads> $CacheDir/merged_population_sample_population.vcf

#############################################################################
##### 04. PCA AND ADMIXTURE
#############################################################################
### 2.1 vcf to bed
plink \
--vcf $CacheDir/merged_population_sample_population.vcf \
--make-bed \
--chr-set 41 --allow-extra-chr \
--out $CacheDir/allsv

### 2.2 PCA
plink \
--bfile $CacheDir/allsv \
--chr-set 41 --allow-extra-chr \
--pca \
--freq \
--out $CacheDir/pca_sample_population

### 2.3 Admixture
# using indepentent snps
mkdir -p $CacheDir/admixture_folder
# 2.3.1 CV error
for K in `seq 2 8`
do
admixture --cv $CacheDir/allsv.bed $K -j${Threads} |tee $CacheDir/admixture_folder/log${K}.out
done
mv $WorkDir/script/allsv* $CacheDir/admixture_folder/
# 2.3.2 Result
for i in {2..8}
do
cat $CacheDir/allsv.fam |awk '{print $1}'|paste - $CacheDir/admixture_folder/allsv.${i}.Q -d ' ' > $CacheDir/admixture_folder/allsv.${i}.Q.Sample
done

#########################################################################################
##### 05 Fst between breeds (group-level SVs)
#########################################################################################
PopList=$WorkDir/data/InfoDir/B.1.Fst_Pop.txt

cat $PopList | while read line;
do
##### Vars
arr=($line)
Pop1=${arr[0]}
Pop2=${arr[1]}
group=${Pop1}_${Pop2}

##### 1 merge to group-level VCF
echo "> $group start merge:"
date
### Attach target sample list
cat $InfoDir/population_${Pop1}.txt $InfoDir/population_${Pop2}.txt > $CacheDir/target_sample.txt
### group-level merge - bcftools merge
bcftools merge \
--threads $Threads -m none \
`ls $ResultDir/*_split_PASS_Q30_SVs.vcf.gz|grep -f $CacheDir/target_sample.txt` \
| bgzip --threads $Threads --force > $CacheDir/merge_${group}.vcf.gz
tabix $CacheDir/merge_${group}.vcf.gz --force
# Remove redundant SVs - Truvari
/pdata1/home/nizijia/download/bin/micromamba run -n MC_graph \
truvari collapse \
-i $CacheDir/merge_${group}.vcf.gz \
-o $CacheDir/truvari_merge_${group}.vcf \
-c $CacheDir/truvari_collapsed_${group}.vcf
# bgzip and tabix
bgzip $CacheDir/truvari_merge_${group}.vcf --threads $Threads
tabix $CacheDir/truvari_merge_${group}.vcf.gz --threads $Threads
# Remove Cache file
rm $CacheDir/target_sample.txt
rm $CacheDir/merge_${group}.vcf*
rm $CacheDir/truvari_collapsed_${group}.vcf*

##### 2 Only SVs, Split, sort
bcftools norm -m- $CacheDir/truvari_merge_${group}.vcf.gz --threads $Threads \
| bcftools sort --max-mem 2G -Oz > $CacheDir/ngs_graph_${group}.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.vcf.gz --force
# Remove cache file
rm $CacheDir/truvari_merge_${group}.vcf*

##### 3 Filter
bcftools +fill-tags $CacheDir/ngs_graph_${group}.vcf.gz \
| bcftools view -i "F_MISSING < 0.2" \
| bcftools view -i "MAF >= 0.01" -Oz > $CacheDir/ngs_graph_${group}.maf.fmissing.filtered.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.maf.fmissing.filtered.vcf.gz --force
EOF

##### 4 SVs Fst caculate
VCF=$CacheDir/ngs_graph_${group}.maf.fmissing.filtered.vcf.gz
vcftools \
--gzvcf $VCF \
--weir-fst-pop $InfoDir/population_${Pop1}.txt \
--weir-fst-pop $InfoDir/population_${Pop2}.txt \
--out $CacheDir/${Pop1}_vs_${Pop2}_FST

done
