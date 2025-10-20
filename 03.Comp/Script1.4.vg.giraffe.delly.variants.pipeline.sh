#!/bin/bash
#SBATCH --job-name=Script1.4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=1000:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1.4_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1.4_error.txt

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
# samtools (v1.11)
# sambamba (v1.0.1)
# freebayes (v1.3.8)
### Ref
# https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/fly

#######################################################################################
############################## Setting params #########################################
##### Params
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script1.4.vg.giraffe.convert.linear.SV
mkdir -p $CacheDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir
ResultDir=$WorkDir/Cache/Script1.4.vg.giraffe.convert.linear.SV/CandidateSVDir
mkdir -p $ResultDir
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.clean.txt"
GraphName="duck_MC_graph_16.d2"
LinearName="Chr_BJ"
Ref="/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies/Chr.BJ.fa"
Threads=40
MAFT="0.01"
FMISST="0.2"

##################################################################
##################### YOUR RUNNING SCRIPTS #######################
##################################################################
############################################################
##### 01 Bowtie align & 02 Manta call SVs
############################################################
cat $TargetList | while read line;
do
##### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

############################################################
##### 01 Bowtie align
############################################################
IndexDir=$CacheDir/Index/
mkdir -p $IndexDir

### Build index
/pdata1/home/nizijia/download/bin/micromamba run -n LinearCallSV_env \
bowtie2-build --threads $Threads $Ref $IndexDir/Chr_BJ

### Mapping
/pdata1/home/nizijia/download/bin/micromamba run -n LinearCallSV_env \
bowtie2 -x $IndexDir/Chr_BJ \
-1 $fq1 \
-2 $fq2 \
--no-mixed -X 2000 -p $Threads \
| samtools view -hb -@ $Threads -o ${CacheDir}/${GraphName}-${sample}.bam

### Sort and Build index
samtools view -h --threads $Threads ${CacheDir}/${GraphName}-${sample}.bam \
|samtools sort --threads $Threads --output-fmt bam > ${CacheDir}/${GraphName}-${sample}.sort.bam
samtools index -@ $Threads ${CacheDir}/${GraphName}-${sample}.sort.bam

### Mark duplicate (by sambamba)
# --overflow-list-size，set temp file size
sambamba markdup \
-t $Threads \
--tmpdir $CacheDir \
--overflow-list-size 10000000 \
${CacheDir}/${GraphName}-${sample}.sort.bam ${CacheDir}/${LinearName}.${sample}.sort.markdup.bam
# Create index
samtools index -@ $Threads ${CacheDir}/${LinearName}-${sample}.sort.markdup.bam

##################################################
##### 03 Manta call SVs
##################################################
rm -rf ${CacheDir}/${LinearName}.${sample}
### 3.1 Set config
/pdata1/home/nizijia/download/bin/micromamba run -n LinearCallSV_env \
configManta.py \
--bam ${CacheDir}/${LinearName}.${sample}.sort.markdup.bam \
--referenceFasta ${Ref} \
--runDir ${CacheDir}/${LinearName}.${sample}

### 3.2 run
/pdata1/home/nizijia/download/bin/micromamba run -n LinearCallSV_env \
python2 ${CacheDir}/${LinearName}.${sample}/runWorkflow.py --jobs $Threads

### 3.3 CandidateSV to ResultDir
cp ${CacheDir}/${LinearName}.${sample}/results/variants/diploidSV.vcf.gz $ResultDir/${LinearName}.${sample}.vcf.gz
cp ${CacheDir}/${LinearName}.${sample}/results/variants/diploidSV.vcf.gz.tbi $ResultDir/${LinearName}.${sample}.vcf.gz.tbi
done

####################################
##### 04 SVs hard filtering
####################################
cat $TargetList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

### Rename sample & Filtering
echo "SAMPLE1   $sample" > $CacheDir/samplename_old_new
bcftools reheader $ResultDir/${LinearName}.${sample}.vcf.gz --sample $CacheDir/samplename_old_new \
| bcftools norm -m -any \
| bcftools view -i 'FILTER=="PASS"' \
| bcftools view -i 'SVLEN >= 50 | SVLEN <= -50' \
| bcftools norm -f $Ref \
| bcftools view -Oz > $ResultDir/${LinearName}.${sample}.filter.vcf.gz
bcftools index $ResultDir/${LinearName}.${sample}.filter.vcf.gz --threads ${Threads} -f -t

### rm file
rm $ResultDir/${LinearName}.${sample}.vcf.gz*
done

#################################################
##### 5 Group-level SVs VCF
#################################################

#################################################
##### 5.1 Merge group-level SVs VCF
#################################################
SampleList="/pdata1/home/nizijia/repository/1.B.6.ncbi_duck_WGS/script/Table1.duck_WGS_SRR_code.tsv"
for group in `cat $SampleList | awk -F'\t' '{print $NF}' | uniq`
do
##### 1 merge to group-level VCF
echo "> $group start merge:"
date
# Get group info
cat $TargetList | grep "$group" | awk -F'\t' '{print $1}' > $CacheDir/target_sample.txt
# merge by SRR - bcftools merge
bcftools merge \
--threads $Threads -m none \
`ls $ResultDir/${LinearName}.*.filter.vcf.gz|grep -f $CacheDir/target_sample.txt` \
| bgzip --threads $Threads --force > $CacheDir/merge_${group}.vcf.gz
tabix $CacheDir/merge_${group}.vcf.gz --force
# remove redundant SVs - Truvari
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
done


###############################################################################################
##### 5.2 Group-level SVs Filter
###############################################################################################
##### Filter
# --maf 0.01: 过滤掉次要等位基因频率（MAF）低于 1% 的变异
# --geno 0.2: 过滤掉缺失率高于 20% 的变异位点
MAFT="0.01"
FMISST="0.5"

for group in `cat $TargetList | awk -F'\t' '{print $NF}' | uniq`
do

bcftools +fill-tags $CacheDir/ngs_graph_${group}.vcf.gz \
| bcftools view -i "F_MISSING < $FMISST" -Oz > $CacheDir/ngs_graph_${group}.fmissing.filtered.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.fmissing.filtered.vcf.gz --force

bcftools +fill-tags $CacheDir/ngs_graph_${group}.vcf.gz \
| bcftools view -e "MAF < $MAFT" -Oz > $CacheDir/ngs_graph_${group}.maf.filtered.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.maf.filtered.vcf.gz --force

bcftools +fill-tags $CacheDir/ngs_graph_${group}.vcf.gz \
| bcftools view -i "F_MISSING < $FMISST" \
| bcftools view -i "MAF >= $MAFT" -Oz > $CacheDir/ngs_graph_${group}.maf.fmissing.filtered.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.maf.fmissing.filtered.vcf.gz --force

done

##### Stats
for i in `ls $CacheDir/*vcf.gz`; do echo "> $i"; bcftools view $i|grep -v "^#"|wc -l;done > $CacheDir/population_SVs_number.txt



###########################################################
##### 06 Merge population-level VCF
###########################################################
##### 1 bcftools merge & Truvari merge
group="population"
cat $SampleList | awk -F'\t' '{print $1}' > $CacheDir/target_sample.txt
# merge by SRR - bcftools merge
bcftools merge \
--threads $Threads -m none \
`ls $ResultDir/${LinearName}.*.filter.vcf.gz | grep -f $CacheDir/target_sample.txt` \
| bgzip --threads $Threads --force > $CacheDir/merge_${group}.vcf.gz
tabix $CacheDir/merge_${group}.vcf.gz --force
# remove redundant SVs - Truvari
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
