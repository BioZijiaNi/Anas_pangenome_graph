#!/bin/bash
#SBATCH --job-name=Script1.1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=1000:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script1_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.1/venv-cactus-v2.9.1/bin/activate
# mamba activate MC_298
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.8/venv-cactus-v2.9.8/bin/activate
### Software
# cactus v2.9.1
# vg v1.60
# bcftools
# vcflib
# samtools (v1.11)
# sambamba (v1.0.1)
# freebayes (v1.3.8)
########################################################################################

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
GraphStatsDir=$CacheDir/GAMstats
mkdir -p $GraphStatsDir
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.clean.txt"
GraphName="duck_MC_graph_16.d2"
Ref="/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies/Chr.BJ.fa"
Threads=72

##################################################################
##################### RUNNING SCRIPTS ############################
##################################################################

##############################################
##### 01 Create Snarls
##############################################
vg snarls --threads $Threads ${GraphDir}/${GraphName}.gbz > ${GraphDir}/${GraphName}.snarls

#############################################
##### 02 vg girraffe mapping
#############################################
# gam to pack
# gam to bam
cat $SampleList | while read line;
do
##### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo "##############"
echo "> $sample"
##### 2.1 Align reads to graph
# https://github.com/vgteam/vg/issues/4368
echo "vg giraffe "
vg giraffe \
-Z ${GraphDir}/${GraphName}.gbz \
-m ${GraphDir}/${GraphName}.min \
-d ${GraphDir}/${GraphName}.dist \
-x ${GraphDir}/${GraphName}.xg \
-f $fq1 -f $fq2 \
-N $sample \
-o gam \
--parameter-preset fast \
--rescue-algorithm none \
-t $Threads --progress > ${CacheDir}/${GraphName}-${sample}.gam

##### 2.2 Pack gam file
echo "vg pack "
vg pack \
-x ${GraphDir}/${GraphName}.xg \
--gam ${CacheDir}/${GraphName}-${sample}.gam \
--min-mapq 5 \
-o ${CacheDir}/${GraphName}-${sample}.pack \
-t ${Threads}

##### 2.3 Extract gam to linear bam format
### Gam to bam
# --sample NAME  set this sample name for all reads
vg surject \
--threads $Threads \
--bam-output \
--sample $sample \
--xg-name ${GraphDir}/${GraphName}.gbz \
${CacheDir}/${GraphName}-${sample}.gam > ${CacheDir}/${GraphName}-${sample}.bam
### Sort and Build index
samtools view -h --threads $Threads ${CacheDir}/${GraphName}-${sample}.bam \
| sed 's/BJ#0#//g' \
|samtools sort --threads $Threads --output-fmt bam > ${CacheDir}/${GraphName}-${sample}.sort.bam
samtools index -@ $Threads ${CacheDir}/${GraphName}-${sample}.sort.bam
### Mark duplicate (by sambamba)
# --overflow-list-size，set temp file size
sambamba markdup \
-t $Threads \
--tmpdir $CacheDir \
--overflow-list-size 10000000 \
${CacheDir}/${GraphName}-${sample}.sort.bam ${CacheDir}/${GraphName}-${sample}.sort.markdup.bam
# Create index
samtools index -@ $Threads ${CacheDir}/${GraphName}-${sample}.sort.markdup.bam

##### 2.4 Remove Tmp File
rm ${CacheDir}/${GraphName}-${sample}.gam
rm ${CacheDir}/${GraphName}-${sample}.bam
rm ${CacheDir}/${GraphName}-${sample}.sort.bam*

####################################
##### 03 freebayes call small snps
####################################
##### 3.1 Freebayes call small variants
freebayes-parallel <(fasta_generate_regions.py ${Ref}.fai 100000) $Threads \
-f $Ref \
${CacheDir}/${GraphName}-${sample}.sort.markdup.bam \
--strict-vcf \
--genotype-qualities \
--skip-coverage 1000  > ${RSDir}/${GraphName}-${sample}.small.variants.vcf

##### 3.2 Compress and build index
bcftools view \
${RSDir}/${GraphName}-${sample}.small.variants.vcf \
--threads $Threads -Oz > ${RSDir}/${GraphName}-${sample}.small.variants.vcf.gz
bcftools index ${RSDir}/${GraphName}-${sample}.small.variants.vcf.gz --threads $Threads -f -t

##### 3.3 Remove Tmp File
rm ${CacheDir}/${GraphName}-${sample}.sort.markdup.bam*
rm ${RSDir}/${GraphName}-${sample}.small.variants.vcf
done

####################################
##### 04 vg call SVs
####################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

##### 4.1 Vg call
# https://github.com/vgteam/vg/issues/4463, mentioned vg pack error
#-r ${GraphDir}/${GraphName}.snarls \
#-x ${GraphDir}/${GraphName}.xg
TimeOutMSG 3h vg call \
${GraphDir}/${GraphName}.gbz \
-r ${GraphDir}/${GraphName}.snarls \
-k ${CacheDir}/${GraphName}-${sample}.pack \
-A -a --gbz --progress \
-t ${Threads} \
-s $sample > ${ResultDir}/${GraphName}-${sample}.vcf

##### 4.2 Compress and construct index
bgzip --force --threads ${Threads} ${ResultDir}/${GraphName}-${sample}.vcf
bcftools index ${ResultDir}/${GraphName}-${sample}.vcf.gz --threads ${Threads} -f -t
done

##### 4.3 Sample-level filter
### Filter
# Reheader, PASS and Qual>=30, Split & Attach variant's length
bcftools view \
-i 'FILTER="PASS" & QUAL>=30' \
${ResultDir}/${GraphName}-${sample}.vcf.gz \
--threads ${Threads} \
|sed 's/BJ#0#//g' \
| bcftools norm -m -any --threads $Thread | vcflength > ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30.vcf.gz
bcftools index ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30.vcf.gz --threads ${Threads} -f -t
# SVs
bcftools view \
${ResultDir}/${GraphName}-${sample}_split_PASS_Q30.vcf.gz \
-i 'length>=50 | length<=-50' -Oz --threads $Thread > ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30_SVs.vcf.gz
bcftools index ${ResultDir}/${GraphName}-${sample}_split_PASS_Q30_SVs.vcf.gz --threads ${Threads} -f -t
done
