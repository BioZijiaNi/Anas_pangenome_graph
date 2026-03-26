#!/bin/bash
######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate ENV_W02_W04_Pop
### Software
# fastp v1.0.1
# fastqc
# multiqc
# bwa (0.7.19)
# samtools (v1.22.1)
# sambamba (v1.0.1)
# freebayes (v1.3.8)
########################################################################################


#######################################################################################
############################## Setting params #########################################
WorkDir="/home1/nizijia/"
CacheDir=$WorkDir/cache/S4.cacje
mkdir -p $CacheDir
RawQCDir=$CacheDir/Raw_reads_fastqc_QC_report
mkdir -p $RawQCDir
CleanDir=$CacheDir/Clean_Reads
mkdir -p $CleanDir
FilteredQCDir=$CacheDir/Filtered_reads_fastqc_QC_report
mkdir -p $FilteredQCDir
AlignmentDir=$CacheDir/Alignments
mkdir -p $AlignmentDir
IndexDir=$CacheDir/Index
mkdir -p $IndexDir
VCFDir=$CacheDir/VCFDir
mkdir -p $VCFDir
FBVCFDir=$CacheDir/Freebayes_VCFDir
mkdir -p $FBVCFDir
PopulationDir=$CacheDir/PopulationDir
mkdir -p $PopulationDir

PopVCFDir=$CacheDir/PopVCFDir
mkdir -p $PopVCFDir

TmpDir=/home1/nizijia/tmp


#SampleList="/home3/nizijia/repository/SampList_1.B.5.txt"
SampleList="/home3/nizijia/repository/SampList_1.B.5_AND_W02_W04_Ma_SB.txt"
QCThreads=6
QCParallel=12

Threads=40
Parallel=1


GATKParallel=20
HCThreads=8
VCFDir=$CacheDir/VCFDir
mkdir -p $VCFDir




TargetList=$CacheDir/Script1.1.samplelist
cat $SampleList > $TargetList
Ref=/home2/nizijia/repository/1.B.9.Pekin_C18_ZJU1.0_W_genome/BJ_ZJU1.0_W.fa


##################################################################
##################### RUNNING SCRIPTS ############################
##################################################################
#####################################
##### 01 Reads QC
#####################################
##### 1.1 Raw reads QC reports
### Parallel running
cat $TargetList | xargs -P $QCParallel -I{} bash -c '
line="{}"
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

/home/nizijia/download/bin/micromamba run -n base fastqc \
    -f fastq \
    -t '"$QCThreads"' \
    -o '"$RawQCDir"' \
    $fq1 \
    $fq2
'
### Multi-sample qc report
multiqc $RawQCDir/ -n Raw_reads_QC_report -o $CacheDir -f

##### 1.2 Filtering reads
cat $TargetList |xargs -P $QCParallel -I{} bash -c '
line="{}"
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

fastp \
  --thread '"$QCThreads"' \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --in1 $fq1 --in2 $fq2 \
  --out1 '"$CleanDir"'/${sample}_R1.fq.gz \
  --out2 '"$CleanDir"'/${sample}_R2.fq.gz \
  --html '"$CleanDir"'/${sample}.html \
  --json '"$CleanDir"'/${sample}.json
'

##### 1.3 Filtered reads QC reports
### Parallel running
cat $TargetList | xargs -P $QCParallel -I{} bash -c '
line="{}"
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

/home/nizijia/download/bin/micromamba run -n base fastqc \
    -f fastq \
    -t '"$QCThreads"' \
    -o '"$FilteredQCDir"' \
    '"$CleanDir"'/${sample}_R1.fq.gz \
    '"$CleanDir"'/${sample}_R2.fq.gz
'
### Multi-sample qc report
multiqc $FilteredQCDir/ -n filtered_reads_QC_report -o $CacheDir -f


###################################
##### 02 Alignment
###################################
##### 2.1 Create index
#bwa index $Ref -p $IndexDir/bwa_index

##### 2.2 - 2.3 Align, Filter, Sort, Markdup
cat $TargetList | xargs -P $Parallel -I{} bash -c '
line="{}"
arr=($line)
sample=${arr[0]}

##### 2.2 BWA mem alignment
### Align
bwa mem \
-R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
-t '"$Threads"' '"$IndexDir"'/bwa_index \
'"$CleanDir"'/${sample}_R1.fq.gz \
'"$CleanDir"'/${sample}_R2.fq.gz \
| samtools view -@ '"$Threads"' -bh > '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.bam

### Output alignments stats
# 步骤 2: 使用 samtools stats 统计上一步生成的 BAM 文件信息
# 统计结果会被重定向到后缀为 .stats 的文件
samtools stats \
-@ '"Threads"' \
'"$AlignmentDir"'/${sample}_ngs_mapped.bwa.bam > '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.bam.stats

##### 2.3 Filter Sort MarkDup
### 2.3 Remove unmapped reads and Sort (by samtools)
# -F 12, remove unmapped PE reads
samtools view -bh \
--threads '"$Threads"' \
-F 12 \
'"$AlignmentDir"'/${sample}_ngs_mapped.bwa.bam \
| samtools sort --threads '"$Threads"' --output-fmt bam  > '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.sort.bam
# Remove bam file
rm '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.bam

### Mark duplicate (by sambamba)
# --overflow-list-size，set temp file size
sambamba markdup \
-t '"$Threads"' \
--tmpdir '"$AlignmentDir"' \
--overflow-list-size 2000000 \
'"$AlignmentDir"'/${sample}_ngs_mapped.bwa.sort.bam '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.markdup.bam
### Create index
samtools index -@ '"$Threads"' '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.markdup.bam
# Remove sort.bam file
rm '"$AlignmentDir"'/${sample}_ngs_mapped.bwa.sort.bam'


##################################################################
##################### 03 PREPARATION (DICT) ######################
##################################################################
# GATK 需要参考基因组的 .dict 文件 (不仅仅是 .fai)
# 检查是否存在，不存在则生成
RefDict="${Ref%.*}.dict"
if [ ! -f "$RefDict" ]; then
    echo "##### Creating Sequence Dictionary for Reference #####"
    gatk CreateSequenceDictionary \
        -R $Ref \
        -O $RefDict
fi

##################################################################
##################### 04 GATK HaplotypeCaller ####################
##################################################################
# 直接使用 markdup.bam 进行变异检测
# 输出为 GVCF 格式，用于后续多样本联合分析

echo "##### Running GATK HaplotypeCaller (Generate GVCFs) #####"

cat $TargetList | xargs -P $GATKParallel -I{} bash -c '
line="{}"
arr=($line)
sample=${arr[0]}

# 输入：上一步生成的 markdup bam
INPUT_BAM="'"$AlignmentDir"'/${sample}_ngs_mapped.bwa.markdup.bam"
# 输出：GVCF
OUTPUT_GVCF="'"$VCFDir"'/${sample}.g.vcf.gz"

# 检查 BAM 是否存在，避免报错
if [ -f "$INPUT_BAM" ]; then
    gatk --java-options "-Xmx8g" HaplotypeCaller \
       -R '"$Ref"' \
       -I $INPUT_BAM \
       -O $OUTPUT_GVCF \
       -ERC GVCF \
       --native-pair-hmm-threads '"$HCThreads"'
else
    echo "Error: BAM file for $sample not found at $INPUT_BAM"
fi
'

##################################################################
##################### 05 JOINT GENOTYPING ########################
##################################################################
# 设置本阶段特有变量
DB_Path="$CacheDir/GenomicsDB"
RawVCF="$PopVCFDir/cohort_raw.vcf.gz"
FinalVCF="$PopVCFDir/cohort_final_PASS.vcf.gz"
# 建议指定一个包含主要染色体的 interval list，或者直接指定染色体名（如 chr1）
# 如果不指定 -L，GenomicsDBImport 会非常慢
INTERVALS="$CacheDir/BJ_ZJU1.0_W_genome_intervals.list" # 请确保你有这个文件，或者改为特定的染色体名


MEM="64G"                     # 给 Java 核心预留足够内存
GATK_THREADS=24                # 几十个样本，联合分型时线程可设为 16-32
BATCH_SIZE=85                  # 设为你的样本总数，让它一次性读入所有样本


##### 5.0
# 在脚本开头运行这一行，动态生成 list 文件
cut -f1 ${Ref}.fai > $INTERVALS

##### 5.1 生成 Sample Map 文件
SampleMap="$CacheDir/cohort_sample_map.txt"

rm -f $SampleMap
cat $TargetList | while read line; do
    arr=($line)
    sample=${arr[0]}
    gvcf="$VCFDir/${sample}.g.vcf.gz"
    echo -e "$sample\t$gvcf" >> $SampleMap
done

##### 5.2 GenomicsDBImport
# 几十个样本建议 batch-size 设为样本总数
echo "##### Importing GVCFs into GenomicsDB #####"
gatk --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport \
    --genomicsdb-workspace-path $DB_Path \
    --sample-name-map $SampleMap \
    -L $INTERVALS \
    --batch-size $BATCH_SIZE \
    --max-num-intervals-to-import-in-parallel 8 \
    --overwrite-existing-genomicsdb-workspace true

##### 5.3 GenotypeGVCFs
echo "##### Running GenotypeGVCFs #####"
gatk --java-options "-Xmx$MEM" GenotypeGVCFs \
    -R $Ref \
    -V gendb://$DB_Path \
    -O $RawVCF \
    --tmp-dir $TmpDir

##################################################################
##################### 06 VARIANT FILTERING #######################
##################################################################
echo "##### Filtering Variants (Hard Filtering) #####"

# 6.1 提取并过滤 SNPs
gatk SelectVariants -V $RawVCF -select-type SNP -O $VCFDir/cohort_snps.vcf.gz

gatk VariantFiltration \
    -V $VCFDir/cohort_snps.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_HARD_FILTER" \
    -O $VCFDir/cohort_snps_filtered.vcf.gz

# 6.2 提取并过滤 Indels
gatk SelectVariants -V $RawVCF -select-type INDEL -O $VCFDir/cohort_indels.vcf.gz

gatk VariantFiltration \
    -V $VCFDir/cohort_indels.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
    --filter-name "INDEL_HARD_FILTER" \
    -O $VCFDir/cohort_indels_filtered.vcf.gz

# 6.3 合并并提取 PASS 位点
gatk MergeVcfs \
    -I $VCFDir/cohort_snps_filtered.vcf.gz \
    -I $VCFDir/cohort_indels_filtered.vcf.gz \
    -O $VCFDir/cohort_combined_flagged.vcf.gz

gatk SelectVariants \
    -V $VCFDir/cohort_combined_flagged.vcf.gz \
    --exclude-filtered \
    -O $FinalVCF

echo "##### Pipeline Finished! Final VCF: $FinalVCF #####"

##################################################################
##################### 06 VARIANT FILTERING (Parallel) #######################
##################################################################
# 变量定义 (确保这些变量在上文已经定义好)
RawVCF="$PopVCFDir/cohort_raw.vcf.gz"
FinalVCF="$PopVCFDir/cohort_final_PASS.vcf.gz"
# 170+ 序列，并行数建议设为 CPU核心数的一半左右
FilteringParallel=4
ChromMem="8g"

echo "##### Starting Parallel Filtering for 170+ Scaffolds #####"
# 1. 获取序列列表
# 确保 Ref 变量指向正确，且 .fai 文件存在
Chroms=$(cut -f1 ${Ref}.fai)

# 2. 检查 VCF 索引
if [ ! -f "${RawVCF}.tbi" ]; then
    echo "Index not found, indexing $RawVCF ..."
    gatk IndexFeatureFile -I $RawVCF
fi

# 3. 定义函数
process_variant_chrom() {
    local chrom=$1
    local ref=$2
    local raw_vcf=$3
    local out_dir=$4
    local mem=$5

    local safe_name=$(echo $chrom | sed 's/[|/]/_/g')

    echo "[$(date +%T)] Processing $chrom ..."

    # 1. 提取该染色体的 SNP 并进行过滤 (VariantFiltration)
    # 第一步：提取并打上过滤标签
    gatk --java-options "-Xmx${mem}" VariantFiltration \
        -R $ref -V $raw_vcf -L $chrom \
        --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "SNP_HARD_FILTER" \
        -O $out_dir/tmp_${safe_name}_snps_flagged.vcf.gz

    # 第二步：只保留通过过滤 (PASS) 的 SNP 位点
    gatk --java-options "-Xmx${mem}" SelectVariants \
        -V $out_dir/tmp_${safe_name}_snps_flagged.vcf.gz \
        -select-type SNP \
        --exclude-filtered \
        -O $out_dir/tmp_${safe_name}_snps_final.vcf.gz

    # 2. 提取该染色体的 Indel 并进行过滤
    gatk --java-options "-Xmx${mem}" VariantFiltration \
        -R $ref -V $raw_vcf -L $chrom \
        --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" \
        --filter-name "INDEL_HARD_FILTER" \
        -O $out_dir/tmp_${safe_name}_indels_flagged.vcf.gz

    gatk --java-options "-Xmx${mem}" SelectVariants \
        -V $out_dir/tmp_${safe_name}_indels_flagged.vcf.gz \
        -select-type INDEL \
        --exclude-filtered \
        -O $out_dir/tmp_${safe_name}_indels_final.vcf.gz

    # 3. 合并该染色体的 SNP 和 Indel
    gatk MergeVcfs \
        -I $out_dir/tmp_${safe_name}_snps_final.vcf.gz \
        -I $out_dir/tmp_${safe_name}_indels_final.vcf.gz \
        -O $out_dir/tmp_${safe_name}_final.vcf.gz

    # 4. 清理所有中间产生的临时文件
    rm $out_dir/tmp_${safe_name}_snps_flagged.vcf.gz* \
       $out_dir/tmp_${safe_name}_snps_final.vcf.gz* \
       $out_dir/tmp_${safe_name}_indels_flagged.vcf.gz* \
       $out_dir/tmp_${safe_name}_indels_final.vcf.gz*
}

export -f process_variant_chrom

# 4. 核心调度：并行执行
# 注意：在传递路径变量时加上了单引号 ''，防止路径中有空格导致参数错位
echo "$Chroms" | xargs -n 1 -P $FilteringParallel -I {} \
    bash -c "process_variant_chrom '{}' '$Ref' '$RawVCF' '$PopVCFDir' '$ChromMem'"

# 5. 合并 170+ 个文件
echo "##### Merging 170+ Scaffolds into final VCF #####"
# 检查是否真的生成了文件，如果上一步全部失败，这里应该停止
count=$(find $PopVCFDir -name "tmp_*_final.vcf.gz" | wc -l)
if [ "$count" -eq "0" ]; then
    echo "ERROR: No intermediate VCFs found! Check the logs above."
    exit 1
fi

find $PopVCFDir -name "tmp_*_final.vcf.gz" | sort > $PopVCFDir/all_scaffolds.list

gatk --java-options "-Xmx16g" MergeVcfs \
    -I $PopVCFDir/all_scaffolds.list \
    -O $FinalVCF

# 6. 清理
rm $PopVCFDir/tmp_*_final.vcf.gz* $PopVCFDir/all_scaffolds.list

echo "##### Hard Filtering Finished! Final File: $FinalVCF #####"

##################################################################
##################### 07 POPULATION VCF QC #######################
##################################################################
echo "##### Starting Population-level Filtering #####"

# 定义输入（接上一步 GATK 的 FinalVCF）
InputVCF=$FinalVCF # 即 $PopVCFDir/cohort_final_PASS.vcf.gz

# 1. 创建染色体名到数字的映射表 (bcftools 专用格式)
# 格式为: OldName NewName
cut -f1 ${Ref}.fai | awk '{print $1, NR}' > $CacheDir/Chromname_Numeric.map

# 2. 增加：提取双等位 SNP 过滤
# -v snps: 只要 SNP
# -m2 -M2: 限制等位基因数量最小为 2，最大为 2 (即严格双等位)
bcftools annotate --rename-chrs $CacheDir/Chromname_Numeric.map $InputVCF --threads $Threads | \
bcftools view -v snps -m2 -M2 \
    -i "F_MISSING <= 0.20 && MAF >= 0.01" \
    --threads $Threads -Oz \
    -o $PopulationDir/merge_population_snps.filtered.vcf.gz

bcftools index --threads $Threads $PopulationDir/merge_population_snps.filtered.vcf.gz --force

# 3. 提取特定染色体 (1-39 为常染色体, 40 为 Z, 135 为 W)
# 注意：你需要确认你的 .fai 文件中第 40 行和 135 行确实是 Z 和 W
R_STRING=$( (seq 1 39; echo 40; echo 135) | paste -s -d ',' )

bcftools view -r $R_STRING -Oz \
-o $PopulationDir/merge_population_snps.filtered.rmscaffold.vcf.gz \
$PopulationDir/merge_population_snps.filtered.vcf.gz

bcftools index -f $PopulationDir/merge_population_snps.filtered.rmscaffold.vcf.gz --force

# 4. 仅保留 1-39 号常染色体用于 ADMIXTURE
R_STRING=$(seq 1 39 | paste -s -d ',')

bcftools view -r $R_STRING -Oz \
-o $PopulationDir/merge_population_snps.filtered.autosomes.vcf.gz \
$PopulationDir/merge_population_snps.filtered.vcf.gz

bcftools index -f $PopulationDir/merge_population_snps.filtered.autosomes.vcf.gz --force

##################################################################
##################### 08 PCA ANALYSIS ############################
##################################################################
echo "##### Converting to PLINK and LD Pruning #####"

# 创建性染色体映射文件 (将数字编号映射到 PLINK 识别的性染色体代码)
# PLINK 默认常染色体 1-39, 40 设为 Z(23), W(24)
#echo "40 23" > $CacheDir/chr_sex_map.txt
#echo "135 24" >> $CacheDir/chr_sex_map.txt

# 0. 使用 --mac 1 剔除单态位点 (和下面all_snp数量一样)
plink --vcf $PopulationDir/merge_population_snps.filtered.autosomes.vcf.gz \
    --make-bed \
    --set-missing-var-ids @:#[BJ] \
    --mac 1 \
    --chr-set 39 \
    --allow-extra-chr \
    --out $PopulationDir/allsnp_rm_nonsnp

# 1. VCF 转 BED 并重命名变异 ID
#plink --vcf $PopulationDir/merge_population_snps.filtered.rmscaffold.vcf.gz \
#--update-chr $CacheDir/chr_sex_map.txt \
plink --vcf $PopulationDir/merge_population_snps.filtered.autosomes.vcf.gz \
    --make-bed \
    --set-missing-var-ids @:#[BJ] \
    --chr-set 39 \
    --allow-extra-chr \
    --out $PopulationDir/allsnp

# 2. LD Pruning (独立性检验)
# 设置阈值 r2 < 0.2 以获得更高质量的独立位点
plink --bfile $PopulationDir/allsnp \
    --indep-pairwise 50 5 0.2 \
    --chr-set 39 \
    --allow-extra-chr \
    --out $PopulationDir/allsnp_pruned

# 3. 提取修剪后的位点
plink --bfile $PopulationDir/allsnp \
    --extract $PopulationDir/allsnp_pruned.prune.in \
    --make-bed \
    --chr-set 39 \
    --allow-extra-chr \
    --out $PopulationDir/allsnp_filter

# 4. 运行 PCA
plink --bfile $PopulationDir/allsnp_filter \
    --pca 10 \
    --chr-set 39 \
    --allow-extra-chr \
    --out $PopulationDir/pca_result

##################################################################
##################### 09 ADMIXTURE ANALYSIS ######################
##################################################################
echo "##### Running Admixture Analysis #####"
mkdir -p $CacheDir/admixture_folder
cd $CacheDir/admixture_folder

# Admixture 输入文件路径
BED_FILE="$PopulationDir/allsnp_filter.bed"

for K in $(seq 1 10)
do
    echo "--- Running K=$K ---"
    /home/nizijia/download/bin/micromamba run -n plink_env admixture \
        --cv $BED_FILE $K -j${Threads} | tee log${K}.out
done

# 整理 CV 错误率，用于确定最佳 K 值
grep "CV" log*.out | awk '{print $3,$4}' | sed 's/):/ /' | sort -k1,1n > cv_errors.txt

# 为 Q 矩阵添加样本 ID 名
cat $PopulationDir/allsnp_filter.fam | awk '{print $1}' > ind2pop.txt

for i in $(seq 1 10)
do
    if [ -f "allsnp_filter.${i}.Q" ]; then
        paste -d ' ' ind2pop.txt allsnp_filter.${i}.Q > allsnp_filter.${i}.Q.with_id
    fi
done

echo "##### Pipeline Analysis Finished! #####"

################################################################
##### 10. 100kb bin admixture
################################################################
# 使用的自定义脚本
# AdmixtureUpdatefam.sh
# AdmixtureUpdatefamSelect.sh
# AdmixtureWindow.sh

##### 1. 参数配置区 (Configuration)
# --- 目录配置 ---
BADir=$CacheDir/W04_BinAdmixtureDir
AMWDir=$CacheDir/W04_Admixturewindow
AMWRDir=$CacheDir/W04_Admixturewindow_result


# --- 阈值与数值参数 ---
WINDOW_SIZE=100000       # 窗口大小 (100kb)
Q_THRESHOLD=0.99         # 有监督学习的 Q 值筛选阈值
MIND_FILTER=0.1          # PLINK 样本缺失率过滤阈值
CHR_COUNT=39             # 物种染色体数量

# --- 样本筛选参数 ---
EXCLUDE_KEYWORD="W04"    # 需要进行分析的sample

# --- 输出文件 ---
FinalResult=$CacheDir/W04_bin_admixture.csv

# 创建必要目录
mkdir -p $BADir $AMWDir $AMWRDir

##### 2. 执行逻辑
# --- 2.1 创建 Genome Windows ---
# 提取染色体长度并生成 bed 窗口
cut -f1,2 ${Ref}.fai | awk -v OFS='\t' '{print NR, $2}' > $BADir/my_genome.chrom.sizes
bedtools makewindows -g $BADir/my_genome.chrom.sizes -w $WINDOW_SIZE > $BADir/100kb_windows.bed

# --- 2.2 准备有监督分析的 FAM 文件 ---
# 将 allsnp_filter.fam 的第 6 列根据 Q 矩阵更新，高于阈值的作为参考群
AdmixtureUpdatefam.sh \
    $CacheDir/admixture_folder/allsnp_filter.2.Q.Sample \
    $PopulationDir/allsnp_filter.fam \
    $BADir/allsnp_filter.supervised.fam \
    $Q_THRESHOLD

# --- 2.3 筛选待剔除样本 (恢复你的原始逻辑) ---
# 1. 找到所有表型为 -9 的样本
# 2. 从中过滤出包含关键字 (如 W04) 的行
# 3. 结果存入 samples_to_remove.txt
awk '$6 == -9 {print $1, $2}' $BADir/allsnp_filter.supervised.fam | \
    grep -v "$EXCLUDE_KEYWORD" > $BADir/samples_to_remove.txt

# 生成专门用于 Bin-Admixture 的基础文件
# 使用 --remove 剔除不需要的样本，保留参考样本和其余待测样本
plink \
    --bfile $PopulationDir/allsnp_filter \
    --remove $BADir/samples_to_remove.txt \
    --mind $MIND_FILTER \
    --make-bed \
    --chr-set $CHR_COUNT \
    --allow-extra-chr \
    --out $BADir/allsnp_for_binAM

# 更新 FAM 文件的表型信息，确保有监督标签正确注入
AdmixtureUpdatefamSelect.sh \
    $BADir/allsnp_filter.supervised.fam \
    $BADir/allsnp_for_binAM.fam \
    $BADir/allsnp_for_binAM.update.fam

# --- 2.4 执行窗口化 Admixture 扫描 ---
echo "##### Running Admixture over Windows #####"
AdmixtureWindow.sh \
    $BADir/allsnp_for_binAM \
    $BADir/allsnp_for_binAM.update.fam \
    $AMWDir \
    $BADir/100kb_windows.bed

# --- 2.5 汇总结果 ---
# 数据都在代码运行的目录下面
mv $(pwd)/window* $AMWRDir

AdmixtureSumm \
    $AMWRDir \
    $BADir/allsnp_for_binAM.update.fam \
    $FinalResult

echo "##### Bin-Admixture Analysis Finished! #####"
