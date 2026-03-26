#!/bin/bash
#SBATCH --job-name=S4.3.2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=48:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script4.3.2_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script4.3.2_error.txt

################## Enviroments setting #######################################################
################## 运行前提前激活环境 ########################################################
### ENV
# mamba activate rnaseq
##############################################################################################

################## Setting input and output directory ###############
########## Params ##########
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
GenomeDir=$WorkDir/genome
CacheDir=$WorkDir/Cache/S4.3
mkdir -p $CacheDir
UnfilterDir=$CacheDir/Unfiltered_reads_QC
mkdir -p $UnfilterDir
FilteredDir=$CacheDir/Filtered_reads_QC
mkdir -p $FilteredDir
CleanDir=/pdata1/home/nizijia/repository/1.B.5.ncbi_duck_rna_clean
mkdir -p $CleanDir

SampleList="/pdata1/home/nizijia/repository/SampleList_RNAseq_Pekin_Mallard_PRJNA645648.txt"
Threads=30
Index=$GenomeDir/Chr.BJ.hisat2.index/index


#######################################
#####  01 reads QC
#######################################
##### 1.1 Unfiltered reads QC report
### Each sample QC report
cat $SampleList | while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}
# run
fastqc \
-f fastq \
-t $Threads \
-o $UnfilterDir \
$fq1 \
$fq2
done
### Multi-sample qc report
multiqc . -n unfiltered_reads_QC_report -o $UnfilterDir -f

##### 1.2 Filter reads
### Each sample QC report
cat $SampleList | while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}
# run
fastp \
--thread $Threads \
--trim_front1 15 --trim_front2 15 \
--trim_tail1 5 --trim_tail2 5 \
--detect_adapter_for_pe \
--in1 $fq1 --in2 $fq2 \
--out1 $CleanDir/${sample}_R1.fq.gz --out2 $CleanDir/${sample}_R2.fq.gz \
--unpaired1 $CleanDir/${sample}_unpaired1.fq.gz --unpaired2 $CleanDir/${sample}_unpaired2.fq.gz
done

##### 1.3 QC reads
### Each sample QC report
cat $SampleList | while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}
# run
fastqc \
-f fastq \
-t $Threads \
-o $FilteredDir \
$CleanDir/${sample}_R1.fq.gz \
$CleanDir/${sample}_R2.fq.gz
done
### Multi-sample qc report
multiqc . -n filtered_reads_QC_report -o $FilteredDir -f

##### 1.4 Batch align RNAseq reads to Chr.BJ genome
cat $SampleList | while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo ">${sample} Start:"
date
### Aligning
hisat2 \
-p ${Threads} \
-x ${Index} \
-S ${CacheDir}/${sample}_clean.sam \
-1 $CleanDir/${sample}_R1.fq.gz \
-2 $CleanDir/${sample}_R2.fq.gz --summary-file ${CacheDir}/${sample}_clean.summary

### sam to bam
# -f 3 : PE, properly alignment
# -F 12: exclude unmapping reads
# -F 256: exclude secondary reads
samtools view \
-@ ${Threads} -h -bS \
-f 3 -F 268 \
${CacheDir}/${sample}_clean.sam > ${CacheDir}/${sample}_clean.bam
### Remove sam file
rm ${CacheDir}/${sample}_clean.sam
done
EOF

##############################
### Quant		  ####
##############################
featureCounts \
-p --countReadPairs -B \
-g gene_name \
-T ${Threads} \
-a $GenomeDir/Chr.BJ.gtf \
-o ${CacheDir}/RNAseq_countmatrix.tsv \
`ls ${CacheDir}/*.bam`

cat ${CacheDir}/RNAseq_countmatrix.tsv \
| grep -v "^#" \
| sed "s|${CacheDir}/||g" \
| sed "s/_clean.bam//g" > ${CacheDir}/RNAseq_countmatrix_fix.tsv

FeatureCTPM.v1 \
-i ${CacheDir}/RNAseq_countmatrix_fix.tsv \
-o ${CacheDir}/RNAseq_TPM.tsv
