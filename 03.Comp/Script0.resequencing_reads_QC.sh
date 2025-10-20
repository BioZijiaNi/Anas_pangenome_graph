#!/bin/bash
#SBATCH --job-name=S0.ncbi_reads_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=72:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script0_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script0_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate reads_qc_env
#######################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
DirName="1.B.6.NCBI.WGS.reads.qc"

Threads=30
SampleList="/pdata1/home/nizijia/repository/SampleList_1.B.6.txt"
GroupName="SampleList_1.B.6"

CacheDir=$WorkDir/Cache/$DirName
mkdir -p $CacheDir
Dir1=$CacheDir/UnFiltered_reads_fastqc
mkdir -p $Dir1
Dir2=$CacheDir/Filtered_reads_qc
mkdir -p $Dir2
Dir3=$CacheDir/Multiqc
mkdir -p $Dir3
CRDir="/pdata1/home/nizijia/repository/1.B.6.ncbi_duck_WGS_clean"
mkdir -p $CRDir

##################################################################
##################### YOUR RUNNING SCRIPTS #######################
##################################################################
###### 01 NCBI download reads Quality
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo "##############"
echo "> $sample"

### Fastqc
SampleFastqcDir=$Dir1/$sample
mkdir -p $SampleFastqcDir

fastqc \
-f fastq \
-t $Threads \
-o $SampleFastqcDir \
$fq1 \
$fq2
done

### Multiqc
multiqc $Dir1/* \
--force \
--filename ${GroupName}_unfiltered_reads \
--outdir $Dir3

##### 02 Reads filter
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo "$sample"
date

#--dedup
#--trim_poly_x
#--dup_calc_accuracy 5
fastp \
--thread $Threads \
--detect_adapter_for_pe \
--qualified_quality_phred 20 \
--in1 $fq1 --in2 $fq2 \
--out1 $CRDir/${sample}_clean_R1.fq.gz --out2 $CRDir/${sample}_clean_R2.fq.gz
done


##### 3 Filtered reads quality control
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

echo "##############"
echo "> $sample"

### Fastqc
SampleFastqcDir=$Dir2/${sample}_clean
mkdir -p $SampleFastqcDir

fastqc \
-f fastq \
-t $Threads \
-o $SampleFastqcDir \
$CRDir/${sample}_clean_R1.fq.gz \
$CRDir/${sample}_clean_R2.fq.gz
done

### Multiqc
multiqc $Dir2/* \
--force \
--filename ${GroupName}_clean_reads \
--outdir $Dir3

### 4 Create clean sample list
#cd $CRDir
#SampList.sh -f ${WorkDir}/data/${GroupName}_clean.txt fq.gz R1 R2
