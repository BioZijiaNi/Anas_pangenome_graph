#!/bin/bash
#SBATCH --job-name=S1.4.SRs_genome_heterozygosity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.4.genome_assessment_SRs_genome_heterozygosity_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.4.genome_assessment_SRs_genome_heterozygosity_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate genome_evaluate_SRs
### Softawre
# jellyfish 2.2.10
# GenomeScope 2.0
#######################################################################################

#######################################################################################
#################################### Setting Params ###################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S1.genome_evaluation
mkdir -p $CacheDir
ResultDir=$WorkDir/Result_genome_evaluation/genome_heterozygosity
mkdir -p $ResultDir

Threads=16
SampleList=/pdata1/home/nizijia/repository/1.B.1.W02_W04_sequcencing_data/SampleList_W02_W04.txt
#SampleList=/pdata1/home/nizijia/project/01.duck_graph_pangenome/data/SampleList_BJ_SX_Ma_SB_NGS_genome.txt
Kmer=19
GENOME_SIZE="1200M"

################# Running script ###################
####################################################
#####  01 reads QC
####################################################
##### 1.1 Unfiltered reads QC report
### Each sample QC report
UnfilterDir=$CacheDir/unfiltered_reads_dir
mkdir -p $UnfilterDir

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
multiqc $UnfilterDir -n unfiltered_reads_QC_report -o $UnfilterDir -f
EOF

##### 1.2 Filter reads, remove low and dup
FilterDir=$CacheDir/Filtered_reads
mkdir -p $FilterDir
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
--dedup \
--detect_adapter_for_pe \
--in1 $fq1 --in2 $fq2 \
--out1 $FilterDir/${sample}_R1.fq.gz --out2 $FilterDir/${sample}_R2.fq.gz \
--unpaired1 $FilterDir/${sample}_unpaired1.fq.gz --unpaired2 $FilterDir/${sample}_unpaired2.fq.gz
done

##### 1.3 QC reads
FilteredQCDir=$CacheDir/Filtered_reads_dir
mkdir -p $FilteredQCDir
echo "$FilteredQCDir"
cat $SampleList |head -n 1| while read line;
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
-o $FilteredQCDir \
$FilterDir/${sample}_R1.fq.gz \
$FilterDir/${sample}_R2.fq.gz
done
### Multi-sample qc report
multiqc $FilteredQCDir -n filtered_reads_QC_report -o $FilteredQCDir -f


###########################################
### 02 Run jellyfish + genomescope2
###########################################
cat $SampleList |while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fq1=${arr[1]}
fq2=${arr[2]}

##### Decompress fq file
gzip -dc $fq1 > $CacheDir/${sample}_R1.fq
gzip -dc $fq2 > $CacheDir/${sample}_R2.fq

###########################################
##### 2.1 jellyfish
###########################################
### jellyfish count
# --canonical: Count both strand, canonical representation
# --size=uint64: *Initial hash size
jellyfish count \
-m $Kmer \
--size $GENOME_SIZE \
--bf-size 10G \
--canonical \
-t $Threads \
$CacheDir/${sample}_R1.fq \
$CacheDir/${sample}_R2.fq \
-o $ResultDir/${sample}_${Kmer}mer.jf
### Remove
#rm $Cache/${sample}_R1.fq
#rm $Cache/${sample}_R2.fq

### jellyfish histo
# Uniq K-mer with a frequency of 1 occurrence was not counted in total K-mer
jellyfish histo \
-t $Threads $ResultDir/${sample}_${Kmer}mer.jf > $ResultDir/${sample}_${Kmer}mer.histo

###########################################
##### 2.2 GenomeScope
###########################################
mkdir -p $ResultDir/genomescope_${sample}_${Kmer}
# max_kmercov，设定高频kmer计数的截止值
genomescope2 \
--input $ResultDir/${sample}_${Kmer}mer.histo \
--ploidy 2 \
--kmer_length $Kmer \
--output $ResultDir/genomescope_${sample}_${Kmer}
done
