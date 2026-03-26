#!/bin/bash
#SBATCH --job-name=S4.3.1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script4.3.1_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script4.3.1_error.txt

################## Enviroments setting #################
################## 运行前提前激活环境 ##################
### ENV
# mamba activate rnaseq
########################################################

################## Setting input and output directory ###############
########## Params ##########
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
GenomeDir=$WorkDir/genome
Genome=$GenomeDir/Chr.BJ.fa
Gff=$GenomeDir/Chr.BJ.gff
CacheDir=$WorkDir/Cache/S4.3
mkdir -p $CacheDir
IndexDir=$GenomeDir/Chr.BJ.hisat2.index
mkdir -p $IndexDir
Threads=1

######################################################
##### Run script
######################################################
##### Alignment
### Build Index
# GFF to GTF
gffread $Gff -T -o $GenomeDir/Chr.BJ.gtf
# Extract exon position info
hisat2_extract_exons.py \
$GenomeDir/Chr.BJ.gtf > $CacheDir/genomic.gtf.exon
# Extract splice sites position info
hisat2_extract_splice_sites.py \
$GenomeDir/Chr.BJ.gtf > $CacheDir/genomic.gtf.ss

# Construct index
hisat2-build \
--ss $CacheDir/genomic.gtf.ss \
--exon $CacheDir/genomic.gtf.exon \
-p $Threads \
$Genome \
$IndexDir/index
