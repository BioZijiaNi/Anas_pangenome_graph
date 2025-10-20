#!/bin/bash
#SBATCH --job-name=A.1.0.pangenome_graph_construction
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.0.pangenome_graph_construction_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.0.pangenome_graph_construction_error.txt


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
### Ref
# https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/fly
#######################################################################################

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
SampleList="/pdata1/home/nizijia/repository/SampleList_anas_level_assemblies.txt"
RefName=BJ
GraphName=duck_MC_graph_16
Threads=72
SplitChromDir=chroms
ChromAlignDir=chrom-alignments
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
ResultDir=$WorkDir/graph/MC_graph_250314
mkdir -p $ResultDir
CacheDir=$WorkDir/Cache/A1.0.MC_graph_250314.cache
mkdir -p $CacheDir
TmpDir=/pdata1/home/nizijia/tmp
mkdir -p $TmpDir

##################### Running script #############################
### 01. Make SV graph
cactus-minigraph \
$CacheDir/jobstore \
$SampleList \
$CacheDir/${GraphName}.sv.gfa.gz \
--reference $RefName \
--workDir ${TmpDir}

### 02. Make the assembly-to-graph alignments
cactus-graphmap \
$CacheDir/jobstore \
$SampleList \
$CacheDir/${GraphName}.sv.gfa.gz \
$CacheDir/${GraphName}.paf \
--reference $RefName \
--outputFasta $CacheDir/${GraphName}.sv.gfa.fa.gz \
--mapCores $Threads \
--workDir ${TmpDir}

### 03. (Optional) Split the input assemblies and PAF into chromosomes
cactus-graphmap-split \
$CacheDir/jobstore \
$SampleList \
$CacheDir/${GraphName}.sv.gfa.gz \
$CacheDir/${GraphName}.paf \
--reference $RefName \
--outDir $CacheDir/$SplitChromDir \
--workDir ${TmpDir}

### 04. Create the Cactus base alignment and "raw" pangenome graph
cactus-align-batch \
$CacheDir/align_batch \
$CacheDir/${SplitChromDir}/chromfile.txt \
$CacheDir/$ChromAlignDir \
--alignOptions "--pangenome --reference $RefName --outVG" \
--workDir ${TmpDir}
EOF

### 05. Create and index the final pangenome graph and produce VCF files and vg-giraffe indexes
# --filter: nodes supported by no less than N assemblies
cactus-graphmap-join \
$CacheDir/graphmap_join \
--vg $CacheDir/${ChromAlignDir}/*.vg \
--hal $CacheDir/${ChromAlignDir}/*.hal \
--outDir $ResultDir \
--outName ${GraphName} \
--reference $RefName \
--vcf full clip filter \
--giraffe full clip filter \
--xg full clip filter \
--gfa full clip filter \
--unchopped-gfa full clip filter \
--filter 2 \
--workDir ${TmpDir}
