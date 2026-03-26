#!/bin/bash

################## Enviroments setting #################
################## 运行前提前激活环境 ##################
### ENV
# mamba activate MC_graph_py39
### Software
########################################################

################## Setting input and output directory ###############
########## Params ##########
WorkDir="/home2/nizijia/projects/1.B.1.graph_genome_evaluation"
CacheDir=$WorkDir/Cache/S10.Truvari_repmask_anno

GraphName="duck_MC_graph_16"
Threads=30

VCF="/home1/nizijia/cache/S8.cache/duck_MC_graph_16.split.SVs.TruvariCollapsed.vcf"
outvcf=$CacheDir/graph_SV_collapsed_repmask_anno.vcf

#######################################
##### 注释流程
########################################
graph="${GraphName}"

repmask_parallel.py \
-i $VCF \
-o $outvcf \
-p '-qq -e NCBI -lib /home2/nizijia/repository/1.B.6.TElib/Final_TElib/Final_TElib_3.fa -lcambig -nocut -div 50 -no_id -s {fasta} -pa {threads} -frag 100000' \
--repmask_path /home/nizijia/script/repmask.py \
-m 50 -M 100000 --threshold 0.5 \
-J 6 \
-T 8 \
-l 5000 \
-w /home1/nizijia/tmp/

bgzip $outvcf --threads $Threads --force
tabix ${outvcf}.gz --threads $Threads --force


### Final output for article
bcftools view -s HL,PZ,LW,LC,ZW,SB,W04,CAU-Wild1_1,Ma,W02,CAU_Laying_1_0,SX,LCWDuck,ZJU1_0,CAU_Pekin_2_0 ${outvcf}.gz \
| truvari anno lcr \
| truvari anno svinfo \
| bcftools query -f '%CHROM\t%POS\t%length.ref\t%length.alt\t%SVTYPE\t%LCR\t%INFO/RM_REF_repeat\t%INFO/RM_REF_clsfam\t%RM_REF_pdiv\t%INFO/RM_ALT_repeat\t%INFO/RM_ALT_clsfam\t%RM_ALT_pdiv\t%FORMAT\n' > ${outvcf}.REF.ALT.anno.tsv
