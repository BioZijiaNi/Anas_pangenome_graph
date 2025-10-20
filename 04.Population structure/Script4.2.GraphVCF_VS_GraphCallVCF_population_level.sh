#!/bin/bash
#SBATCH --job-name=S4.2.comp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script4.2_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script4.2_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
### Ref
# https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/fly
##############################################################################################

################## Setting input and output directory ###############
########## Params ##########
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script4.cache/GraphVCF_vs_GraphCallVCF
mkdir -p $CacheDir
ResultDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_vg_pipeline
mkdir -p $ResultDir
RSDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_linear_small_variants
mkdir -p $RSDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir


SampleList="/pdata1/home/nizijia/repository/1.B.6.ncbi_duck_WGS/script/Table1.duck_WGS_SRR_code.tsv"

InputVcf="duck_MC_graph_16.d2.vcf.gz"
PopVCF="/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/Script3.cache/Graph_SVs/ngs_graph_population.vcf.gz"
FilteredPopVCF="/pdata1/home/nizijia/project/01.duck_graph_pangenome/Cache/Script3.cache/Graph_SVs/merged_population.vcf.gz"
group="population"
Threads=16
CompResDir=$CacheDir/${group}_comp

################## Running script ###################
<<EOF
####################################################
##### Graph VCF
####################################################
### Split VCF
# 会提取出: . 和 0, 需要进一步去除
# --samples-file $GraphDir/AnasPla_SB_list
bcftools view $GraphDir/$InputVcf --threads $Threads \
| bcftools norm -m -any --threads $Threads \
| vcflength \
| bcftools view -Oz --threads $Threads > $CacheDir/graph.select.split.vcf.gz

### Only SVs
# | bcftools view --genotype ^miss -e 'FORMAT/GT="0"' \
bcftools view -i 'length >= 50 | length <= -50' $CacheDir/graph.select.split.vcf.gz --threads $Threads \
| bcftools sort -Oz > $CacheDir/graph.select.split.SVs.vcf.gz
tabix $CacheDir/graph.select.split.SVs.vcf.gz --force

# 进一步合并 - Truvari
truvari collapse \
-i $CacheDir/graph.select.split.SVs.vcf.gz \
-o $CacheDir/truvari_merge_GraphVCF.vcf \
-c $CacheDir/truvari_collapsed_${group}.vcf
# bgzip and tabix
bcftools view $CacheDir/truvari_merge_GraphVCF.vcf --threads $Threads -Oz > $CacheDir/truvari_merge_GraphVCF.vcf.gz
tabix $CacheDir/truvari_merge_GraphVCF.vcf.gz --force

# Remove Cache file
rm $CacheDir/target_sample.txt
rm $CacheDir/truvari_merge_GraphVCF.vcf
rm $CacheDir/merge_${group}.vcf*
rm $CacheDir/truvari_collapsed_${group}.vcf*

### Stats
#bcftools +counts $CacheDir/graph.select.split.SVs.vcf.gz
#bcftools view cache/MC_graph_19.SB.split.SVs.vcf.gz |grep -v "^#"|awk '$10!= 0{print $10}'|grep -v "\."|wc -l
EOF

<<EOF
##### 对于 MC_graph_16 已经有了
###################################################
##### NGS call graph VCF
##################################################
### Only SVs, Split, sort
bcftools norm -m none $WorkDir/Cache/truvari_merge_SVs_${group}.vcf --threads $Threads \
| bcftools sort --max-mem 2G -Oz > $CacheDir/ngs_graph_${group}.vcf.gz
bcftools index --threads $Threads $CacheDir/ngs_graph_${group}.vcf.gz --force
EOF

<<EOF
#################################################
##### 1.1 Graph_VCF VS NGS_call_graph_VCF  Comp
#################################################

CompResDir=$CacheDir/FilteredSVs_${group}_comp

##### 1 Comp
# If occur error: permission denied
export TMPDIR=/pdata1/home/nizijia/tmp/
# -c $CacheDir/ngs_graph_${group}.vcf.gz \
# -c $WorkDir/Cache/merge_${group}.vcf.gz
# -b $CacheDir/graph.select.split.SVs.vcf.gz \
rm -rf ${CompResDir}
truvari bench \
-b $CacheDir/truvari_merge_GraphVCF.vcf.gz \
-c $FilteredPopVCF \
-o ${CompResDir} \
--sizemax 100000
EOF


<<EOF
#################################################
##### 1.2 Graph_VCF VS NGS_call_graph_VCF  Comp
#################################################

CompResDir=$CacheDir/${group}_comp

### Comp
# If occur error: permission denied
export TMPDIR=/pdata1/home/nizijia/tmp/
# -c $CacheDir/ngs_graph_${group}.vcf.gz \
# -c $WorkDir/Cache/merge_${group}.vcf.gz
# -b $CacheDir/graph.select.split.SVs.vcf.gz \
rm -rf $CompResDir
truvari bench \
-b $CacheDir/truvari_merge_GraphVCF.vcf.gz \
-c $PopVCF \
-o $CompResDir \
--sizemax 100000
#truvari refine $CompResDir --threads ${Threads} --reference $WorkDir/genome/Chr.BJ.fa
EOF

<<EOF
##### Extract Info
##### 1 SVLEN SVTYPE FROM AC
# 用于计算不同长度的recall情况
### NGS Call Graph
truvari anno svinfo $CompResDir/tp-base.vcf.gz \
| bcftools query -f '%CHROM\t%POS\t%length\t%SVTYPE\tNGS_CALL_Graph\t%AC\n' > $CompResDir/cache1
### Graph
#truvari anno svinfo $CacheDir/graph.select.split.SVs.vcf.gz \
truvari anno svinfo $CacheDir/truvari_merge_GraphVCF.vcf.gz \
| bcftools query -f '%CHROM\t%POS\t%length\t%SVTYPE\tGraph_Contain\t%AC\n' > $CompResDir/cache2
### Merge
cat $CompResDir/cache2 $CompResDir/cache1 > $CompResDir/GraphVCF_VS_${group}.tsv
### remove
rm $CompResDir/cache?

#####################################################
##### 2 Allele number in graph and AF(ngs called SVs)
#####################################################
for num in `seq 1 1 15`
do
### Graph VCF (AC)
bcftools view -i "AC=$num" $CacheDir/graph.select.split.SVs.vcf.gz -Oz > $CacheDir/graph.select.split.SVs.${num}.vcf.gz
bcftools index --threads $Threads $CacheDir/graph.select.split.SVs.${num}.vcf.gz --force
### Comp
# If occur error: permission denied
# Try: export TMPDIR=/pdata1/home/nizijia/tmp/
# -c $CacheDir/ngs_graph_${group}.vcf.gz \
# -c $WorkDir/Cache/merge_${group}.vcf.gz
CompResDir=$CacheDir/${group}_${num}
rm -rf $CompResDir
truvari bench \
-b $CacheDir/graph.select.split.SVs.${num}.vcf.gz \
-c $PopVCF \
-o $CompResDir \
--sizemax 100000
### Exract
truvari anno svinfo $CompResDir/tp-comp.vcf.gz \
| bcftools +fill-tags \
| bcftools query -f '%CHROM\t%POS\t%length\t%SVTYPE\tNGS_CALL_Graph\t%AC\t%AF\n' \
| awk -v var=$num '{print $0"\t"var}' > $CacheDir/GraphVCF_VS_${group}_AC_${num}.tsv
done

##### Merge
cat $CacheDir/GraphVCF_VS_${group}_AC_*.tsv > $CacheDir/GraphVCF_VS_${group}_AC.tsv

#### Remove
for num in `seq 1 1 15`
do
rm $CacheDir/GraphVCF_VS_${group}_AC_${num}.tsv
rm $CacheDir/graph.select.split.SVs.${num}.vcf.gz*
done


#for i in `seq 1 1 17`; do echo ">AC:$i";cat wild_CV_BJ_SX_${i}/log.txt|grep -E "recall|base|comp" ;done
EOF

###################################################################
###### 3 Group level VCF comp
###################################################################
# 从这一步开始group变量名变化
for group in `cat $TargetList | awk -F'\t' '{print $NF}' | uniq`
do

GroupVCFDir=$WorkDir/Cache/Script3.cache/Graph_SVs

### Comp
# If occur error: permission denied
export TMPDIR=/pdata1/home/nizijia/tmp/

CompResDir=$CacheDir/${group}_comp
rm -rf $CompResDir
truvari bench \
-c $GroupVCFDir/ngs_graph_${group}.maf.fmissing.filtered.vcf.gz \
-b $CacheDir/truvari_merge_GraphVCF.vcf.gz \
-o $CompResDir \
--sizemax 100000
### Exract
truvari anno svinfo $CompResDir/tp-comp.vcf.gz \
| bcftools +fill-tags \
| bcftools query -f '%CHROM\t%POS\t%length\t%SVTYPE\tNGS_CALL_Graph\t%AC\t%AF\t\n' \
| awk -v var=$group '{print $0"\t"var}' > $CacheDir/GraphVCF_VS_${group}_AC.tsv
done

cat $CacheDir/GraphVCF_VS_*@*_AC.tsv > $CacheDir/GraphVCF_VS_GroupGraphCallVCF.tsv
