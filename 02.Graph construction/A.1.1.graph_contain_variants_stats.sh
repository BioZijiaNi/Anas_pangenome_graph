#!/bin/bash
#SBATCH --job-name=A.1.1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.1_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/A.1.1_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate MC_graph
#######################################################################################

#######################################################################################
##### Setting Params #####
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/A.1.1.cache
mkdir -p $CacheDir
GraphDir=$WorkDir/graph/MC_graph_250314/
GraphName="duck_MC_graph_16"
Threads=16

################## Running script ####################################
######################################################################
##### 01 Graph VCF of different variants number
### Output: $CacheDir/${graph}.split.${type}.vcf.stats、A.1.1_out.txt
### multiqc in MC_graph exist bug
######################################################################
TargetGraph="${GraphName}.d2 ${GraphName} ${GraphName}.full"
for graph in $TargetGraph
do

echo "> $graph start:"
date
##### 1.1 Pre-processing of graph vcf
# Split
bcftools norm -m -any $GraphDir/${graph}.vcf.gz --threads $Threads \
| vcflength \
| bcftools sort --max-mem 2G -Oz > $CacheDir/${graph}.split.vcf.gz
tabix --threads $Threads --force $CacheDir/${graph}.split.vcf.gz

##### 1.2 Separate graph vcf
# SNPs (contain MNPs)
bcftools view -i 'length = 0' $CacheDir/${graph}.split.vcf.gz --threads $Threads -Oz > $CacheDir/${graph}.split.SNPs.vcf.gz
# InDel
bcftools view -i 'length < 50 | length > -50' $CacheDir/${graph}.split.vcf.gz --threads $Threads \
| bcftools view -e 'length=0' --threads $Threads -Oz > $CacheDir/${graph}.split.InDels.vcf.gz
# SVs
bcftools view -i 'length >= 50|length <= -50' $CacheDir/${graph}.split.vcf.gz --threads $Threads -Oz > $CacheDir/${graph}.split.SVs.vcf.gz

##### 1.3 Stats (in A.1.1_out.txt )
### Bcftools stats
  for type in SNPs InDels SVs
  do
  bcftools stats --threads $Threads $CacheDir/${graph}.split.${type}.vcf.gz > $CacheDir/${graph}.split.${type}.vcf.stats
  done

### Variants number
  for type in SNPs InDels SVs
  do
  Number=`bcftools view --threads $Threads $CacheDir/${graph}.split.${type}.vcf.gz|grep -v "^#"|wc -l`
  echo -e "@ ${graph}\t${type}\t${Number}"
  done
done

####################################################
##### 02 Graph VCF (SVs site length distribution)
### Output: $CacheDir/${graph}.SVs.stats
####################################################
TargetGraph="${GraphName}.d2 ${GraphName} ${GraphName}.full"
for graph in $TargetGraph
do
truvari anno svinfo $CacheDir/${graph}.split.SVs.vcf.gz \
| truvari anno lcr \
| bcftools query -f '%length.alt\t%length.ref\t%length\t%SVTYPE\tGraph_Contain\t%AC\t%AF\t%LCR\n' > $CacheDir/${graph}.SVs.stats
done

### Comp
# If occur error: permission denied
# Try: export TMPDIR=/pdata1/home/nizijia/tmp/
CompResDir=$CacheDir/SVs_Default_vs_Filtered_graph
mkdir -p $CompResDir
rm -rf $CompResDir

tabix $CacheDir/${GraphName}.d2.split.SVs.vcf.gz
tabix $CacheDir/${GraphName}.split.SVs.vcf.gz

truvari bench \
-b $CacheDir/${GraphName}.d2.split.SVs.vcf.gz \
-c $CacheDir/${GraphName}.split.SVs.vcf.gz \
-o $CompResDir \
--typeignore

####################################################
##### 03 Graph VCF stats
# Output: $CacheDir/${graph}.SVs.repmask.anno.vcf.query
####################################################
TargetGraph="${GraphName}.d2 ${GraphName} ${GraphName}.full"
for graph in $TargetGraph
do
# -p '-qq -e hmmer -species human -lcambig -nocut -div 50 -no_id -s {fasta}'
truvari anno repmask \
-i $CacheDir/${graph}.split.SVs.vcf.gz \
--min-length 50 --max-length=100000 -e RepeatMasker --threads $Threads > $CacheDir/${graph}.SVs.repmask.anno.vcf
bgzip $CacheDir/${graph}.SVs.repmask.anno.vcf --threads $Threads --force
tabix $CacheDir/${graph}.SVs.repmask.anno.vcf.gz --threads $Threads --force

bcftools view $CacheDir/${graph}.SVs.repmask.anno.vcf.gz \
| truvari anno lcr \
| truvari anno svinfo \
| bcftools query -f '%length.alt\t%length.ref\t%length\t%SVTYPE\t%LCR\t%RM_clsfam\t%RM_repeat\t%RM_score\n' > $CacheDir/${graph}.SVs.repmask.anno.vcf.query
done

##################################################################
##### 04 Variants source
# Output: $CacheDir/${graph}.split.SVs.variants.source.Header.txt
##################################################################
TargetGraph="${GraphName}.d2 ${GraphName} ${GraphName}.full"
for graph in $TargetGraph
do
echo "> $graph"

bcftools view $CacheDir/${graph}.split.SVs.vcf.gz \
| bcftools query -l \
| tr '\n' '\t' \
| awk '{print "Length\tAC\tGT\t"$0}' > $CacheDir/samplename.txt

bcftools view $CacheDir/${graph}.split.SVs.vcf.gz \
| bcftools query -f '%length\t%AC\t%FORMAT\n' > $CacheDir/${graph}.split.SVs.variants.source.txt

cat \
$CacheDir/samplename.txt \
$CacheDir/${graph}.split.SVs.variants.source.txt > $CacheDir/${graph}.split.SVs.variants.source.Header.txt

rm $CacheDir/samplename.txt
rm $CacheDir/${graph}.split.SVs.variants.source.txt
done
EOF

###################################################################
##### 05 Share and common variants
# Output: $CacheDir/${graph}.split.SVs.variants.source.Header.txt
###################################################################
TargetGraph="${GraphName}.d2 ${GraphName} ${GraphName}.full"
for graph in $TargetGraph
do
echo "> $graph"

bcftools view $CacheDir/${graph}.split.SVs.vcf.gz \
| bcftools query -l \
| tr '\n' '\t' \
| awk '{print "CHROM\tPOS\tLength\tAC\tGT\t"$0}' > $CacheDir/samplename.txt

bcftools view $CacheDir/${graph}.split.SVs.vcf.gz \
| bcftools query -f '%CHROM\t%POS\t%length\t%AC\t%FORMAT\n' > $CacheDir/${graph}.split.SVs.variants.source.txt

cat \
$CacheDir/samplename.txt \
$CacheDir/${graph}.split.SVs.variants.source.txt > $CacheDir/${graph}.split.SVs.variants.txt

rm $CacheDir/samplename.txt
rm $CacheDir/${graph}.split.SVs.variants.source.txt
done
