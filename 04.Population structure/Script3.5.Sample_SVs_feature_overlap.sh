#!/bin/bash
#SBATCH --job-name=S3.5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.5_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script3.5_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate base
#######################################################################################

#######################################################################################
############################## Setting params #########################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/Script3.cache/Graph_SVs_overlap_genome_feature
mkdir -p $CacheDir
ResultDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_vg_pipeline
mkdir -p $ResultDir
RSDir=$WorkDir/Cache/Script1.vg.giraffe.varirant.pipeline/Result_linear_small_variants
mkdir -p $RSDir
GraphDir=$WorkDir/graph/MC_graph_250314
mkdir -p $GraphDir

GFDir=$WorkDir/genome
GFF=$GFDir/Chr.BJ.fix.ChromName.gff
FA=/pdata1/home/nizijia/repository/1.B.2.anas_pangenome_assemblies/Chr.BJ.fa
Length=5000

GraphName="duck_MC_graph_16.d2"
Threads=16

SampleList="/pdata1/home/nizijia/repository/1.B.6.ncbi_duck_WGS/script/Table1.duck_WGS_SRR_code.tsv"


################## Running script ###################
##########################################
##### 01. Create Genome Feature File #####
##########################################
### 5kb upstream (considered strand +-)
cat $GFF \
| awk '$3=="mRNA" {print $0}' \
| bedtools flank -i - -g ${FA}.fai -l $Length -r 0 -s \
| awk '$3 == "mRNA" {
    # Extract gene from COL9
    if (match($9, /gene "([^"]+)"/, arr)) {
        gene = arr[1]
    } else if (match($9, /gene=([^;]+)/, arr)) {
        gene = arr[1]
    }
    else {
        gene = "Unknown"
    }

    # Output bed
    print $1"\t"$4 - 1"\t"$5"\t"gene"\t5kb_upstream\t"$7
}' > $CacheDir/cache.bed

merge_intervals.py \
--input $CacheDir/cache.bed \
--output $CacheDir/5kb_upstream.bed \
--group_columns 1,4,5,6 \
--output_columns 1,2,3,4,5,6

### 5kb downstream (considered strand +-)
cat $GFF \
| awk '$3=="mRNA" {print $0}' \
| bedtools flank -i - -g ${FA}.fai -l 0 -r $Length -s \
| awk '$3 == "mRNA" {
    # Extract gene from COL9
    if (match($9, /gene "([^"]+)"/, arr)) {
        gene = arr[1]
    } else if (match($9, /gene=([^;]+)/, arr)) {
        gene = arr[1]
    }
    else {
        gene = "Unknown"
    }

    # Output bed
    print $1"\t"$4 - 1"\t"$5"\t"gene"\t5kb_downstream\t"$7
}' > $CacheDir/cache.bed

merge_intervals.py \
--input $CacheDir/cache.bed \
--output $CacheDir/5kb_downstream.bed \
--group_columns 1,4,5,6 \
--output_columns 1,2,3,4,5,6

### gene-exon
cat $GFF \
| awk '$3 == "exon" {
    #  Extract gene_name from COL9
    if (match($9, /gene_name "([^"]+)"/, arr)) {
        gene_name = arr[1]
    } else if (match($9, /gene_name=([^;]+)/, arr)) {
        gene_name = arr[1]
    }
    else {
        gene_name = "Unknown"

    # Output bed
    print $1"\t"$4 - 1"\t"$5"\t"gene_name"\texon\t"$7
}' > $CacheDir/exon.bed
### gene
cat $GFF \
| awk '$3 == "gene" {
    if (match($9, /gene "([^"]+)"/, arr)) {
        gene = arr[1]
    } else if (match($9, /gene=([^;]+)/, arr)) {
        gene = arr[1]
    }
    else {
        gene = "Unknown"
    }

    print $1"\t"$4 - 1"\t"$5"\t"gene"\tgene\t"$7
}' > $CacheDir/gene.bed
### intergenic
# Some genes in the annotation file are overlapping, which is why there are ultimately only 14,532 intergenic regions.
# cat Chr.BJ.gff | awk '$3=="gene" {print $0}'|bedtools sort -i - -g Chr.BJ.fa.fai|bedtools cluster -i -|tail
cat $GFF | awk '$3=="gene" {print $0}' \
| bedtools sort -i - -g ${FA}.fai \
| bedtools complement -i - -g ${FA}.fai |awk '{print $0"\tNA\tintergenic\tNA"}' > $CacheDir/intergenic.bed
### merge
cat \
$CacheDir/5kb_upstream.bed \
$CacheDir/5kb_downstream.bed \
$CacheDir/exon.bed \
$CacheDir/gene.bed \
$CacheDir/intergenic.bed > $GFDir/Chr.BJ.genome.feature.bed
# Remove
rm $CacheDir/5kb_upstream.bed $CacheDir/5kb_downstream.bed $CacheDir/exon.bed $CacheDir/gene.bed $CacheDir/intergenic.bed


####################################
##### 02. SVs in genome feature
####################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
group=${arr[4]}
##### 2.1 overlapping per sample 
echo "> $sample start merge:"
date
### Overlaping
bedtools intersect \
-b ${ResultDir}/${GraphName}-${sample}_clean_split_PASS_Q30_SVs.vcf.gz \
-a $WorkDir/genome/Chr.BJ.genome.feature.bed -c \
| awk -v var=$sample '{print $0"\t"var}' \
| awk -v var=$group '{print $0"\t"var}' > $CacheDir/${sample}_genome_feature.txt
done
### Merge
cat $CacheDir/SRR*_genome_feature.txt > $CacheDir/merged_genome_feature.txt
gzip $CacheDir/merged_genome_feature.txt -f
### Remove cache
rm $CacheDir/SRR*_genome_feature.txt
