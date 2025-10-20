#!/bin/bash
#SBATCH --job-name=S1.2.assembly_stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.2.genome_assessment_assembly_stats_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.2.genome_assessment_assembly_stats_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate base
### Softawre
# assembly-stats
# bedtools
#######################################################################################

#######################################################################################
##### Setting Params #####
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache
mkdir -p $CacheDir
ResultDir=$WorkDir/Result_genome_evaluation/

Threads=4
SampleList=$WorkDir/data/genome_list

################## Running script ##################
####################################################
##### 01 Chromosome-level assembly stats
####################################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fa=${arr[1]}

### Input and output
genome=$fa
out_dir=$ResultDir/${sample}_assembly_stats

### Running
assembly-stats -t  $fa > $out_dir/assembly_stats
done

####################################################
##### 02 Contig-level assembly stats
####################################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fa=${arr[1]}

### Input and output
genome=$fa
out_dir=$ResultDir/${sample}_assembly_stats
mkdir -p $out_dir

### Running
# 1. Extract Gap region
samtools faidx $genome
seqkit locate -i -p '"N+"' -r $genome --non-greedy --bed -P --threads $Threads > $out_dir/gaps.bed
# 2. Break N+ and generate contigs
bedtools complement -i $out_dir/gaps.bed -g <(cut -f1,2 ${genome}.fai) > $out_dir/non_gaps.bed
bedtools getfasta -fi $genome -bed $out_dir/non_gaps.bed -fo $out_dir/contigs.fasta -name
# Remove cache
rm $out_dir/non_gaps.bed
# 3. stats
assembly-stats -t  $out_dir/contigs.fasta > $out_dir/contig_stats
# 4 Enzyme site in gap
bedtools complement \
-i $out_dir/gaps.bed \
-g ${genome}.fai \
| awk '{print $0"\t"$3-$2}' | awk '$4 <= 6 {print $0}' \
| bedtools getfasta -fi $genome -bed - > $out_dir/gaps_contain_enzyme_sites.fa
done
EOF


############################################################################
##### 03 Contig-level assembly stats (After merge bionano DLE-1 enzyme site)
############################################################################
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fa=${arr[1]}

### Input and output
genome=$fa
out_dir=$ResultDir/${sample}_assembly_stats
mkdir -p $out_dir

# 1. convert seq to one line, match CTTAGG in N
seqkit seq -w 0 $genome \
| perl -pe 's/(?<=[Nn])(CTTAAG)(?=[Nn])/"N" x length($1)/gei' \
| seqkit seq -w 60 > $out_dir/${sample}_merged_enzyme_site.fa
# 2. Extract gap region
samtools faidx $out_dir/${sample}_merged_enzyme_site.fa
seqkit locate -i -p '"N+"' -r $out_dir/${sample}_merged_enzyme_site.fa --non-greedy --bed -P --threads $Threads > $out_dir/gaps_merged_enzyme_site.bed
# 3. Break N+ and generate contigs
bedtools complement -i $out_dir/gaps_merged_enzyme_site.bed -g <(cut -f1,2 $out_dir/${sample}_merged_enzyme_site.fa.fai) > $out_dir/non_gaps.bed
bedtools getfasta -fi $out_dir/${sample}_merged_enzyme_site.fa -bed $out_dir/non_gaps.bed -fo $out_dir/contigs_merged_enzyme_site.fasta -name
# Remove cache
rm $out_dir/non_gaps.bed
# 4. stats
assembly-stats -t  $out_dir/contigs_merged_enzyme_site.fasta > $out_dir/contig_merged_enzyme_site.stats
done


###########################################################################
##### 04. Contig Length for Treemap
##########################################################################
cat $SampleList | while read line;
do
# Vars
arr=($line)
sample=${arr[0]}
fa=${arr[1]}

genome=$fa
out_dir=$ResultDir/${sample}_assembly_stats
mkdir -p $out_dir

seqkit fx2tab $out_dir/contigs_merged_enzyme_site.fasta --name --length -j $Threads \
| awk -v var=$sample '{print $1"\t"$2"\t"var}'> $out_dir/${sample}_contigs_merged_enzyme_site_name_length.txt
done
