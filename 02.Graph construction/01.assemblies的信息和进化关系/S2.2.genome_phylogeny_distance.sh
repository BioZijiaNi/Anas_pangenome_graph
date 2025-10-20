#!/bin/bash
#SBATCH --job-name=S2.2.assemblies_mash_phylogeny_tree
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S2.2.assemblies_mash_phylogeny_tree_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S2.2.assemblies_mash_phylogeny_tree_error.txt

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
### ENV
# mamba activate base
### Softawre
# mash v2.3
#######################################################################################

######################## Enviroments setting ##########################################
################## Activate environment before running ################################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S2.cache
mkdir -p $CacheDir
ResultDir=$WorkDir/Result_genome_evaluation/panduck_assemblies_mash_phylogeny_tree_autosome_ChrZ
mkdir -p $ResultDir
Threads=4
SampleList=/pdata1/home/nizijia/repository/SampleList_for_phylo_tree_autosome_chrZ.txt

################# Running Script #################
##### 01 Mash sketch
#-k <int>: k-mer length (default 21; longer k-mers increase specificity but reduce sensitivity).
#-s <int>: Sketch size (number of hash values, default 1000; larger values improve accuracy).
#-m <int>: Minimum occurrence count for each k-mer (filters low-frequency k-mers, default 1).
cat $SampleList | while read line;
do
### Vars
arr=($line)
sample=${arr[0]}
fa=${arr[1]}
### Run
mash sketch \
-s 10000 \
-k 21 \
$fa \
-o $ResultDir/${sample}.msh
done

##### 02 Merge sketch
rm $ResultDir/panduck_genome.msh
mash paste \
$ResultDir/panduck_genome \
`ls $ResultDir/*msh`

##### 03 Caculate dist
mash dist \
-s 10000 \
$ResultDir/panduck_genome.msh $ResultDir/panduck_genome.msh > $ResultDir/panduck_genome_dist.tbl
EOF

##### 04 Handle result file format
cat $ResultDir/panduck_genome_dist.tbl \
| sed 's/\/pdata1\/home\/nizijia\/repository\/1.B.2.anas_pangenome_assemblies\///g' \
| sed 's/\/pdata1\/home\/nizijia\/repository\/1.B.2.phylo_tree_out_group\///g' \
| sed 's/.fa//g' \
| sed 's/GCA_017639305.1_CAU_Laying_1.0_genomic/CAU_Laying_1_0/g' \
| sed 's/GCA_017639285.1_CAU_Pekin_2.0_genomic/CAU_Pekin_2_0/g' \
| sed 's/GCA_008746955.1_CAU_Wild1.0_genomic/CAU_Wild1_0/g' \
| sed 's/GCA_015476345.1_ZJU1.0_genomic/ZJU1_0/g' \
| sed 's/GCA_039998725.1_IASCAAS_LianchengWhiteDuck_genomic/LCWDuck/g'  > $ResultDir/panduck_genome_dist_rename_autosome_ChrZ.tbl
