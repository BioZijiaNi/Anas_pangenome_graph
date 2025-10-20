#!/bin/bash
#SBATCH --job-name=S2.4.sim_reads_mapping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.4_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/Script2.4_error.txt

######################## Enviroments setting ###########################################
################## Activate environment before running #################################
### ENV
# mamba activate MC_graph
# source /pdata1/home/nizijia/download/cactus-bin-v2.9.1/venv-cactus-v2.9.1/bin/activate
### Software
# vg v1.60.0
# minimap2 v2.28-r1209
### Reference
# Rice ES et al. BMC Biol. 2023. https://doi.org/10.1186/s12915-023-01758-0.
########################################################################################

###################################################################
############################### Params ############################
WorkDir="/pdata1/home/nizijia/project/01.duck_graph_pangenome"
CacheDir=$WorkDir/Cache/S2.sim_reads_mapping
mkdir -p $CacheDir
GraphDir=$WorkDir/graph/MC_graph_250314/
mkdir -p $GraphDir
GraphName="duck_MC_graph_16"
GraphType=".d2."
Threads=10

##################### YOUR RUNNING SCRIPTS #######################
##################################################################
##### 1 fastq_simulation
##################################################################
cpus=$Threads
InputGraph=$GraphDir/${GraphName}${GraphType}gbz

vg sim --progress \
    --xg-name ${InputGraph} \
    --num-reads 1000000 \
    --read-length 150 \
    --align-out \
    --random-seed 12345 \
    --sub-rate 0.0024 \
    --indel-rate 0.00029 \
    --frag-len 570 \
    --frag-std-dev 165 \
-t $Threads > $CacheDir/simulated.gam
vg view -X -a $CacheDir/simulated.gam > $CacheDir/simulated_1M.interleaved.fastq


################################################################################
##### 2 giraffe mapping
################################################################################
cpus=$Threads
fastq=$CacheDir/simulated_1M.interleaved.fastq
ResultDir=$CacheDir/sim${GraphType}reads.mapping
mkdir -p $ResultDir

for TGN in ${GraphName}.full ${GraphName} ${GraphName}.d2
do
vg giraffe -t $cpus -p \
-Z ${GraphDir}/${TGN}.gbz \
-m ${GraphDir}/${TGN}.min \
-d ${GraphDir}/${TGN}.dist \
-i -f ${fastq} \
--parameter-preset fast > $ResultDir/${TGN}.mapped.giraffe.gam

vg annotate -t $cpus -a $ResultDir/${TGN}.mapped.giraffe.gam -m -x ${GraphDir}/${TGN}.gbz \
| vg gamcompare -T -t $cpus -s -r 100 - $CacheDir/simulated.gam  > $CacheDir/${TGN}.giraffe_scores.tsv

echo -e "min_mapQ\treads_aligned\treads_aligned_correctly" \
    > $ResultDir/${TGN}.giraffe_cumulative_rates.tsv
awk 'NR>1 {aligned[$2]++} NR>1&&$1==1{correct[$2]++}
     END{OFS="\t"; for (i in aligned) print i, aligned[i], correct[i]}' \
     $CacheDir/${TGN}.giraffe_scores.tsv \
    | sort -nrk1,1 | awk '{OFS="\t"; cum_aligned+=$2; cum_correct+=$3;
         print $1, cum_aligned, cum_correct}' \
    >> $ResultDir/${TGN}.giraffe_cumulative_rates.tsv
done

################################################################################
##### 3 minimap2 mapping
################################################################################
cpus=$Threads
fastq=$CacheDir/simulated_1M.interleaved.fastq
ResultDir=$CacheDir/sim${GraphType}reads.mapping
mkdir -p $ResultDir

##### 3.1 Extract BJ path
vg paths -F -Q 'BJ#0#' -x ${WorkDir}/${GraphName}.gbz > $CacheDir/BJ.paths.fa

##### 3.2 Mapping
### 3.2.1 Mapping
# The direct alignment performed by minimap2 includes both index construction and mapping.
# This comparison is not entirely reasonable because vg giraffe uses a pre-built index.
minimap2 -ax sr -t $cpus --secondary=no \
$GraphDir/BJ.paths.fa \
$fastq \
| awk '/^@/ || $3 != "*"' > $CacheDir/mapped.no_unmapped.minimap.sam
EOF

### 3.2.2 Index and Mapping
minimap2 -x sr -t $cpus -d $CacheDir/index.mmi $GraphDir/BJ.paths.fa
minimap2 -ax sr -t $cpus --secondary=no \
$CacheDir/index.mmi \
$fastq \
| awk '/^@/ || $3 != "*"' > $CacheDir/mapped.no_unmapped.minimap.sam

vg inject -t $cpus -x ${InputGraph} $CacheDir/mapped.no_unmapped.minimap.sam > $CacheDir/mapped.minimap.gam

vg annotate -t $cpus \
    -a $CacheDir/mapped.minimap.gam \
    -m -x ${InputGraph} \
    | vg gamcompare -t $cpus \
    -s -r 100 -T - $CacheDir/simulated.gam \
    > $CacheDir/minimap_scores.tsv

echo -e "min_mapQ\treads_aligned\treads_aligned_correctly" \
    > $ResultDir/minimap_cumulative_rates.tsv
awk 'NR>1 {aligned[$2]++} NR>1&&$1==1{correct[$2]++}
     END{OFS="\t"; for (i in aligned) print i, aligned[i], correct[i]}' \
    $CacheDir/minimap_scores.tsv \
    | sort -nrk1,1 | awk '{OFS="\t"; cum_aligned+=$2; cum_correct+=$3;
         print $1, cum_aligned, cum_correct}' \
    >> $ResultDir/minimap_cumulative_rates.tsv
