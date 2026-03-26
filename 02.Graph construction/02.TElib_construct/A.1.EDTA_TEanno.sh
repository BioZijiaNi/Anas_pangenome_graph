################################
### ENVs
# mamba activate TEanno
### Software
# EDTA v2.2.2 (Using panEDTA.sh in the folder replace the origin file)
################################

##### Parameters #####
SampleList=/home2/nizijia/repository/1.B.5.MC_graph_16_assemblies/SampleList_whole_genome.txt
Threads=48

##### Running panEDTA
### Fix tmp
export TMPDIR=/home3/nizijia/tmp
panEDTA.sh -g $SampleList -f 2 -t $Threads
