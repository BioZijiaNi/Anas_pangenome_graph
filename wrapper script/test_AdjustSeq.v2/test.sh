minimap2 -x asm5 --cs -t 2 genomeA.fa genomeB.fa > minimap2.paf
AdjustSeq.v3 --paf minimap2.paf --genomeB genomeB.fa --output seq_rc.fa --log progress.log

### Check
# 对比序列方向（原序列第一个碱基 vs 调整后最后一个碱基）
original_first=$(seqkit head -n 1 genomeA.fa | seqkit seq -s -w 0 | cut -c 1)
adjusted_last=$(seqkit head -n 1 genomeB.fa | seqkit seq -s -w 0 | rev | cut -c 1)

echo "原序列第一个碱基: $original_first"
echo "调整后最后一个碱基: $adjusted_last"
