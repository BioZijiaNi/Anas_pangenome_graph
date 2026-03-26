#!/bin/bash
#SBATCH --job-name=S1.3.synteny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --output=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.3.genome_assessment_synteny_out.txt
#SBATCH --error=/pdata1/home/nizijia/project/01.duck_graph_pangenome/log/S1.3.genome_assessment_synteny_error.txt


################## Enviroments setting #######################################################

################## 运行前提前激活环境 ########################################################
### ENV
# mamba activate genome_evaluate_synteny
### Softawre
# plotsr == 1.1.1
# syri == 1.7.0
# minimap2 == 2.28-r1209
# /pdata1/home/nizijia/download/NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl == 1.41
#############################################################################################



################## Setting input and output directory ###############
########## Params ##########
WorkDir="/home1/nizijia/cache"
CacheDir="/home1/nizijia/tmp"
mkdir -p $CacheDir
ResultDir=$WorkDir/genome_synteny_16_assemblies
mkdir -p $ResultDir
Threads=20

# 核心修改：定义输入的表格文件路径
# 表格格式：样本名  路径 (空格或Tab分隔)
#SampleTable="/home2/nizijia/repository/1.B.5.MC_graph_16_assemblies/SampleList_chrom_genome.txt"
SampleTable="/home2/nizijia/repository/1.B.5.MC_graph_16_assemblies/SampleList_chrom_genome.reversed.reorder.txt"

# 1. 读取表格到数组
samples=()
paths=()

# 使用 while 循环读取，去除 "Chr." 前缀（如果需要）
while read -r name path; do
    # 这里的 name 会变成 BJ, CAU_Laying_1.0 等
    clean_name=$(echo $name | sed 's/Chr.//')
    samples+=("$clean_name")
    paths+=("$path")
done < "$SampleTable"

num_samples=${#samples[@]}
echo "Loaded $num_samples genomes from table."

AFDir=$ResultDir/adjusted_fasta
mkdir -p $AFDir

################# Running script ###################
##### 1.1 动态提取染色体并处理
echo "Running 1.1 Extract target chromosomes"

for i in $(seq 0 $((num_samples-1))); do
    s_name=${samples[$i]}
    s_path=${paths[$i]}

    echo "Processing $s_name ..."

    # 1. 自动从 .fai 索引中抓取符合规律的 ID
    # 逻辑：查找以 Chr 开头，中间有数字或 Z/W，并带有样本后缀的 ID
    # 比如：Chr01_Ma, Chr32_Z_Ma 等
    if [ ! -f "${s_path}.fai" ]; then samtools faidx "$s_path"; fi

    # 提取常染色体 1-30 和 性染色体 Z/W 的实际 ID
    # 这里使用 grep 匹配：Chr(数字)_(后缀) 或 Chr(Z/W)_(后缀) 或 Chr(32_Z)_(后缀)
    target_ids=$(cut -f1 "${s_path}.fai" | grep -E "Chr([0-9]+|Z|W|32_Z|33_W)_$s_name" | xargs)

    if [ -z "$target_ids" ]; then
        echo "Warning: No chromosomes found for $s_name using pattern Chr.*_$s_name"
        # 备选方案：如果不带后缀，尝试直接提取 Chr01 等
        target_ids=$(cut -f1 "${s_path}.fai" | grep -E "Chr([0-9]+|Z|W|32_Z|33_W)$" | xargs)
    fi

    # 2. 提取、更名、补零
    # seqkit replace 的逻辑：
    # - 第一个：去掉下划线及其后的样本名 (Chr01_Ma -> Chr01)
    # - 第二个：处理特殊的 Chr32_Z / Chr33_W (Chr32_Z -> ChrZ)
    # - 第三个：将 Chr1 这种单位数补全为 Chr01 (如果已经是 Chr01 则保持不变)
    samtools faidx "$s_path" $target_ids | \
    seqkit replace -p "_${s_name}$" -r "" | \
    seqkit replace -p "^Chr32_Z" -r "ChrZ" | \
    seqkit replace -p "^Chr33_W" -r "ChrW" | \
    seqkit replace -p "Chr(\d)$" -r "Chr0\${1}" \
    > $CacheDir/${s_name}_cache.fa

    echo "  $s_name done. Output: $CacheDir/${s_name}_cache.fa"
done

##### 1.2 Reorder and reverse (以数组第一个基因组为参考)
echo "Running 1.2 Reorder"
ref_name=${samples[0]}
# 第一个基因组直接作为参考
cat $CacheDir/${ref_name}_cache.fa > $AFDir/${ref_name}.adjust.fa

# 从第二个开始进行 FnaReorder
# 注意：原脚本这里写的是 seq 6，我将其改回了 seq 1 (从第2个样本开始)，请根据实际需求调整
for i in $(seq 1 $((num_samples-1))); do
    curr_name=${samples[$i]}

    # 如果文件已存在且非空，可以选择跳过（可选）
    if [ -s "$AFDir/${curr_name}.adjust.fa" ]; then
        echo "  $AFDir/${curr_name}.adjust.fa exists, skipping reorder..."
    else
        FnaReorder.sh --ref $AFDir/${ref_name}.adjust.fa \
                      --query $CacheDir/${curr_name}_cache.fa \
                      -o $AFDir/${curr_name} -t $Threads
    fi
done


##### 1.2.5 [新增] 提取指定染色体 (生成 .final.fa)
# 这一步只提取 Chr01-Chr30, ChrZ，并且不覆盖 adjust.fa
echo "Running 1.2.5 Subsampling target chromosomes to .final.fa"

# 正则逻辑：
# ^Chr       : 以 Chr 开头
# (0[1-9]    : 匹配 01-09
# |[12][0-9] : 匹配 10-29
# |30        : 匹配 30
# |Z|W)      : 匹配 Z 或 W
# $          : 字符串结尾 (防止匹配到 Chr1_random 等)
Target_Regex="^Chr(0[1-9]|[12][0-9]|Z)$"

for i in $(seq 0 $((num_samples-1))); do
    s_name=${samples[$i]}
    echo "  Subsampling $s_name ..."

    # 输入用 adjust.fa，输出为 final.fa
    seqkit grep -r -p "$Target_Regex" $AFDir/${s_name}.adjust.fa > $AFDir/${s_name}.final.fa

    # 建立索引，确保后续程序能正确读取长度
    samtools faidx $AFDir/${s_name}.final.fa
done


##### 1.3 Create Pairwise List
# 生成相邻两两比对的列表
> $AFDir/list
for i in $(seq 0 $((num_samples-2))); do
    echo -e "${samples[$i]}\t${samples[$((i+1))]}" >> $AFDir/list
done

##### 1.6 Running NGenomeSyn Pairwise Alignment
echo "Running 1.6 NGenomeSyn Alignments"
cat $AFDir/list | while read s1 s2; do
    echo "  Aligning $s1 vs $s2 ..."
    perl /home/nizijia/download/NGenomeSyn-1.41/bin/GetTwoGenomeSyn.pl \
    -InGenomeA $AFDir/${s1}.final.fa \
    -InGenomeB $AFDir/${s2}.final.fa \
    -MinLenA 100000 -MinLenB 100000 -MinAlnLen 500000 \
    -MappingBin minimap2 -BinDir minimap2 \
    -OutPrefix $CacheDir/${s1}_vs_${s2} -NumThreads $Threads
done


################# 1.7 自动化生成配置文件 #################
echo "Generating NGenomeSyn configuration file..."
Config_File="$ResultDir/Pangenome_16_assemblies.conf"

# --- 核心布局参数（调整这里即可控制紧凑度）---
Start_Y=100       # 第一个基因组距离顶部的距离
L_Width=40        # 【核心】Linker 的高度。设为 60-80 会非常紧凑，Linker 变短
C_Width=8        # 染色体线条的粗细
# ------------------------------------------

# 自动计算画布总高度 (对应参数 body)
# 公式：起始位置 + (所有基因组层数 * 间距) + 底部留白
num_samples=${#samples[@]}
Dynamic_Body=$(bc <<< "$Start_Y + ($num_samples * $L_Width) + 150")

# 1. 写入 Global Setting
cat << EOF > $Config_File
SetParaFor = global
body = $Dynamic_Body       # 对应你提供的参数，动态调整画布高度
CanvasWidth = 1200         # 对应你提供的参数
# 如果需要调整左右边距，可以取消下面注释
# left = 100
# right = 120
# up = 55
EOF

# 2. 写入 GenomeInfoFile (16个)
for i in $(seq 0 $((num_samples-1))); do
    idx=$((i+1))
    if [ $i -eq 0 ]; then
        echo "GenomeInfoFile${idx}=$CacheDir/${samples[0]}_vs_${samples[1]}.A.chr.len" >> $Config_File
    else
        echo "GenomeInfoFile${idx}=$CacheDir/${samples[$((i-1))]}_vs_${samples[$i]}.B.chr.len" >> $Config_File
    fi
done

# 3. 写入 LinkFile (15组)
for i in $(seq 0 $((num_samples-2))); do
    idx1=$((i+1))
    idx2=$((i+2))
    echo "LinkFileRef${idx1}VsRef${idx2}=$CacheDir/${samples[$i]}_vs_${samples[$((i+1))]}.link" >> $Config_File
done

# 4. 写入 Link 样式（灰色）
cat << EOF >> $Config_File

SetParaFor = LinkALL
fill = grey
stroke = grey
fill-opacity = 0.5
stroke-opacity = 0.5
OutsideWidthDeta = 1    # 隐藏参数，防止灰色盖住染色体边框
EOF

# 5. 写入 16 个 Genome 的详细坐标设置
for i in $(seq 0 $((num_samples-1))); do
    idx=$((i+1))
    # 计算每个基因组确切的 Y 轴位置
    current_y=$(bc <<< "$Start_Y + $i * $L_Width")

    cat << EOF >> $Config_File
SetParaFor = Genome${idx}
MoveToY = $current_y    # 强制移动到计算出的 Y 坐标
ChrWidth = $C_Width     # 染色体高度
LinkWidth = $L_Width    # 与下一个基因组之间的连线宽度
NormalizedScale = 1
GenomeName = ${samples[$i]}
LabelFontSize = 10      # 16个基因组建议每层单独设小一点字体

EOF
done

echo "Done! Config saved to $Config_File (body=$Dynamic_Body, LinkWidth=$L_Width)"

################# 1.8 运行绘图 #################
/home/nizijia/download/NGenomeSyn-1.41/bin/NGenomeSyn \
-InConf $Config_File \
-OutPut $ResultDir/Pangenome_16_assemblies.pdf
