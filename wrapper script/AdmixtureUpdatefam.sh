#!/bin/bash

# ==============================================================================
# 脚本功能: 根据 ADMIXTURE Q 文件结果更新 PLINK FAM 文件以启用监督模式
# 使用方法: bash update_fam_with_params.sh <Q_FILE> <ORIGINAL_FAM> <OUTPUT_FAM> <THRESHOLD>
# ==============================================================================

# --- 1. 参数接收与检查 ---

Q_FILE="$1"          # 第 1 个参数：祖源比例文件 (.Q 文件)
ORIGINAL_FAM="$2"    # 第 2 个参数：原始 PLINK .fam 文件
SUPERVISED_FAM="$3"  # 第 3 个参数：输出的监督模式 .fam 文件
THRESHOLD="$4"       # 第 4 个参数：祖源比例阈值 (例如 0.95)

# 检查参数数量
if [ "$#" -ne 4 ]; then
    echo "使用方法 (Usage):"
    echo "--------------------------------------------------------------------------------------"
    echo "bash $0 <Q_FILE> <ORIGINAL_FAM> <OUTPUT_FAM> <THRESHOLD>"
    echo ""
    echo "参数说明:"
    echo "  <Q_FILE>: 包含祖源比例的文件 (例如: my_data.2.Q 或您提供的 admixture_q_file.txt)"
    echo "  <ORIGINAL_FAM>: 原始的 PLINK .fam 文件 (第六列通常全为 -9)"
    echo "  <OUTPUT_FAM>: 生成的监督模式 .fam 文件名 (例如: supervised.fam)"
    echo "  <THRESHOLD>: 确定亲本的祖源比例阈值 (例如: 0.95, 表示 Q1 或 Q2 超过此值则固定)"
    echo "--------------------------------------------------------------------------------------"
    exit 1
fi

# 检查依赖文件是否存在
if [ ! -f "$Q_FILE" ]; then
    echo "错误: 找不到祖源比例文件 (${Q_FILE})。"
    exit 1
fi

if [ ! -f "$ORIGINAL_FAM" ]; then
    echo "错误: 找不到原始 FAM 文件 (${ORIGINAL_FAM})。"
    exit 1
fi

# 检查阈值是否是有效的数字 (简单检查)
if ! [[ "$THRESHOLD" =~ ^[0-9]*(\.)?[0-9]+$ ]]; then
    echo "错误: 阈值 (${THRESHOLD}) 必须是有效的数字。"
    exit 1
fi

echo "--- 开始更新 FAM 文件 ---"
echo "祖源比例文件 (Q): ${Q_FILE}"
echo "原始 FAM 文件:    ${ORIGINAL_FAM}"
echo "阈值 (Threshold): ${THRESHOLD}"
echo "输出文件:         ${SUPERVISED_FAM}"
echo "---------------------------"

# --- 2. AWK 核心逻辑 ---

# AWK 脚本：
# 1. 以阈值 ($THRESHOLD) 作为 AWK 变量 THRESH 传入。
# 2. 首先读取 Q 文件，将 K1, K2 比例映射到个体 ID。
# 3. 然后读取原始 FAM 文件，根据 K1/K2 比例修改第 6 列。
awk -v THRESH="$THRESHOLD" \
    'BEGIN { OFS="\t" } 
     
     # 第一遍处理 Q 文件 (Q_FILE)
     # 文件记录编号 (FNR) 等于总记录编号 (NR) 时，表示正在处理第一个文件 (Q_FILE)
     FNR==NR {
         # $1: Individual ID
         # $2: K1 (Parent 1) proportion
         # $3: K2 (Parent 2) proportion
         K1[$1] = $2
         K2[$1] = $3
         next
     }

     # 第二遍处理 FAM 文件 (ORIGINAL_FAM)
     {
         # PLINK FAM 文件中个体 ID 是第二列 ($2)
         ind_id = $2
         
         # 默认值保持为原始的第六列的值 (通常是 -9)
         new_pheno = $6 
         
         # 检查该 ID 是否有祖源比例数据
         if (ind_id in K1) {
             k1_prop = K1[ind_id]
             k2_prop = K2[ind_id]

             # 逻辑 1: K1 比例 >= 阈值 -> 设为 1 (Parent 1)
             if (k1_prop >= THRESH) {
                 new_pheno = 1
             } 
             # 逻辑 2: K2 比例 >= 阈值 -> 设为 2 (Parent 2)
             else if (k2_prop >= THRESH) {
                 new_pheno = 2
             }
             # 否则 (两者都低于阈值)，保持 new_pheno = -9
         }

         # 替换第 6 列并打印新行
         $6 = new_pheno
         print
     }' "$Q_FILE" "$ORIGINAL_FAM" > "$SUPERVISED_FAM"

echo "--- FAM 文件更新完毕，已保存到 ${SUPERVISED_FAM} ---"
