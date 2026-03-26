#!/bin/bash

# ==============================================================================
# 局部祖源分析 (Windowed ADMIXTURE) 脚本 (参数化版本)
# ==============================================================================

# --- 1. 定义变量和参数检查 ---

# $1: PLINK 二进制文件前缀 (例如: my_data)
MAIN_BFILE="$1"
# $2: 监督模式的 FAM 文件 (例如: supervised.fam)
SUPERVISED_FAM="$2"
# $3: 结果输出目录 (例如: admixture_results_100kb)
OUT_DIR="$3"
# 窗口文件 (假设已经生成，且名称固定，或者可以作为第四个参数传入)
WINDOWS_BED="$4" 

# K值固定为 2 (可以根据需要改为参数 $4)
K_VALUE=2

# 检查参数数量
if [ "$#" -lt 3 ]; then
    echo "使用方法 (Usage):"
    echo "--------------------------------------------------------"
    echo "bash $0 <BFILE_PREFIX> <SUPERVISED_FAM> <OUTPUT_DIRECTORY>"
    echo ""
    echo "参数说明:"
    echo "  <BFILE_PREFIX>: 主 PLINK 文件的前缀 (例如: my_data，对应 my_data.bed/.bim/.fam)"
    echo "  <SUPERVISED_FAM>: 包含 K1=1, K2=2, Admixed=-9 的 .fam 文件路径"
    echo "  <OUTPUT_DIRECTORY>: 存放每个窗口 ADMIXTURE 结果的目录名称"
    echo ""
    echo "注意: 确保当前目录下存在 ${WINDOWS_BED} 文件。"
    echo "--------------------------------------------------------"
    exit 1
fi

# 检查依赖文件是否存在
if [ ! -f "${MAIN_BFILE}.bed" ] || [ ! -f "${MAIN_BFILE}.bim" ] || [ ! -f "${MAIN_BFILE}.fam" ]; then
    echo "错误: 找不到主 PLINK 文件 (${MAIN_BFILE}.bed/.bim/.fam)。请检查 <BFILE_PREFIX> 参数。"
    exit 1
fi

if [ ! -f "${SUPERVISED_FAM}" ]; then
    echo "错误: 找不到监督 FAM 文件 (${SUPERVISED_FAM})。请检查 <SUPERVISED_FAM> 参数。"
    exit 1
fi

if [ ! -f "${WINDOWS_BED}" ]; then
    echo "错误: 找不到窗口 BED 文件 (${WINDOWS_BED})。请先用 bedtools 生成。"
    exit 1
fi

# --- 2. 准备输出目录 ---
echo "创建结果输出目录: ${OUT_DIR}"
mkdir -p "$OUT_DIR"

# --- 3. 循环运行 ADMIXTURE ---
echo "--- 开始局部祖源 (K=${K_VALUE}) 分析 ---"

# 读取 100kb_windows.bed 文件的每一行
# 格式为: <chr> <start> <end>
while read -r chr start end; do

    # 为此窗口定义一个唯一的文件名
    win_name="window_${chr}_${start}_${end}"
    win_path="${OUT_DIR}/${win_name}"
    
    echo "--- 正在处理窗口: ${win_name} (${chr}:${start}-${end}) ---"

    # 1. 使用 PLINK 提取此窗口的 SNP
    plink --bfile "$MAIN_BFILE" \
          --chr "$chr" \
          --from-bp "$start" \
          --to-bp "$end" \
          --make-bed --chr-set 39 --allow-extra-chr \
          --out "$win_path" \
          --silent  # 减少 PLINK 输出信息

    # 2. 检查此窗口中是否有 SNP
    if [ ! -s "${win_path}.bim" ]; then
        echo "警告: 在 ${win_name} 中未找到 SNP (bim文件为空)。跳过。"
        # 删除生成的空文件，保持目录整洁
        rm -f ${win_path}.* continue
    fi

    # 3. 将此窗口的 .fam 文件替换为我们的监督文件
    cp "$SUPERVISED_FAM" "${win_path}.fam"

    # 4. 运行 ADMIXTURE (K=2)
    # 使用 -j 选项利用多核 (例如: -j4 使用 4 个核心，请根据你的服务器资源调整)
    # ADMIXTURE 的输出文件将是 ${win_path}.2.Q 和 ${win_path}.2.P
    admixture "${win_path}.bed" "$K_VALUE" -j4 > "${win_path}.log"

    # (可选) 清理中间文件以节省空间，只保留重要的 .Q, .P, .log
    # rm -f ${win_path}.bed ${win_path}.bim ${win_path}.fam

done < "$WINDOWS_BED"

echo "--- 所有窗口处理完毕 ---"
echo "结果文件 (.Q 和 .P) 已保存到目录: ${OUT_DIR}"
