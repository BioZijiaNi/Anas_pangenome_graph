#!/bin/bash

# 脚本名称: update_fam_pheno.sh
# 描述: 根据源FAM文件的第6列（表型），更新目标FAM文件的第6列。
# 用法: ./update_fam_pheno.sh <源FAM文件> <目标FAM文件> <输出FAM文件>

# --- 检查参数数量 ---
if [ "$#" -ne 3 ]; then
    echo "使用错误！"
    echo "用法: $0 <源FAM文件> <目标FAM文件> <输出FAM文件>"
    echo "  <源FAM文件>: 包含新的第6列表型的文件 (例如 allsnp_filter.supervised.fam)"
    echo "  <目标FAM文件>: 需要更新第6列表型的文件 (例如 allsnp_filter.fam)"
    echo "  <输出FAM文件>: 更新后的新文件名"
    exit 1
fi

FAM_SOURCE="$1"
FAM_TARGET="$2"
FAM_OUTPUT="$3"

# --- 检查文件是否存在 ---
if [ ! -f "$FAM_SOURCE" ] || [ ! -f "$FAM_TARGET" ]; then
    echo "错误: 至少有一个输入文件不存在。"
    echo "请检查文件路径: $FAM_SOURCE 和 $FAM_TARGET"
    exit 1
fi

echo "--- 正在执行表型更新 ---"
echo "源文件 (提供新表型): $FAM_SOURCE"
echo "目标文件 (待更新): $FAM_TARGET"
echo "输出文件: $FAM_OUTPUT"

# --- AWK 核心逻辑 ---
# FS和OFS都使用空格/制表符的混合模式（默认行为），以确保能正确处理您的输入格式
awk '
BEGIN {
    # 启用混合分隔符模式，适用于大多数FAM文件
}
# 1. 读取源文件 (FAM_SOURCE)
FNR==NR {
    # 将 $1 (FID) 和 $2 (IID) 作为键，第6列 ($6) 作为值存储
    P[$1 " " $2] = $6; 
    next
}
# 2. 读取目标文件 (FAM_TARGET)
{
    key = $1 " " $2;
    if (key in P) {
        # 如果样本ID匹配，则用源文件中的新值更新第6列
        $6 = P[key];
    }
    
    # 打印整行。由于我们没有明确设置OFS，awk会使用其默认的单个空格或tab来分隔字段，
    # 确保了输出格式是标准的fam格式。
    print
}' "$FAM_SOURCE" "$FAM_TARGET" > "$FAM_OUTPUT"

echo "--- 更新完成！新文件已保存至 $FAM_OUTPUT ---"
