#!/bin/bash

#################################################
# Date： 2025/05
# Author： ZijiaNi + DeepSeek R1 Thinking mode
#################################################

# ================= 参数说明 =================
# 功能：封装序列调整流程，输入参考和查询文件，输出调整后的序列
# 示例：./script.sh --ref ref.fa --query query.fa -o output -t 8

usage() {
  echo "用法: $0 -r|--ref <参考文件> -q|--query <查询文件> -o <输出前缀> [-t 线程数]"
  echo "选项:"
  echo "  -r, --ref   参考基因组文件（必须）"
  echo "  -q, --query 待处理的查询文件（必须）"
  echo "  -o          输出前缀（生成 [前缀].adjust.fa 和 [前缀].log）"
  echo "  -t          线程数（默认: 4）"
  exit 1
}

# ================= 参数解析 =================
Threads=4
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--ref)
      RefFile="$2"
      shift 2
      ;;
    -q|--query)
      QueryFile="$2"
      shift 2
      ;;
    -o)
      OutputPrefix="$2"
      shift 2
      ;;
    -t)
      Threads="$2"
      shift 2
      ;;
    *)
      usage
      ;;
  esac
done

# 检查必要参数
[ -z "${RefFile:-}" ] && { echo "错误：缺少 --ref 参数"; usage; }
[ -z "${QueryFile:-}" ] && { echo "错误：缺少 --query 参数"; usage; }
[ -z "${OutputPrefix:-}" ] && { echo "错误：缺少 -o 参数"; usage; }
[ ! -f "$RefFile" ] && { echo "错误：参考文件不存在 $RefFile"; exit 1; }
[ ! -f "$QueryFile" ] && { echo "错误：查询文件不存在 $QueryFile"; exit 1; }

# ================= 主流程 =================
echo "==== 开始处理: $(basename "$QueryFile") ===="

# 临时中间文件（自动清理）
TMP_CHR="${OutputPrefix}.tmp.chr.fa"
TMP_PAF="${OutputPrefix}.tmp.paf"

# Step 1: 排序并提取ID
echo "1. 处理查询文件..."
seqkit sort -n -N -j "$Threads" "$QueryFile" | \
  seqkit seq --only-id -j "$Threads" > "$TMP_CHR"

# Step 2: minimap2比对
echo "2. 运行 minimap2..."
minimap2 -x asm5 --cs -t "$Threads" \
  "$RefFile" \
  "$TMP_CHR" \
  > "$TMP_PAF"

# Step 3: 调整序列
echo "3. 运行 AdjustSeq..."
AdjustSeq.v3 \
  --paf "$TMP_PAF" \
  --genomeB "$TMP_CHR" \
  --output "${OutputPrefix}.adjust.fa" \
  --log "${OutputPrefix}.log"

# 清理中间文件
echo "4. 清理临时文件..."
rm -f "$TMP_CHR" "$TMP_PAF"

echo "==== 处理完成 → 输出文件: ${OutputPrefix}.adjust.fa 和 ${OutputPrefix}.log ===="
