#!/usr/bin/env python3

import sys
from collections import defaultdict
import argparse

def merge_intervals(input_file, output_file, group_columns, output_columns):
    # 读取输入文件
    with open(input_file, 'r') as f:
        data = f.readlines()

    # 按指定的多列分组并合并区间
    group_intervals = defaultdict(list)
    for line in data:
        parts = line.strip().split()
        # 提取分组键
        group_key = tuple(parts[int(col) - 1] for col in group_columns)  # 列索引从 1 开始
        # 提取区间信息
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        group_intervals[group_key].append((chrom, start, end, parts))  # 保存原始行数据

    # 合并每个分组的区间
    merged_results = []
    for group_key, intervals in group_intervals.items():
        # 按起始位置排序
        intervals.sort(key=lambda x: x[1])
        merged = []
        current_chrom, current_start, current_end, current_parts = intervals[0]

        for chrom, start, end, parts in intervals[1:]:
            if start <= current_end:  # 区间重叠
                current_end = max(current_end, end)
            else:
                merged.append((current_chrom, current_start, current_end, current_parts))
                current_chrom, current_start, current_end, current_parts = chrom, start, end, parts
        merged.append((current_chrom, current_start, current_end, current_parts))

        # 将合并后的区间添加到结果中
        for chrom, start, end, parts in merged:
            # 提取输出列
            output_values = [parts[int(col) - 1] for col in output_columns]
            # 将分组键和区间信息组合成一行
            group_values = " ".join(group_key)
            merged_results.append('\t'.join(output_values) + '\n')

    # 写入输出文件
    with open(output_file, 'w') as f:
        f.writelines(merged_results)

if __name__ == "__main__":
    # 使用 argparse 解析命令行参数
    parser = argparse.ArgumentParser(description="Merge intervals based on grouping columns.")
    parser.add_argument("--input", required=True, help="Input file path")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--group_columns", required=True, help="Grouping columns (comma-separated, 1-based)")
    parser.add_argument("--output_columns", required=True, help="Output columns (comma-separated, 1-based)")
    args = parser.parse_args()

    # 解析分组列和输出列
    group_columns = list(map(int, args.group_columns.split(',')))
    output_columns = list(map(int, args.output_columns.split(',')))

    # 调用合并函数
    merge_intervals(args.input, args.output, group_columns, output_columns)
    print(f"Results written to {args.output}")
