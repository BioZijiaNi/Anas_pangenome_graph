#!/usr/bin/env python3

import os
import argparse
import subprocess
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor

# 定义 repmask.py 新增的 Header INFO 字段
RM_HEADER_EXTRAS = [
    '##INFO=<ID=RM_REF_score,Number=1,Type=Integer,Description="RepMask bit score for REF sequence">',
    '##INFO=<ID=RM_REF_pdiv,Number=1,Type=Float,Description="RepMask % divergence for REF sequence">',
    '##INFO=<ID=RM_REF_repeat,Number=1,Type=String,Description="RepMask matching repeat for REF sequence">',
    '##INFO=<ID=RM_REF_clsfam,Number=1,Type=String,Description="RepMask repeat class/family for REF sequence">',
    '##INFO=<ID=RM_ALT_score,Number=1,Type=Integer,Description="RepMask bit score for ALT sequence">',
    '##INFO=<ID=RM_ALT_pdiv,Number=1,Type=Float,Description="RepMask % divergence for ALT sequence">',
    '##INFO=<ID=RM_ALT_repeat,Number=1,Type=String,Description="RepMask matching repeat for ALT sequence">',
    '##INFO=<ID=RM_ALT_clsfam,Number=1,Type=String,Description="RepMask repeat class/family for ALT sequence">'
]

def get_args():
    parser = argparse.ArgumentParser(description="并行拆分运行 repmask.py (支持 -M 最大长度限制)")
    # 基础输入输出
    parser.add_argument("-i", "--input", required=True, help="输入 VCF")
    parser.add_argument("-o", "--output", required=True, help="最终合并的 VCF")
    parser.add_argument("-w", "--workdir", required=True, help="临时文件存放目录")
    
    # 脚本路径
    parser.add_argument("--repmask_path", default="./repmask.py", help="repmask.py 脚本路径")
    
    # 资源分配
    parser.add_argument("-J", "--jobs", type=int, default=4, help="并行进程数")
    parser.add_argument("-T", "--threads_per_job", type=int, default=8, help="单任务线程数")
    parser.add_argument("-l", "--lines", type=int, default=5000, help="每块变异行数")
    
    # 传递给子脚本的长度控制 (新增 -M)
    parser.add_argument("-m", "--min-length", default="50", help="最小长度限制")
    parser.add_argument("-M", "--max-length", default="100000", help="最大长度限制 (默认 100kb)")
    
    # 其他子脚本参数
    parser.add_argument("-t", "--threshold", default="0.4")
    parser.add_argument("-p", "--params", 
                        default="-qq -e NCBI -lib /home2/nizijia/repository/1.B.6.TElib/Final_TElib/Final_TElib_3.fa -lcambig -nocut -div 50 -no_id -s {fasta} -pa {threads} -frag 100000")
    
    parser.add_argument("--keep-cache", action="store_true")
    return parser.parse_args()

def run_repmask(vcf_part, args):
    """调用子脚本处理分块"""
    output_part = vcf_part + ".anno"
    current_p = args.params.replace("{threads}", str(args.threads_per_job))
    script_abs_path = os.path.abspath(args.repmask_path)
    
    # 构造命令：加入 -m 和 -M 参数
    cmd = (
        f"python {script_abs_path} -i {vcf_part} -o {output_part} "
        f"-m {args.min_length} -M {args.max_length} "
        f"-p '{current_p}' -T {args.threads_per_job} -t {args.threshold} "
        f"-w {args.workdir}"
    )
    
    try:
        env = os.environ.copy()
        env["TMPDIR"] = os.path.abspath(args.workdir)
        # 执行子脚本
        subprocess.run(cmd, shell=True, check=True, capture_output=True, env=env)
        return output_part
    except subprocess.CalledProcessError as e:
        print(f"[错误] 分块 {vcf_part} 失败: {e.stderr.decode()}")
        return None

def main():
    args = get_args()
    work_dir = os.path.abspath(args.workdir)
    
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    
    try:
        # 1. 解析 Header 并拆分
        original_header, column_header, parts = [], "", []
        
        print(f"[*] 正在拆分 VCF，过滤参数: Min={args.min_length}, Max={args.max_length}")
        with open(args.input, 'r') as f:
            for line in f:
                if line.startswith('##'): original_header.append(line)
                elif line.startswith('#CHROM'):
                    column_header = line
                    break
            
            sub_header = original_header + [column_header]
            data_buffer, part_idx = [], 0
            for line in f:
                data_buffer.append(line)
                if len(data_buffer) >= args.lines:
                    p_path = os.path.join(work_dir, f"part_{part_idx:04d}.vcf")
                    with open(p_path, 'w') as out:
                        out.writelines(sub_header + data_buffer)
                    parts.append(p_path); data_buffer = []; part_idx += 1
            if data_buffer:
                p_path = os.path.join(work_dir, f"part_{part_idx:04d}.vcf")
                with open(p_path, 'w') as out:
                    out.writelines(sub_header + data_buffer)
                parts.append(p_path)

        # 2. 并行调度
        print(f"[*] 任务数: {len(parts)}, 并行进程: {args.jobs}")
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            results = list(executor.map(run_repmask, parts, [args]*len(parts)))

        # 3. 合并
        print(f"[*] 正在合并至: {args.output}")
        valid_res = sorted([r for r in results if r is not None])
        with open(args.output, 'w') as final_f:
            final_f.writelines(original_header)
            for extra in RM_HEADER_EXTRAS: final_f.write(extra + "\n")
            final_f.write(column_header)
            for res in valid_res:
                with open(res, 'r') as rf:
                    for line in rf:
                        if not line.startswith('#'): final_f.write(line)
        print(f"[√] 全部分块处理并合并完成")

    finally:
        if not args.keep_cache:
            shutil.rmtree(work_dir)
            print("[*] 临时工作目录已清理")

if __name__ == "__main__":
    main()
