#!/bin/bash

##############################
# Author: Zijia Ni
# Date: 2023/1/10
##############################

####################################
### Params 
####################################
# Global params
WorkDir="$(pwd)/"
CacheDir="/pdata1/home/nizijia/project/02.fanya/A.1.BJ_FY_graph_rna/tmp"

# 默认值
FinalTxt=""
Suffix=""
Fq1Suffix=""
Fq2Suffix=""

# 显示帮助信息
show_help() {
    echo "USAGE: $0 [OPTIONS]"
    echo "Need in target workdir"
    echo "OPTIONS:"
    echo "  -f, --final=<文件名>    设置最终文件名"
    echo "  -s, --suffix=<后缀>     设置后缀"
    echo "  -1, --fq1-suffix=<后缀> 设置Fq1的后缀"
    echo "  -2, --fq2-suffix=<后缀> 设置Fq2的后缀"
    echo "  -h, --help              显示此帮助信息"
}

# 处理命名参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -f|--final)
            FinalTxt="$2"
            shift 2
            ;;
        -s|--suffix)
            Suffix="$2"
            shift 2
            ;;
        -1|--fq1-suffix)
            Fq1Suffix="$2"
            shift 2
            ;;
        -2|--fq2-suffix)
            Fq2Suffix="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "错误: 无效选项: $1" >&2
            show_help
            exit 1
            ;;
    esac
done

# 输出参数值
echo "FinalTxt: $FinalTxt"
echo "Suffix: $Suffix"
echo "Fq1Suffix: $Fq1Suffix"
echo "Fq2Suffix: $Fq2Suffix"

######################################################
### Running script
###################################################### 
echo "$WorkDir"

ls *${Suffix} |cut -d "." -f1|cut -d "_" -f 1|sort|uniq > ${CacheDir}/cache1
ls *${Fq1Suffix}.${Suffix} > ${CacheDir}/cache2
ls *${Fq2Suffix}.${Suffix} > ${CacheDir}/cache3

paste -d "\t" `ls ${CacheDir}/cache?` > ${CacheDir}/cache4
cat ${CacheDir}/cache4 | awk -v dir=$WorkDir '{print $1"\t"dir$2"\t"dir$3}' > $FinalTxt

# 存一份样本文件在上级目录
# cp $FinalTxt ../$FinalTxt
# rm $FinalTxt

# remove cache file
rm ${CacheDir}/cache?

