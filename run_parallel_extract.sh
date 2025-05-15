#!/bin/bash

# 运行并行提取 SNP 的脚本
# 依赖：valid_regions.tsv 文件，extract_snps.sh 可执行脚本

# 参数说明：
# - --bar 显示进度条
# - -j 16 并行 16 个任务
# - --colsep '\t' 指定以 TAB 分隔列
# - :::: valid_regions.tsv 表示读取输入文件

echo "✅ 启动并行提取 SNP ..."

parallel --bar -j 16 --colsep '\t' ./extract_snps.sh {1} {2} {3} {4} :::: valid_regions.tsv

echo "✅ 所有任务已完成。"
