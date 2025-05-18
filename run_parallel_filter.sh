#!/bin/bash
# === 并行批量运行所有 VCF 文件 ===

VCF_DIR="/projects/YangLabData/Ruilong/APE_project/snp_by_region_vcftools_parallel" # 输入 VCF 文件路径
SCRIPT_PATH="/projects/YangLabData/Ruilong/APE_project/filter_snp.sh"

# 导出变量以便 parallel 使用
export SCRIPT_PATH

# 开始并行处理（最多并行 8 个任务，可根据节点调整 -j）
find "$VCF_DIR" -name "*.vcf" | parallel -j 16 --bar 'bash "$SCRIPT_PATH" {} 0.05 0.2'
