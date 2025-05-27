#!/bin/bash
# === Ruilong Chen's SNP Preprocessing Script with Duplicate Removal ===

# ===== 参数 =====
VCF_FILE=$1                                     # 输入 VCF 文件路径
MAF_THRESHOLD=${2:-0.01}                        # 默认 MAF 阈值为 0.01
LD_R2_THRESHOLD=${3:-0.95}                      # 默认 LD r² 阈值为 0.95

# ===== 自动生成输出文件夹和前缀 =====
VCF_BASENAME=$(basename "$VCF_FILE" .vcf)
OUT_DIR="/projects/YangLabData/Ruilong/APE_project/plink_result/$VCF_BASENAME"
OUT_PREFIX="$VCF_BASENAME"
mkdir -p "$OUT_DIR"

echo "Output directory: $OUT_DIR"
echo "Output file prefix: $OUT_PREFIX"

# ===== Step 0: 去除重复变异（相同 CHROM+POS+REF+ALT）=====
echo "Step 0: Remove duplicate variants from VCF by CHROM+POS+REF+ALT"
DEDUP_VCF="${OUT_DIR}/${OUT_PREFIX}_dedup.vcf"
awk -F'\t' '
BEGIN { OFS = FS }
# 保留 VCF header
/^#/ { print; next }

{
  # 字段CHROM=$1, POS=$2, REF=$4, ALT=$5
  ref = $4
  alt = $5

  # 创建唯一标识CHROM_POS_min(REF,ALT)_max(REF,ALT)
  key = $1 "_" $2 "_" (ref < alt ? ref : alt) "_" (ref > alt ? ref : alt)

  # 去除完全重复或对称等位重复的变异
  if (!(key in seen)) {
    seen[key] = 1
    print
  }
}
' "$VCF_FILE" > "$DEDUP_VCF"

echo " VCF de-duplication completed: $DEDUP_VCF"

# ===== Step 1: VCF → PLINK =====
echo "Step 1: Convert VCF to PLINK binary format (initial conversion)"
plink --vcf "$DEDUP_VCF" \
      --vcf-half-call missing \
      --set-missing-var-ids @:#:\$1:\$2 \
      --make-bed \
      --out "${OUT_DIR}/${OUT_PREFIX}_step1_deduped"

if [ $? -ne 0 ]; then
    echo "❌ Error during PLINK Step 1 (initial conversion). Exiting."
    exit 1
fi
echo "✅ PLINK Step 1 completed."

# ===== Step 2: MAF 过滤 =====
echo "Step 2: Filter by MAF >= $MAF_THRESHOLD"
plink --bfile "${OUT_DIR}/${OUT_PREFIX}_step1_deduped" \
      --maf "$MAF_THRESHOLD" \
      --make-bed \
      --out "${OUT_DIR}/${OUT_PREFIX}_maf_filtered"

if [ $? -ne 0 ]; then
    echo "❌ Error during PLINK Step 2 (MAF filtering). Exiting."
    exit 1
fi
echo "✅ PLINK Step 2 completed."

# ===== Step 3: LD Pruning =====
echo "Step 3: LD pruning with r² <= $LD_R2_THRESHOLD"
plink --bfile "${OUT_DIR}/${OUT_PREFIX}_maf_filtered" \
      --indep-pairwise 50 5 "$LD_R2_THRESHOLD" \
      --out "${OUT_DIR}/${OUT_PREFIX}_pruned"

if [ $? -ne 0 ]; then
    echo "❌ Error during PLINK Step 3 (LD pruning). Exiting."
    exit 1
fi
echo "✅ PLINK Step 3 completed."

# ===== Step 4: 提取 Pruned SNPs =====
echo "Step 4: Extract pruned SNPs"
plink --bfile "${OUT_DIR}/${OUT_PREFIX}_maf_filtered" \
      --extract "${OUT_DIR}/${OUT_PREFIX}_pruned.prune.in" \
      --make-bed \
      --out "${OUT_DIR}/${OUT_PREFIX}_final"

if [ $? -ne 0 ]; then
    echo "❌ Error during PLINK Step 4 (extract pruned SNPs). Exiting."
    exit 1
fi
echo "✅ PLINK Step 4 completed."



# ===== Step 5: Convert final binary PLINK files to .raw genotype matrix =====
echo "Step 5: Generate .raw genotype matrix (0/1/2 format)"
plink --bfile "${OUT_DIR}/${OUT_PREFIX}_final" \
      --recode A \
      --out "${OUT_DIR}/${OUT_PREFIX}_final"

if [ $? -ne 0 ]; then
    echo "❌ Error during PLINK Step 5 (recode to .raw). Exiting."
    exit 1
fi
echo "✅ PLINK Step 5 (.raw file generation) completed."


# ===== Step 6: 清理中间文件（只保留 final 输出） =====
echo "Step 6: Clean up intermediate files..."

cd "$OUT_DIR"

# 删除中间 VCF
rm -f "${OUT_PREFIX}_dedup.vcf"

# 删除 step1、maf_filtered、pruned 中间产物
rm -f "${OUT_PREFIX}_step1_deduped".*
rm -f "${OUT_PREFIX}_maf_filtered".*
rm -f "${OUT_PREFIX}_pruned".*

echo "✅ Cleanup complete. Only final PLINK files retained:"
ls -lh "${OUT_PREFIX}_final".*
