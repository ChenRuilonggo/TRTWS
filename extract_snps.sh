#!/bin/bash

# 参数：chrom, start, end, region_id
chrom=$1
start=$2
end=$3
region_id=$4

VCF="/projects/YangLabData/rparrish/ROSMAP_WGS/genotype/CHR${chrom}_ROSMAP_WGS_b38_lifted_TIGAR.vcf.gz"
SAMPLES="sample_list.txt"
OUTDIR="snp_by_region_vcftools_parallel"
mkdir -p "$OUTDIR"

if [[ ! -f "$VCF" ]]; then
    echo " VCF not found for chr${chrom}: $VCF"
    exit 1
fi

tmp_vcf="${OUTDIR}/${region_id}.raw.vcf"
out_vcf="${OUTDIR}/${region_id}.vcf"

# ✅ 跳过已完成
if [[ -f "$out_vcf" ]]; then
    echo " Skipping $region_id (already exists)"
    exit 0
fi

start_time=$(date +%s)

echo " Extracting ${start}-${end} → $region_id"

vcftools \
  --gzvcf "$VCF" \
  --chr "$chrom" \
  --from-bp "$start" \
  --to-bp "$end" \
  --keep "$SAMPLES" \
  --recode \
  --stdout > "$tmp_vcf"

awk 'BEGIN{OFS="\t"} /^#/ {print; next} length($4)==1 && length($5)==1 && $5 != "*" {print}' "$tmp_vcf" > "$out_vcf"

rm "$tmp_vcf"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo " $region_id done in ${elapsed}s"
echo -e "$region_id\t$elapsed" >> "${OUTDIR}/region_timings.tsv"
