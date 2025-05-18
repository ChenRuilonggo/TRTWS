#!/bin/bash

# ========== CONFIG ==========

VCF_DIR="/projects/YangLabData/Ruilong/APE_project/snp_by_region_vcftools_parallel"
CACHE_DIR="/projects/YangLabData/Ruilong/APE_project/tools/vep_data/cache"
# FASTA path confirmed to be without .gz
FASTA="$CACHE_DIR/homo_sapiens/110_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"

OUT_TXT_DIR="$VCF_DIR/rsid_txt"
mkdir -p "$OUT_TXT_DIR"

# ========== ACTIVATE CONDA ENVIRONMENT ==========
# Ensure this is the correct name of your VEP Conda environment for VEP 110
CONDA_ENV_NAME="vep_v110_env" # This should be the environment where VEP 110 was successfully installed.

echo "Activating Conda environment: $CONDA_ENV_NAME"
if [ -f "/projects/YangLabData/Ruilong/APE_project/tools/miniconda3/etc/profile.d/conda.sh" ]; then
    source "/projects/YangLabData/Ruilong/APE_project/tools/miniconda3/etc/profile.d/conda.sh"
else
    echo "ERROR: Conda activation script not found. Please check Miniconda installation path."
    exit 1
fi

conda activate "$CONDA_ENV_NAME"

if [ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV_NAME" ]; then
    echo "ERROR: Failed to activate Conda environment: $CONDA_ENV_NAME. Exiting."
    exit 1
fi
echo "Conda environment $CONDA_DEFAULT_ENV activated."

# ========== SET VEP PATH (now uses Conda's VEP) ==========
VEP_PATH="vep"

# ========== PARALLEL FUNCTION ==========
process_one_vcf() {
  vcf_file="$1"  # Corrected variable assignment
  filename=$(basename "$vcf_file" .vcf)  # Corrected variable assignment
  out_txt="$OUT_TXT_DIR/${filename}.txt"  # Corrected variable assignment

  echo "Annotating: $filename.vcf"

  # Use the VEP_PATH variable, which now refers to Conda's VEP.
  "$VEP_PATH" \
    -i "$vcf_file" \
    -o STDOUT \
    --vcf \
    --offline --cache \
    --assembly GRCh38 \
    --species homo_sapiens \
    --dir_cache "$CACHE_DIR" \
    --fasta "$FASTA" \
    --symbol \
    --force_overwrite \
  | grep -v "^##" \
  | awk -F'\t' 'BEGIN{OFS="\t"} {
      # $1 is CHROM, $2 is POS, $8 is INFO field
      info_field = $8;  # Corrected variable reference
      split(info_field, all_info_tags, ";");  # Corrected variable reference
      for (tag_idx in all_info_tags) {
          if (all_info_tags[tag_idx] ~ /^CSQ=/) {  # Corrected variable reference
              csq_full_string = all_info_tags[tag_idx];  # Corrected variable reference
              sub(/^CSQ=/, "", csq_full_string);  # Corrected variable reference
              split(csq_full_string, csq_array, ",");  # Corrected variable reference
              for (csq_entry_idx in csq_array) {
                  current_csq_entry = csq_array[csq_entry_idx];  # Corrected variable reference
                  split(current_csq_entry, csq_fields, "|");  # Corrected variable reference
                  # Existing_variation is the 18th field in your CSQ format
                  if (length(csq_fields) >= 18 && csq_fields[18] ~ /^rs[0-9]+/) {  # Corrected regex
                      print $1, $2, csq_fields[18];  # Output CHROM, POS, rsID
                      next; 
                  }
              }
              break; 
          }
      }
  }' | sort -u > "$out_txt"  # Corrected variable reference

  # Check if awk produced any output (i.e., if any rsIDs were found) before adding the header.
  if [ -s "$out_txt" ]; then
      (echo -e "CHROM\tPOS\trsID"; cat "$out_txt") > "${out_txt}.tmp" && mv "${out_txt}.tmp" "$out_txt"
  else
      echo -e "CHROM\tPOS\trsID" > "$out_txt"  # If no rsIDs found, create file with only header
      echo "No rsIDs extracted by awk for $filename.vcf (output file created with header only)."
  fi

  echo "Done: $out_txt"
}
# Export the function and necessary variables so 'parallel' can access them.
export -f process_one_vcf
export VEP_PATH CACHE_DIR FASTA OUT_TXT_DIR

# ========== RUN IN PARALLEL ==========

# Find all .vcf files in VCF_DIR and process them using GNU Parallel with 24 jobs.
# find "$VCF_DIR" -name "*.vcf" | parallel -j 24 process_one_vcf {}

process_one_vcf "/projects/YangLabData/Ruilong/APE_project/snp_by_region_vcftools_parallel/18_54027805_57546606.vcf"

echo "All VCFs processed. rsID files saved in: $OUT_TXT_DIR"

# ========== DEACTIVATE CONDA ENVIRONMENT (good practice) ==========
conda deactivate
echo "Conda environment $CONDA_DEFAULT_ENV (was $CONDA_ENV_NAME) deactivated."