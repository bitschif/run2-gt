#!/bin/bash
#===============================================================================
# STEP 06: Variant Calling - FreeBayes
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

CALLER="freebayes"
log_info "===== STEP 06: ${CALLER} ====="
start_timer

# Input
source "${PREPROC_DIR}/bam_path.sh"
check_file "${FINAL_BAM}" || exit 1
check_tool freebayes || exit 1
check_tool hap.py || exit 1

OUT_DIR="${VARIANT_DIR}/${CALLER}"

#-------------------------------------------------------------------------------
# 1. Run FreeBayes
#-------------------------------------------------------------------------------
log_info "Running FreeBayes..."

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf"

freebayes \
    -f "${REF_FASTA}" \
    -b "${FINAL_BAM}" \
    --min-alternate-count "${FB_MIN_ALT_COUNT}" \
    --min-alternate-fraction "${FB_MIN_ALT_FRACTION}" \
    --min-mapping-quality "${MIN_MAPPING_QUALITY}" \
    --min-base-quality "${MIN_BASE_QUALITY}" \
    --genotype-qualities \
    > "${RAW_VCF}" \
    2> "${LOG_DIR}/${CALLER}.log"

check_exit "FreeBayes"

# Compress and index
bgzip -f "${RAW_VCF}"
tabix -f -p vcf "${RAW_VCF}.gz"

#-------------------------------------------------------------------------------
# 2. Filter and normalize
#-------------------------------------------------------------------------------
log_info "Filtering and normalizing variants..."

FILTERED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz"

# Quality filter
bcftools filter \
    -i 'QUAL>30 && INFO/DP>10' \
    "${RAW_VCF}.gz" | \
bcftools norm \
    -f "${REF_FASTA}" \
    -m -both \
    -Oz -o "${FILTERED_VCF}"

tabix -f -p vcf "${FILTERED_VCF}"

#-------------------------------------------------------------------------------
# 3. Extract high-quality variants
#-------------------------------------------------------------------------------
log_info "Extracting high-quality variants..."

PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"

# Add PASS tag to high-quality variants
bcftools filter \
    -i 'QUAL>30 && INFO/DP>10' \
    -s LowQual \
    "${FILTERED_VCF}" | \
bcftools view -f "PASS,." -Oz -o "${PASS_VCF}"

tabix -f -p vcf "${PASS_VCF}"

# Normalize for benchmarking consistency
NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
normalize_vcf "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

# Split by type
bcftools view -v snps "${NORMALIZED_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${NORMALIZED_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -f -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 4. Stats
#-------------------------------------------------------------------------------
bcftools stats "${NORMALIZED_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${NORMALIZED_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${NORMALIZED_VCF}" | wc -l)

log_info "Results: $((N_SNP + N_INDEL)) variants (${N_SNP} SNPs, ${N_INDEL} INDELs)"

#-------------------------------------------------------------------------------
# 5. Benchmark with hap.py
#-------------------------------------------------------------------------------
log_info "Benchmarking ${CALLER} with hap.py..."
HAPPY_PREFIX="${OUT_DIR}/benchmark/${PREFIX}_${CALLER}"
run_happy "${TRUTH_VCF}" "${NORMALIZED_VCF}" "${HAPPY_PREFIX}"

SUMMARY_CSV="${HAPPY_PREFIX}.summary.csv"
SUMMARY_TSV="${OUT_DIR}/${PREFIX}_${CALLER}_happy_summary.tsv"
echo -e "Caller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1" > "${SUMMARY_TSV}"
if [[ -f "${SUMMARY_CSV}" ]]; then
    parse_happy_summary "${SUMMARY_CSV}" "${CALLER}" >> "${SUMMARY_TSV}"
else
    log_warn "hap.py summary not found for ${CALLER}"
fi

end_timer "06_${CALLER}"
log_info "===== ${CALLER} Complete ====="
