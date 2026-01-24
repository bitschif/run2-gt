#!/bin/bash
#===============================================================================
# STEP 03: Variant Calling - GATK HaplotypeCaller
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

CALLER="gatk"
log_info "===== STEP 03: ${CALLER} HaplotypeCaller ====="
start_timer

# Input
source "${PREPROC_DIR}/bam_path.sh"
check_file "${FINAL_BAM}" || exit 1
check_tool hap.py || exit 1

OUT_DIR="${VARIANT_DIR}/${CALLER}"

#-------------------------------------------------------------------------------
# 1. Run HaplotypeCaller
#-------------------------------------------------------------------------------
log_info "Running GATK HaplotypeCaller..."

RAW_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_raw.vcf.gz"

gatk HaplotypeCaller \
    --java-options "${JAVA_OPTS}" \
    -R "${REF_FASTA}" \
    -I "${FINAL_BAM}" \
    -O "${RAW_VCF}" \
    --standard-min-confidence-threshold-for-calling "${GATK_STAND_CALL_CONF}" \
    --native-pair-hmm-threads "${THREADS}" \
    2>&1 | tee "${LOG_DIR}/${CALLER}.log"

check_exit "HaplotypeCaller"

#-------------------------------------------------------------------------------
# 2. Hard filtering (GATK Best Practices)
#-------------------------------------------------------------------------------
log_info "Applying hard filters..."

# SNPs
gatk SelectVariants -V "${RAW_VCF}" -select-type SNP -O "${OUT_DIR}/snps_raw.vcf.gz"
gatk VariantFiltration \
    -V "${OUT_DIR}/snps_raw.vcf.gz" \
    -O "${OUT_DIR}/snps_filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
    --filter-name "SNP_FILTER"

# INDELs
gatk SelectVariants -V "${RAW_VCF}" -select-type INDEL -O "${OUT_DIR}/indels_raw.vcf.gz"
gatk VariantFiltration \
    -V "${OUT_DIR}/indels_raw.vcf.gz" \
    -O "${OUT_DIR}/indels_filtered.vcf.gz" \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0" \
    --filter-name "INDEL_FILTER"

# Merge
gatk MergeVcfs \
    -I "${OUT_DIR}/snps_filtered.vcf.gz" \
    -I "${OUT_DIR}/indels_filtered.vcf.gz" \
    -O "${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz"

#-------------------------------------------------------------------------------
# 3. Extract PASS variants
#-------------------------------------------------------------------------------
log_info "Extracting PASS variants..."

PASS_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.vcf.gz"
bcftools view -f PASS "${OUT_DIR}/${PREFIX}_${CALLER}_filtered.vcf.gz" -Oz -o "${PASS_VCF}"
tabix -p vcf "${PASS_VCF}"

# Normalize for benchmarking consistency
NORMALIZED_VCF="${OUT_DIR}/${PREFIX}_${CALLER}_pass.norm.vcf.gz"
normalize_vcf "${PASS_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"

# Split by type for benchmarking
bcftools view -v snps "${NORMALIZED_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
bcftools view -v indels "${NORMALIZED_VCF}" -Oz -o "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"
tabix -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_snp.vcf.gz"
tabix -p vcf "${OUT_DIR}/${PREFIX}_${CALLER}_indel.vcf.gz"

#-------------------------------------------------------------------------------
# 4. Stats
#-------------------------------------------------------------------------------
bcftools stats "${NORMALIZED_VCF}" > "${OUT_DIR}/${PREFIX}_${CALLER}_stats.txt"

N_SNP=$(bcftools view -H -v snps "${NORMALIZED_VCF}" | wc -l)
N_INDEL=$(bcftools view -H -v indels "${NORMALIZED_VCF}" | wc -l)

log_info "Results:  $((N_SNP + N_INDEL)) variants (${N_SNP} SNPs, ${N_INDEL} INDELs)"

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

end_timer "03_${CALLER}"
log_info "===== ${CALLER} Complete ====="
