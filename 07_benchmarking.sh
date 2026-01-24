#!/bin/bash
#===============================================================================
# STEP 07: Benchmarking
# Compare variant calls against truth set
# Tools: hap.py (Illumina)
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 07: Benchmarking ====="
start_timer

check_tool "hap.py" || exit 1

CALLERS=("gatk" "deepvariant" "strelka2" "freebayes")
TRUTH_SNP="${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
TRUTH_INDEL="${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"

check_file "${TRUTH_VCF}" || exit 1
check_file "${HIGH_CONF_BED}" || exit 1
check_file "${REF_FASTA}" || exit 1

#-------------------------------------------------------------------------------
# Main benchmarking
#-------------------------------------------------------------------------------
SUMMARY="${BENCH_DIR}/benchmark_summary.tsv"
echo -e "Caller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1" > "${SUMMARY}"

for caller in "${CALLERS[@]}"; do
    log_info "Benchmarking ${caller}..."
    
    QUERY_VCF="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_pass.vcf.gz"
    NORMALIZED_VCF="${VARIANT_DIR}/${caller}/${PREFIX}_${caller}_pass.norm.vcf.gz"
    if [[ !  -f "${QUERY_VCF}" ]]; then
        log_warn "  VCF not found, skipping..."
        continue
    fi
    
    BENCH_CALLER="${BENCH_DIR}/${caller}"
    ensure_dir "${BENCH_CALLER}"
    
    log_info "  hap.py..."
    HAPPY_PREFIX="${BENCH_CALLER}/happy/${PREFIX}_${caller}"
    if [[ -f "${NORMALIZED_VCF}" ]]; then
        run_happy "${TRUTH_VCF}" "${NORMALIZED_VCF}" "${HAPPY_PREFIX}"
    else
        normalize_vcf "${QUERY_VCF}" "${NORMALIZED_VCF}" "${REF_FASTA}"
        run_happy "${TRUTH_VCF}" "${NORMALIZED_VCF}" "${HAPPY_PREFIX}"
    fi

    SUMMARY_CSV="${HAPPY_PREFIX}.summary.csv"
    if [[ -f "${SUMMARY_CSV}" ]]; then
        parse_happy_summary "${SUMMARY_CSV}" "${caller}" >> "${SUMMARY}"
    else
        log_warn "  hap.py summary not found for ${caller}"
    fi
done

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
log_info "===== Benchmarking Summary ====="
if command -v column &>/dev/null; then
    column -t -s$'\t' "${SUMMARY}"
else
    cat "${SUMMARY}"
fi

end_timer "07_benchmarking"
log_info "===== Benchmarking Complete ====="
