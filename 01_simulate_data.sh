#!/bin/bash
#===============================================================================
# STEP 01: Simulate Data
# Tools: simutator (iqbal-lab-org/simutator), art_illumina
#===============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config/config.sh"
source "${SCRIPT_DIR}/scripts/helper_functions.sh"

log_info "===== STEP 01: Simulate Data ====="
start_timer

#-------------------------------------------------------------------------------
# 1. Check tools and reference
#-------------------------------------------------------------------------------
check_tool simutator || exit 1
check_tool art_illumina || exit 1
check_file "${REF_FASTA}" || exit 1

#-------------------------------------------------------------------------------
# 2. Generate mutations with simutator
#-------------------------------------------------------------------------------
log_info "Generating mutations with simutator..."
log_info "  SNPs: every ${SNP_DIST} bp"
log_info "  Deletions: ${DEL_LEN} bp every ${DEL_DIST} bp"
log_info "  Insertions: ${INS_LEN} bp every ${INS_DIST} bp"

SIM_PREFIX="${SIM_DIR}/${PREFIX}"

simutator mutate_fasta \
    --snps ${SNP_DIST} \
    --dels ${DEL_DIST}:${DEL_LEN} \
    --ins ${INS_DIST}:${INS_LEN} \
    --seed ${SEED} \
    "${REF_FASTA}" \
    "${SIM_PREFIX}"

check_exit "simutator"

#-------------------------------------------------------------------------------
# 3. Find output files from simutator
#-------------------------------------------------------------------------------
log_info "Processing simutator output..."

MUTATED_FASTA=$(ls ${SIM_PREFIX}*.fa 2>/dev/null | grep -v ".original" | head -1)
ORIGINAL_VCF=$(ls ${SIM_PREFIX}*.original.vcf 2>/dev/null | head -1)

if [[ -z "${MUTATED_FASTA}" ]] || [[ -z "${ORIGINAL_VCF}" ]]; then
    log_error "simutator output files not found"
    ls -la "${SIM_DIR}/"
    exit 1
fi

log_info "  Mutated FASTA: ${MUTATED_FASTA}"
log_info "  Truth VCF: ${ORIGINAL_VCF}"

samtools faidx "${MUTATED_FASTA}"

#-------------------------------------------------------------------------------
# 4. Process truth VCF - Fix missing header before bcftools
#-------------------------------------------------------------------------------
log_info "Processing truth VCF..."

FIXED_VCF="${SIM_DIR}/${PREFIX}_truth_fixed.vcf"
HEADER_FILE="${SIM_DIR}/header_fix.txt"

# Create new header with FORMAT GT definition
cat > "${HEADER_FILE}" << 'EOF'
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
EOF

# Build fixed VCF: new header + contig line + column header + data
{
    cat "${HEADER_FILE}"
    grep "^##contig" "${ORIGINAL_VCF}" || echo "##contig=<ID=chr22>"
    grep "^#CHROM" "${ORIGINAL_VCF}"
    grep -v "^#" "${ORIGINAL_VCF}" | sort -k1,1 -k2,2n
} > "${FIXED_VCF}"

# Compress and index using bgzip (not bcftools sort to avoid parsing issues)
bgzip -c "${FIXED_VCF}" > "${TRUTH_VCF}"
tabix -p vcf "${TRUTH_VCF}"

# Create separate SNP and INDEL files
bcftools view -v snps "${TRUTH_VCF}" -Oz -o "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
bcftools view -v indels "${TRUTH_VCF}" -Oz -o "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz"
tabix -p vcf "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz"

# Count variants
TOTAL_VARS=$(zcat "${TRUTH_VCF}" | grep -v "^#" | wc -l)
SNP_COUNT=$(zcat "${SIM_DIR}/${PREFIX}_truth_snp.vcf.gz" | grep -v "^#" | wc -l)
INDEL_COUNT=$(zcat "${SIM_DIR}/${PREFIX}_truth_indel.vcf.gz" | grep -v "^#" | wc -l)

log_info "  Total variants: ${TOTAL_VARS}"
log_info "  SNPs: ${SNP_COUNT}"
log_info "  Indels: ${INDEL_COUNT}"

# Cleanup temp files
rm -f "${FIXED_VCF}" "${HEADER_FILE}"

#-------------------------------------------------------------------------------
# 5. Generate reads with ART
#-------------------------------------------------------------------------------
log_info "Generating reads with ART Illumina..."
log_info "  Coverage: ${COVERAGE}x"
log_info "  Read length: ${READ_LENGTH} bp"
log_info "  Fragment: ${FRAGMENT_MEAN} +/- ${FRAGMENT_SD} bp"

READ_PREFIX="${SIM_DIR}/${PREFIX}"

art_illumina \
    -ss ${ART_PLATFORM} \
    -i "${MUTATED_FASTA}" \
    -p \
    -l ${READ_LENGTH} \
    -f ${COVERAGE} \
    -m ${FRAGMENT_MEAN} \
    -s ${FRAGMENT_SD} \
    -rs ${SEED} \
    -o "${READ_PREFIX}_" \
    -na

check_exit "art_illumina"

# Rename and compress output files
mv "${READ_PREFIX}_1.fq" "${READ_PREFIX}_R1.fastq"
mv "${READ_PREFIX}_2.fq" "${READ_PREFIX}_R2.fastq"
gzip -f "${READ_PREFIX}_R1.fastq"
gzip -f "${READ_PREFIX}_R2.fastq"

log_info "  R1: ${READ_PREFIX}_R1.fastq.gz"
log_info "  R2: ${READ_PREFIX}_R2.fastq.gz"

#-------------------------------------------------------------------------------
# 6. Create callable regions BED
#-------------------------------------------------------------------------------
log_info "Creating callable regions BED..."

awk -v OFS='\t' '{print $1, 0, $2}' "${REF_FAI}" > "${HIGH_CONF_BED}"

log_info "  BED: ${HIGH_CONF_BED}"

end_timer "01_simulate_data"
log_info "===== Simulation Complete ====="
