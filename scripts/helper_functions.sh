#!/bin/bash
#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $1" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        return 1
    fi
    return 0
}

check_tool() {
    if ! command -v "$1" &> /dev/null; then
        log_error "Tool not found: $1"
        return 1
    fi
    return 0
}

ensure_dir() {
    [[ ! -d "$1" ]] && mkdir -p "$1"
}

start_timer() {
    export STEP_START_TIME=$(date +%s)
}

end_timer() {
    local step_name="$1"
    local end_time=$(date +%s)
    local duration=$((end_time - STEP_START_TIME))
    log_info "${step_name} completed in ${duration} seconds"
    echo "${step_name},${duration}" >> "${LOG_DIR}/runtime.csv"
}

check_exit() {
    if [[ $? -ne 0 ]]; then
        log_error "$1 failed"
        exit 1
    fi
    log_info "$1 completed"
}

get_python_bin() {
    if [[ -n "${PYTHON_BIN:-}" ]]; then
        return 0
    fi
    if command -v python3 &>/dev/null; then
        PYTHON_BIN="python3"
    elif command -v python &>/dev/null; then
        PYTHON_BIN="python"
    else
        log_error "Python not found for parsing hap.py summary"
        return 1
    fi
    export PYTHON_BIN
}

normalize_vcf() {
    local input_vcf="$1"
    local output_vcf="$2"
    local ref_fasta="$3"

    bcftools norm \
        -f "${ref_fasta}" \
        -m -both \
        -Oz -o "${output_vcf}" \
        "${input_vcf}"
    tabix -f -p vcf "${output_vcf}"
}

run_happy() {
    local truth_vcf="$1"
    local query_vcf="$2"
    local out_prefix="$3"

    ensure_dir "$(dirname "${out_prefix}")"
    hap.py "${truth_vcf}" "${query_vcf}" \
        -f "${HIGH_CONF_BED}" \
        -r "${REF_FASTA}" \
        -o "${out_prefix}"
}

parse_happy_summary() {
    local summary_csv="$1"
    local caller="$2"

    get_python_bin || return 1

    "${PYTHON_BIN}" - "${summary_csv}" "${caller}" <<'PY'
import csv
import re
import sys

summary_csv, caller = sys.argv[1:3]

def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())

with open(summary_csv, newline="") as handle:
    reader = csv.DictReader(handle)
    if reader.fieldnames is None:
        sys.exit(0)
    normalized = {norm(name): name for name in reader.fieldnames}

    def pick(*names):
        for name in names:
            key = norm(name)
            if key in normalized:
                return normalized[key]
        return None

    type_key = pick("type")
    tp_key = pick("tp")
    fp_key = pick("fp")
    fn_key = pick("fn")
    precision_key = pick("precision", "prec")
    recall_key = pick("recall")
    f1_key = pick("f1", "f1score")

    if not all([type_key, tp_key, fp_key, fn_key, precision_key, recall_key, f1_key]):
        sys.exit(0)

    wanted = {
        "ALL": {"all", "total", "overall"},
        "SNP": {"snp"},
        "INDEL": {"indel"},
    }

    for row in reader:
        row_type = row.get(type_key, "").strip().lower()
        for out_type, labels in wanted.items():
            if row_type in labels:
                print(
                    f"{caller}\t{out_type}\t"
                    f"{row.get(tp_key, '')}\t{row.get(fp_key, '')}\t{row.get(fn_key, '')}\t"
                    f"{row.get(precision_key, '')}\t{row.get(recall_key, '')}\t{row.get(f1_key, '')}"
                )
                break
PY
}
