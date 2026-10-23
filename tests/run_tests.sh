#!/bin/bash
#
# Automated tests of the CUT&RUN pipeline (about 2 minutes).
#
# Usage: run_tests.sh [test_dir]
# Software: the conda environments of the pipeline, given by name or path:
#   TEST_CONDA_ENV=cutrun TEST_CONDA_ENV_R=cutrun_r bash tests/run_tests.sh
# Leave them unset if all tools and R packages are already on PATH.
#
# 1. Syntax of all bash and R scripts
# 2. Helper functions: sample list and peak settings
# 3. build_reference.sh on a small simulated genome
# 4. Input checks (validate.sh): good input passes, mistakes are caught
# 5. Whole pipeline on simulated reads (make_test_data.sh), FASTQ to DESeq2,
#    compared with the known simulated peaks and differences
# 6. Reruns: finished steps are skipped; the diff part runs on its own

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly TESTS_DIR
PIPELINE_DIR="$(dirname "${TESTS_DIR}")"
readonly PIPELINE_DIR

export CONDA_ENV="${TEST_CONDA_ENV:-}"
export CONDA_ENV_R="${TEST_CONDA_ENV_R:-}"

n_pass=0
n_fail=0
test_dir=""
config=""
peak_callers="macs2"

#######################################
# Records a passed or failed test.
# Globals:
#   n_pass, n_fail
# Arguments:
#   Test description.
#######################################
pass() {
  echo "  PASS  $1"
  n_pass=$(( n_pass + 1 ))
}
fail() {
  echo "  FAIL  $1"
  n_fail=$(( n_fail + 1 ))
}

#######################################
# Passes if a command succeeds (its output is hidden).
# Arguments:
#   Test description, command and its arguments.
#######################################
check() {
  if "${@:2}" > /dev/null 2>&1; then pass "$1"; else fail "$1"; fi
}

#######################################
# Passes if a value is the expected value.
# Arguments:
#   Test description, value, expected value.
#######################################
same() {
  if [[ "$2" == "$3" ]]; then
    pass "$1"
  else
    fail "$1  (got '$2', expected '$3')"
  fi
}

#######################################
# Passes if a number is in a range.
# Arguments:
#   Test description, value, minimum, maximum.
#######################################
between() {
  if [[ "$2" =~ ^[0-9]+$ ]] && (( $2 >= $3 && $2 <= $4 )); then
    pass "$1 ($2)"
  else
    fail "$1  (got '$2', expected $3 to $4)"
  fi
}

#######################################
# Runs a helper function of scripts/common.sh with a config.
# Arguments:
#   Config, function, its arguments.
#######################################
call() {
  bash -c 'source "$0"; read_config "$1"; shift; "$@"' \
    "${PIPELINE_DIR}/scripts/common.sh" "$@" 2> /dev/null
}

#######################################
# Runs a command in the software environment of the pipeline.
# Globals:
#   config
# Arguments:
#   "tools" or "R", command, its arguments.
#######################################
run_in() {
  bash -c 'source "$0"; read_config "$1"; activate_env "$2"; shift 2; "$@"' \
    "${PIPELINE_DIR}/scripts/common.sh" "${config}" "$@"
}

#######################################
# Number of lines of a file (0 if it does not exist).
# Arguments:
#   File.
#######################################
count_lines() {
  if [[ -f "$1" ]]; then wc -l < "$1" | tr -d ' '; else echo 0; fi
}

#######################################
# Prints one column of a row of a tab-separated table with a header.
# Arguments:
#   File, value of column 1 or 2 that selects the row, column name.
#######################################
table_value() {
  awk -F'\t' -v row="$2" -v column="$3" '
      NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
      $1 == row || $2 == row { print $col[column] }' "$1" 2> /dev/null
}

#######################################
# Writes the test peak settings: FlexPeak without its tuning sweep (fast).
# Globals:
#   test_dir
#######################################
make_peak_settings() {
  awk 'BEGIN { FS = OFS = "\t" }
       !/^#/ && $1 != "mark" { $3 = $3 " --no-tune" } { print }' \
    "${PIPELINE_DIR}/config/peak_settings.txt" > "${test_dir}/peak_settings.txt"
}

#######################################
# Writes a test config that uses a given sample list.
# Globals:
#   test_dir, peak_callers, CONDA_ENV, CONDA_ENV_R
# Arguments:
#   Sample list, output config.
#######################################
make_config() {
  local comparisons="${test_dir}/comparison_list.txt"
  sed -e "s#^WORKDIR=.*#WORKDIR=\"${test_dir}/Data\"#" \
    -e "s#^PEAK_SETTINGS=.*#PEAK_SETTINGS=\"${test_dir}/peak_settings.txt\"#" \
    -e "s#^PEAK_CALLERS=.*#PEAK_CALLERS=\"${peak_callers}\"#" \
    -e "s#^DATADIR=.*#DATADIR=\"${test_dir}/Raw\"#" \
    -e "s#^SAMPLE_LIST=.*#SAMPLE_LIST=\"$1\"#" \
    -e "s#^COMPARISON_LIST=.*#COMPARISON_LIST=\"${comparisons}\"#" \
    -e "s#^LIBRARIES=.*#LIBRARIES=\"${test_dir}/indexes\"#" \
    -e "s#^REFSIZE=.*#REFSIZE=\"1e6\"#" \
    -e "s#^CONDA_ENV=.*#CONDA_ENV=\"${CONDA_ENV}\"#" \
    -e "s#^CONDA_ENV_R=.*#CONDA_ENV_R=\"${CONDA_ENV_R}\"#" \
    -e "s#^SCHEDULER=.*#SCHEDULER=\"local\"#" \
    -e "s#^THREADS=.*#THREADS=2#" \
    -e 's#^DESIGN=.*#DESIGN="~ batch + condition"#' \
    "${PIPELINE_DIR}/config/config.sh" > "$2"
}

#######################################
# Passes if validate.sh rejects a sample list with a mistake.
# Globals:
#   test_dir
# Arguments:
#   Mistake, expected error message, extra sample list row (or "header").
#######################################
check_catches() {
  local bad_list="${test_dir}/sample_list_bad.txt"
  local bad_config="${test_dir}/config_bad.sh"
  if [[ "$3" == header ]]; then
    tr '\t' ',' < "${test_dir}/sample_list.txt" > "${bad_list}"
  else
    { cat "${test_dir}/sample_list.txt"; printf '%b\n' "$3"; } > "${bad_list}"
  fi
  make_config "${bad_list}" "${bad_config}"

  local output status
  output="$(bash "${PIPELINE_DIR}/scripts/validate.sh" "${bad_config}" 2>&1)"
  status=$?
  if (( status != 0 )) && [[ "${output}" == *"$2"* ]]; then
    pass "catches: $1"
  else
    fail "catches: $1  (${output##*ERROR: })"
  fi
}

#######################################
# Checks the HTML report: data block is valid JSON with all peaks and
# comparisons, a gene can be found, and the JavaScript has no syntax error.
# Arguments:
#   Report file, number of union peaks.
#######################################
check_report() {
  check "HTML report" test -s "$1"
  same "report: peaks and comparisons" "$(python3 -c '
import json, re, sys
html = open(sys.argv[1]).read()
block = re.search(r"<script id=\"report-data\"[^>]*>(.*?)</script>", html, re.S)
mark = json.loads(block.group(1))["marks"]["H3K27ac"]
print(len(mark["peaks"]["id"]), len(mark["comparisons"]))' "$1" \
    2> /dev/null)" "$2 3"
  check "report: gene names to search" grep -q '"Gene1' "$1"
  if command -v node > /dev/null; then
    local report_js="${test_dir}/report.js"
    sed -n '/^<script>$/,/^<\/script>$/p' "$1" | sed '1d;$d' > "${report_js}"
    check "report: JavaScript syntax" node --check "${report_js}"
  fi
}

test_syntax() {
  echo "1. Syntax"
  local script
  for script in "${PIPELINE_DIR}"/*.sh "${PIPELINE_DIR}"/scripts/*.sh \
    "${PIPELINE_DIR}"/reference/*.sh "${PIPELINE_DIR}"/tests/*.sh; do
    check "bash syntax: ${script#"${PIPELINE_DIR}"/}" bash -n "${script}"
  done
  for script in "${PIPELINE_DIR}"/scripts/R/*.R; do
    check "R syntax: ${script#"${PIPELINE_DIR}"/}" \
      run_in R Rscript -e "invisible(parse('${script}'))"
  done
}

test_helpers() {
  echo "2. Helper functions"
  same "all_samples" "$(call "${config}" all_samples | wc -l | tr -d ' ')" 8
  same "sample_info: condition" \
    "$(call "${config}" sample_info Arm_H3K27ac_R1 condition)" "Arm"
  same "sample_info: two lanes" \
    "$(call "${config}" sample_info Naive_H3K27ac_R1 fastq1)" \
    "Naive_H3K27ac_R1_L1_R1.fastq.gz,Naive_H3K27ac_R1_L2_R1.fastq.gz"
  same "samples_of mark and condition" \
    "$(call "${config}" samples_of H3K27ac Naive | xargs)" \
    "Naive_H3K27ac_R1 Naive_H3K27ac_R2"
  same "conditions_of (sample list order)" \
    "$(call "${config}" conditions_of H3K27ac | xargs)" "Naive Arm Cl13"
  same "all_marks (without IgG)" "$(call "${config}" all_marks | xargs)" \
    "H3K27ac"
  same "peak_setting of a mark" \
    "$(call "${config}" peak_setting H3K9me3 min_fold)" "1.2"
  same "peak_setting default row" \
    "$(call "${config}" peak_setting NewMark min_width)" "100"

  # Sample list saved by Excel on Windows (CRLF line ends), with a comment
  local windows_list="${test_dir}/sample_list_windows.txt"
  { echo "# comment"; sed 's/$/\r/' "${test_dir}/sample_list.txt"; } \
    > "${windows_list}"
  make_config "${windows_list}" "${test_dir}/config_windows.sh"
  same "Windows line ends and comments" \
    "$(call "${test_dir}/config_windows.sh" sample_info IgG_R2 control)" "-"
}

test_reference() {
  echo "3. Reference (build_reference.sh)"
  check "build_reference.sh runs" bash -c "bash \
'${PIPELINE_DIR}/reference/build_reference.sh' mm10 '${test_dir}/indexes' 2 \
> '${test_dir}/build_reference.log' 2>&1"
  local gnm_idx="${test_dir}/indexes/mm10"
  check "bowtie2 index" test -f "${gnm_idx}/bowtie2/mm10.1.bt2"
  same "chr1toY has no chrM" \
    "$(cut -f 1 "${gnm_idx}/mm10.chr1toY.bed" 2> /dev/null | xargs)" \
    "chr1 chr2"
  same "genome sizes" "$(count_lines "${gnm_idx}/mm10.genomesize")" 3
  same "gene table" "$(count_lines "${gnm_idx}/genes.txt")" 21
  check "TxDb" test -s "${gnm_idx}/genes.txdb.sqlite"
}

test_input_checks() {
  echo "4. Input checks (validate.sh)"
  check "good input passes" bash "${PIPELINE_DIR}/scripts/validate.sh" \
    "${config}"
  local igg_fastq="IgG_R1_R1.fastq.gz\tIgG_R1_R2.fastq.gz"
  check_catches "comma-separated file" "header must be" header
  check_catches "duplicate sample ID" "used more than once" \
    "IgG_R1\t${igg_fastq}\tIgG\tIgG\t1\t-\tb1"
  check_catches "space in sample ID" "may only have" \
    "Bad name\t${igg_fastq}\tIgG\tIgG\t3\t-\tb1"
  check_catches "missing fastq" "not found" \
    "New_R1\tmissing_R1.fastq.gz\tmissing_R2.fastq.gz\tH3K27ac\tArm\t3\t-\tb1"
  check_catches "unknown control" "not in the sample list" \
    "New_R1\t${igg_fastq}\tH3K27ac\tArm\t3\tIgG_R9\tb1"
  check_catches "missing column" "8 tab-separated columns" \
    "New_R1\t${igg_fastq}\tH3K27ac\tArm\t3"
}

test_pipeline() {
  echo "5. Whole pipeline on simulated data"
  check "run_pipeline.sh finishes" bash -c "bash \
'${PIPELINE_DIR}/run_pipeline.sh' '${config}' > '${test_dir}/run.log' 2>&1"
  local data="${test_dir}/Data"
  local stats="${data}/BAM/Arm_H3K27ac_R1.align_stats.txt"

  # FastQC and MultiQC of raw fastq (8 samples, one with 2 lanes: 18 files)
  same "FastQC of every raw fastq file" "$(find "${data}/QC/FastQC_raw" \
    -name '*_fastqc.html' 2> /dev/null | wc -l | tr -d ' ')" 18
  check "MultiQC of raw fastq" \
    test -s "${data}/QC/MultiQC/multiqc_raw_fastq.html"
  check "MultiQC of all steps" test -s "${data}/QC/MultiQC/multiqc_report.html"

  # Trimming, lanes, alignment, bigwig
  same "fastp report of every sample" \
    "$(find "${data}/QC" -name '*.fastp.json' | wc -l | tr -d ' ')" 8
  same "FastQC of trimmed reads" "$(find "${data}/QC/FastQC_trimmed" \
    -name '*_fastqc.html' 2> /dev/null | wc -l | tr -d ' ')" 16
  local raw_pairs
  raw_pairs="$(gzip -dc "${test_dir}"/Raw/Naive_H3K27ac_R1_L*_R1.fastq.gz \
    | awk 'END { print NR / 4 }')"
  same "2 lanes merged (read pairs)" "$(table_value \
    "${data}/BAM/Naive_H3K27ac_R1.align_stats.txt" Naive_H3K27ac_R1 \
    input_pairs)" "${raw_pairs}"
  between "alignment rate %" \
    "$(table_value "${stats}" Arm_H3K27ac_R1 align_rate | cut -d. -f1)" 90 100
  same "no chrM reads after filtering" "$(run_in tools samtools idxstats \
    "${data}/BAM/Arm_H3K27ac_R1.map.dedup.bam" 2> /dev/null \
    | awk '$1 == "chrM" { print $3 }')" 0
  between "fragment size peak (bp)" "$(sort -k2,2nr \
    "${data}/BAM/Arm_H3K27ac_R1.fragment_size.txt" 2> /dev/null \
    | awk 'NR == 1 { print $1 }')" 120 180
  same "bigwig of every sample" \
    "$(find "${data}/Wig" -name '*.rpm.bw' | wc -l | tr -d ' ')" 8

  # Peaks: 100 simulated sites; sites 1-10 are also in IgG
  between "MACS2 peaks of Arm_R1" \
    "$(count_lines "${data}/Peak/peak_Arm_H3K27ac_R1.macs2.bed")" 90 110
  between "peaks for the analysis (Arm_R1)" \
    "$(count_lines "${data}/Peak/peak_Arm_H3K27ac_R1.final.bed")" 90 110
  if [[ "${peak_callers}" == *flexpeak* ]]; then
    between "FlexPeak peaks of Arm_R1" \
      "$(count_lines "${data}/Peak/peak_Arm_H3K27ac_R1.flexpeak.bed")" 90 110
    between "MACS2 vs FlexPeak overlap (Jaccard %)" "$(table_value \
      "${data}/Peak/peak_Arm_H3K27ac_R1.callers.txt" Arm_H3K27ac_R1 jaccard \
      | awk '{ printf "%d", $1 * 100 }')" 50 100
  else
    echo "  SKIP  FlexPeak (not installed in the test environment)"
  fi
  between "union peaks (sites 11-100)" "$(count_lines \
    "${data}/Consensus/H3K27ac/peak_H3K27ac_union.bed")" 86 94
  between "annotated union peaks" "$(awk -F'\t' 'NR == 1 {
      for (i = 1; i <= NF; i++) if ($i == "Anno") a = i; next }
      $a != "NA" { n++ } END { print n + 0 }' \
    "${data}/Annotation/H3K27ac/peak_H3K27ac_union.anno.txt")" 86 94

  # Count matrices: union peaks x 6 samples
  local counts="${data}/ReadCount/H3K27ac"
  local n_union
  n_union="$(count_lines "${data}/Consensus/H3K27ac/peak_H3K27ac_union.bed")"
  same "raw matrix: peaks" "$(count_lines "${counts}/RawRC_H3K27ac.txt")" \
    "$(( n_union + 1 ))"
  check "raw matrix: sample columns" grep -q \
    $'Naive_H3K27ac_R1\tNaive_H3K27ac_R2\tArm_H3K27ac_R1' \
    "${counts}/RawRC_H3K27ac.txt"
  check "normalized matrix" test -s "${counts}/NormRC_H3K27ac.txt"
  check "VST matrix" test -s "${counts}/VST_H3K27ac.txt"
  between "FRiP of Arm_R1 (%)" "$(table_value \
    "${counts}/sizeFactors_H3K27ac.txt" Arm_H3K27ac_R1 FRiP \
    | awk '{ printf "%d", $1 * 100 }')" 30 90

  # Differential peaks: every 7th site lower in Naive (13 sites), every 10th
  # site higher in Cl13 (9 sites), site 70 both
  local summary="${data}/DiffPeak/H3K27ac/deseq2_summary_H3K27ac.txt"
  between "Arm vs Naive: Up" "$(table_value "${summary}" Arm_vs_Naive Up)" \
    11 15
  between "Cl13 vs Arm: Up" "$(table_value "${summary}" Cl13_vs_Arm Up)" 7 11
  between "Cl13 vs Naive: Up" "$(table_value "${summary}" Cl13_vs_Naive Up)" \
    18 24
  same "no Down peaks (none simulated)" "$(awk -F'\t' 'NR > 1 { print $NF }' \
    "${summary}" 2> /dev/null | xargs)" "0 0 0"

  # Figures and HTML report
  local figures="${data}/Figures/H3K27ac"
  check "signal matrix (deepTools)" \
    test -s "${figures}/signal_H3K27ac.matrix.gz"
  check "summary figure" test -s "${figures}/Figure_H3K27ac_summary.pdf"
  check_report "${data}/Report/CUTRUN_report.html" "${n_union}"

  # QC summary and reproducibility records
  same "QC summary of all samples" \
    "$(count_lines "${data}/QC/Summary/sample_QC_summary.txt")" 9
  check "software versions" grep -q bowtie2 "${data}/Log/software_versions.txt"
  check "settings of the run" compgen -G "${data}/Log/run_*/config.sh"
  check "R sessionInfo" \
    test -s "${data}/DiffPeak/H3K27ac/sessionInfo_H3K27ac.txt"
}

test_reruns() {
  echo "6. Reruns"
  local data="${test_dir}/Data"
  check "rerun finishes" bash -c "bash '${PIPELINE_DIR}/run_pipeline.sh' \
'${config}' > '${test_dir}/rerun.log' 2>&1"
  check "rerun skips finished alignment" grep -q "align: already done, skip" \
    "${data}/Log/MyCutRun_align_Arm_H3K27ac_R1.err"

  local summary="${data}/DiffPeak/H3K27ac/deseq2_summary_H3K27ac.txt"
  rm -f "${summary}"
  check "diff part runs on its own" bash -c "bash \
'${PIPELINE_DIR}/run_pipeline.sh' '${config}' diff \
> '${test_dir}/diff.log' 2>&1"
  between "diff part: Cl13 vs Arm Up" \
    "$(table_value "${summary}" Cl13_vs_Arm Up)" 7 11
}

main() {
  test_dir="${1:-${TESTS_DIR}/test_output}"
  # Only remove a test folder made by this script
  if [[ -f "${test_dir}/.cutrun_test" ]]; then rm -rf "${test_dir}"; fi
  if [[ -e "${test_dir}" ]]; then
    echo "${test_dir} exists and is not a test folder" >&2
    exit 1
  fi
  mkdir -p "${test_dir}"
  touch "${test_dir}/.cutrun_test"
  echo "Test folder: ${test_dir}"

  bash "${TESTS_DIR}/make_test_data.sh" "${test_dir}" > /dev/null \
    || { echo "Could not make the test data" >&2; exit 1; }
  config="${test_dir}/config.sh"
  make_peak_settings
  make_config "${test_dir}/sample_list.txt" "${config}"
  # Test FlexPeak too if it is installed in the tools environment
  if run_in tools command -v flexpeak > /dev/null 2>&1; then
    peak_callers="macs2 flexpeak"
    make_config "${test_dir}/sample_list.txt" "${config}"
  fi
  echo "Peak callers: ${peak_callers}"

  test_syntax
  test_helpers
  test_reference
  test_input_checks
  test_pipeline
  test_reruns

  echo
  echo "${n_pass} passed, ${n_fail} failed"
  if (( n_fail > 0 )); then
    echo "Logs: ${test_dir}/*.log and ${test_dir}/Data/Log/"
    exit 1
  fi
}

main "$@"
