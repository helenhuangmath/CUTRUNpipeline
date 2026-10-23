#!/bin/bash
#
# Checks config.sh, the sample list, the comparison list and the reference
# files before any job is submitted (run_pipeline.sh runs it first).
#
# Usage: validate.sh <config.sh>
# Exit status: 0 if the inputs are good, 1 if there are errors.

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

readonly HEADER="sample fastq1 fastq2 mark condition replicate control batch"
n_errors=0

#######################################
# Prints an input error and counts it.
# Globals:
#   n_errors
# Arguments:
#   Message.
#######################################
input_error() {
  err "$*"
  n_errors=$(( n_errors + 1 ))
}

#######################################
# Prints a warning.
# Arguments:
#   Message.
#######################################
warning() {
  log "WARNING: $*"
}

#######################################
# Returns 0 if a name has only letters, numbers, _ . and -.
# Arguments:
#   Name.
#######################################
is_good_name() {
  [[ "$1" =~ ^[A-Za-z0-9_.-]+$ ]]
}

#######################################
# Checks the header, columns and sample IDs of the sample list.
# Globals:
#   SAMPLE_LIST
#######################################
check_sample_list_format() {
  local header
  header="$(read_table "${SAMPLE_LIST}" | head -1 | tr '\t' ' ')"
  if [[ "${header}" != "${HEADER}" ]]; then
    die "The sample list header must be (tab-separated): ${HEADER}"
  fi

  local row
  while read -r row; do
    input_error "Sample ${row} does not have 8 tab-separated columns"
  done < <(read_table "${SAMPLE_LIST}" \
    | awk -F'\t' 'NR > 1 && NF != 8 { print $1 }')

  local sample_id
  while IFS= read -r sample_id; do
    is_good_name "${sample_id}" \
      || input_error "Sample ID '${sample_id}' may only have letters," \
        "numbers, _ . -"
  done < <(read_table "${SAMPLE_LIST}" | awk -F'\t' 'NR > 1 { print $1 }')

  while read -r sample_id; do
    input_error "Sample ID ${sample_id} is used more than once"
  done < <(all_samples | sort | uniq -d)

  (( n_errors == 0 )) || die "Please fix the sample list first."
}

#######################################
# Checks names, fastq files and the control of every sample.
# Globals:
#   SAMPLE_LIST, DATADIR
#######################################
check_samples() {
  local samples=() sample_id column value read fastq control
  while read -r sample_id; do samples+=("${sample_id}"); done < <(all_samples)

  for sample_id in "${samples[@]}"; do
    for column in mark condition; do
      value="$(sample_info "${sample_id}" "${column}")"
      is_good_name "${value}" || input_error "${sample_id}: ${column}" \
        "'${value}' may only have letters, numbers, _ . -"
    done

    for read in 1 2; do
      local fastqs
      IFS=',' read -ra fastqs <<< "$(sample_info "${sample_id}" "fastq${read}")"
      for fastq in "${fastqs[@]}"; do
        [[ "${fastq}" == *.gz ]] \
          || input_error "${sample_id}: fastq must be gzipped: ${fastq}"
        [[ -f "$(full_path "${fastq}")" ]] || input_error "${sample_id}:" \
          "fastq${read} not found: $(full_path "${fastq}")"
      done
    done

    control="$(sample_info "${sample_id}" control)"
    if [[ "${control}" != "-" && " ${samples[*]} " != *" ${control} "* ]]; then
      input_error "${sample_id}: control ${control} is not in the sample list"
    fi
  done
}

#######################################
# Checks marks, replicates and IgG samples; prints a summary.
# Globals:
#   PEAK_SETTINGS, MIN_REPS, IGG_MODE, IGG_MARK
#######################################
check_marks() {
  local marks=() mark cond n_reps summary
  while read -r mark; do marks+=("${mark}"); done < <(all_marks)
  (( ${#marks[@]} > 0 )) || input_error "No samples to analyse (only IgG)"

  for mark in "${marks[@]}"; do
    grep -q "^${mark}"$'\t' "${PEAK_SETTINGS}" || warning "${mark} is not" \
      "in peak_settings.txt; the default settings are used"
    summary=""
    while read -r cond; do
      n_reps="$(samples_of "${mark}" "${cond}" | wc -l | tr -d ' ')"
      summary+=" ${cond}=${n_reps}"
      if (( n_reps < MIN_REPS )); then
        warning "${mark} ${cond} has ${n_reps} replicate(s);" \
          "its consensus peaks use all of them"
      fi
    done < <(conditions_of "${mark}")
    log "${mark}:${summary}"
  done

  if [[ "${IGG_MODE}" != none && -z "$(samples_of "${IGG_MARK}")" ]]; then
    warning "No IgG samples (mark ${IGG_MARK}); the IgG step is skipped"
  fi
}

#######################################
# Checks that both conditions of each comparison exist for the mark.
# Globals:
#   COMPARISON_LIST
#######################################
check_comparisons() {
  if [[ -z "${COMPARISON_LIST}" ]]; then
    warning "COMPARISON_LIST is empty: no differential peaks"
    return 0
  fi
  if [[ ! -f "${COMPARISON_LIST}" ]]; then
    input_error "Comparison list not found: ${COMPARISON_LIST}"
    return 0
  fi

  local comparison_mark treatment reference mark conditions
  while read -r comparison_mark treatment reference; do
    local marks=("${comparison_mark}")
    if [[ "${comparison_mark}" == "*" ]]; then
      marks=()
      while read -r mark; do marks+=("${mark}"); done < <(all_marks)
    fi
    for mark in "${marks[@]}"; do
      conditions=" $(conditions_of "${mark}" | tr '\n' ' ')"
      if [[ "${conditions}" != *" ${treatment} "* \
        || "${conditions}" != *" ${reference} "* ]]; then
        warning "${mark}: ${treatment} vs ${reference} is skipped," \
          "a condition is not found"
      fi
    done
  done < <(read_table "${COMPARISON_LIST}" \
    | awk -F'\t' 'NR > 1 { print $1, $2, $3 }')
}

#######################################
# Checks the peak callers and which peaks are used for the analysis.
# Globals:
#   PEAK_CALLERS, ANALYSIS_PEAKS
#######################################
check_peak_callers() {
  local caller
  for caller in ${PEAK_CALLERS}; do
    [[ "${caller}" == macs2 || "${caller}" == flexpeak ]] \
      || input_error "PEAK_CALLERS: unknown peak caller '${caller}'"
  done
  case "${ANALYSIS_PEAKS}" in
    macs2 | flexpeak)
      [[ " ${PEAK_CALLERS} " == *" ${ANALYSIS_PEAKS} "* ]] || input_error \
        "ANALYSIS_PEAKS=${ANALYSIS_PEAKS} but it is not in PEAK_CALLERS"
      ;;
    both)
      [[ " ${PEAK_CALLERS} " == *" macs2 "* \
        && " ${PEAK_CALLERS} " == *" flexpeak "* ]] \
        || input_error "ANALYSIS_PEAKS=both needs both peak callers in" \
          "PEAK_CALLERS"
      ;;
    *)
      input_error "ANALYSIS_PEAKS must be macs2, flexpeak or both"
      ;;
  esac
}

#######################################
# Checks the reference files and the scheduler.
# Globals:
#   BOWTIE_IDX, GNM_SIZES, BLACKLIST, CHR1TOY_BED, TXDB, SPIKEIN_IDX,
#   BW_NORM, SCHEDULER
#######################################
check_reference_and_cluster() {
  [[ -f "${BOWTIE_IDX}.1.bt2" || -f "${BOWTIE_IDX}.1.bt2l" ]] \
    || input_error "bowtie2 index not found: ${BOWTIE_IDX}"
  [[ -f "${GNM_SIZES}" ]] || input_error "GNM_SIZES not found: ${GNM_SIZES}"
  [[ -f "${BLACKLIST}" ]] || input_error "BLACKLIST not found: ${BLACKLIST}"
  [[ -f "${CHR1TOY_BED}" ]] \
    || input_error "CHR1TOY_BED not found: ${CHR1TOY_BED}"
  [[ -f "${TXDB}" ]] || warning "TXDB not found: no ChIPseeker annotation"

  if [[ -n "${SPIKEIN_IDX}" && ! -f "${SPIKEIN_IDX}.1.bt2" ]]; then
    input_error "spike-in index not found: ${SPIKEIN_IDX}"
  fi
  if [[ "${BW_NORM}" == spikein && -z "${SPIKEIN_IDX}" ]]; then
    input_error "BW_NORM=spikein needs SPIKEIN_IDX"
  fi

  case "${SCHEDULER}" in
    lsf) command -v bsub > /dev/null || input_error "bsub not found" ;;
    slurm) command -v sbatch > /dev/null || input_error "sbatch not found" ;;
    local) ;;
    *) input_error "SCHEDULER must be lsf, slurm or local" ;;
  esac
}

main() {
  read_config "$1"
  [[ -f "${SAMPLE_LIST}" ]] || die "Sample list not found: ${SAMPLE_LIST}"
  [[ -f "${PEAK_SETTINGS}" ]] || die "Peak settings not found: ${PEAK_SETTINGS}"

  check_sample_list_format
  check_samples
  check_marks
  check_comparisons
  check_peak_callers
  check_reference_and_cluster

  (( n_errors == 0 )) || die "${n_errors} error(s). Please fix them."
  log "Input check passed: $(all_samples | wc -l | tr -d ' ') samples"
}

main "$@"
