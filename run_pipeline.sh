#!/bin/bash
#
# Runs the CUT&RUN pipeline: checks the inputs, then submits the jobs.
#
# Usage: run_pipeline.sh <config.sh> [part]
#   part: all (default) | qc | align | peaks | matrix | diff
#     qc     : FastQC and MultiQC of the raw fastq files       (one job)
#     align  : trimming, alignment, bigwig     (one job per sample)
#     peaks  : peak calling, MACS2 / FlexPeak  (one job per sample)
#     matrix : consensus and union peaks, annotation, raw and normalized
#              count matrices, QC summary      (one job)
#     diff   : differential peaks, figures, HTML report        (one job)
# Finished steps are skipped: after fixing a problem, run the same command.
# Check progress: check_status.sh <config.sh>

source "$(dirname "${BASH_SOURCE[0]}")/scripts/common.sh"

#######################################
# Prints the usage (the header comment of this file).
#######################################
usage() {
  sed -n '3,15p' "$0" | sed 's/^# \{0,1\}//'
}

#######################################
# Saves the settings and the pipeline version of this run in Log/run_<date>/.
# Globals:
#   CONFIG, SAMPLE_LIST, PEAK_SETTINGS, COMPARISON_LIST, WORKDIR, PIPELINE_DIR
#######################################
save_run_settings() {
  local run_log
  run_log="${WORKDIR}/Log/run_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "${run_log}"
  cp "${CONFIG}" "${SAMPLE_LIST}" "${PEAK_SETTINGS}" "${run_log}/"
  if [[ -f "${COMPARISON_LIST}" ]]; then
    cp "${COMPARISON_LIST}" "${run_log}/"
  fi

  local version commit
  version="$(cat "${PIPELINE_DIR}/VERSION")"
  commit="$(git -C "${PIPELINE_DIR}" rev-parse --short HEAD 2> /dev/null \
    || echo "not a git copy")"
  echo "CUT&RUN pipeline ${version} (${commit})" \
    > "${run_log}/pipeline_version.txt"
  log "Pipeline ${version} (${commit}). Settings saved in ${run_log}"
}

#######################################
# Prints the job ID of a sample.
# Arguments:
#   Sample ID, list of "sample_id:job_id" (space-separated).
#######################################
job_of() {
  local pair
  for pair in $2; do
    if [[ "${pair%%:*}" == "$1" ]]; then echo "${pair#*:}"; fi
  done
}

main() {
  if (( $# < 1 )); then
    usage
    exit 1
  fi
  read_config "$1"
  local part="${2:-all}"
  case "${part}" in
    all | qc | align | peaks | matrix | diff) ;;
    *) die "Unknown part '${part}'. Use all, qc, align, peaks, matrix, diff." ;;
  esac

  bash "${CODE_DIR}/validate.sh" "${CONFIG}" \
    || die "Please fix the errors above. Nothing was submitted."
  save_run_settings

  local samples=() sample_id
  while read -r sample_id; do samples+=("${sample_id}"); done < <(all_samples)

  local align_jobs="" peak_jobs="" matrix_job="" job_id wait_for control

  # 0. One job: FastQC and MultiQC of the raw fastq (runs next to the others)
  if [[ "${part}" == all || "${part}" == qc ]]; then
    job_id="$(submit_job "${PROJECT_NAME}_qc" "" run_all_samples.sh \
      "${CONFIG}" qc)"
    log "Submitted FastQC/MultiQC of raw fastq (${job_id})"
  fi

  # 1. One job per sample: trimming, alignment, bigwig
  if [[ "${part}" == all || "${part}" == align ]]; then
    for sample_id in "${samples[@]}"; do
      job_id="$(submit_job "${PROJECT_NAME}_align_${sample_id}" "" \
        run_sample.sh "${CONFIG}" "${sample_id}" align)"
      align_jobs+=" ${sample_id}:${job_id}"
      log "Submitted alignment of ${sample_id} (${job_id})"
    done
  fi

  # 2. One job per sample: peak calling. Waits for the alignment of the
  # sample, and of its IgG with IGG_MODE=control.
  if [[ "${part}" == all || "${part}" == peaks ]]; then
    for sample_id in "${samples[@]}"; do
      wait_for="$(job_of "${sample_id}" "${align_jobs}")"
      control="$(sample_info "${sample_id}" control)"
      if [[ "${IGG_MODE}" == control && "${control}" != "-" ]]; then
        wait_for+=" $(job_of "${control}" "${align_jobs}")"
      fi
      job_id="$(submit_job "${PROJECT_NAME}_peaks_${sample_id}" \
        "${wait_for}" run_sample.sh "${CONFIG}" "${sample_id}" peaks)"
      peak_jobs+=" ${job_id}"
      log "Submitted peak calling of ${sample_id} (${job_id})"
    done
  fi

  # 3. One job: count matrices of all samples (waits for all peak calling)
  if [[ "${part}" == all || "${part}" == matrix ]]; then
    matrix_job="$(submit_job "${PROJECT_NAME}_matrix" "${peak_jobs}" \
      run_all_samples.sh "${CONFIG}" matrix)"
    log "Submitted count matrices (${matrix_job})"
  fi

  # 4. One job: differential peaks, figures and HTML report (waits for the
  # count matrices)
  if [[ "${part}" == all || "${part}" == diff ]]; then
    job_id="$(submit_job "${PROJECT_NAME}_diff" "${matrix_job}" \
      run_all_samples.sh "${CONFIG}" diff)"
    log "Submitted differential peaks, figures and report (${job_id})"
  fi

  log "All submitted. Logs: ${WORKDIR}/Log"
  log "Check progress: bash ${PIPELINE_DIR}/check_status.sh ${CONFIG}"
}

main "$@"
