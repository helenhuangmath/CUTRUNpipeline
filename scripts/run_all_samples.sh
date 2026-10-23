#!/bin/bash
#
# Runs the steps that use all samples. These steps always run again, so new
# samples or changed settings are picked up.
#
# Usage: run_all_samples.sh <config.sh> qc
#          FastQC and MultiQC of the raw fastq files
#        run_all_samples.sh <config.sh> matrix
#          consensus and union peaks, annotation, count matrices, QC summary
#        run_all_samples.sh <config.sh> diff
#          differential peaks (DESeq2), summary figures, HTML report

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

#######################################
# Prints the first line with a number of a program's version message.
# Arguments:
#   Program and its version option, e.g. "bowtie2 --version".
#######################################
version_of() {
  if ! command -v "$1" > /dev/null; then
    echo "not found"
    return 0
  fi
  local message
  message="$("$@" 2>&1 || true)"
  echo "${message}" | grep '[0-9]' | head -1 || true
}

#######################################
# Writes the software versions used, for reproducibility and methods text.
# Globals:
#   WORKDIR, PIPELINE_DIR
#######################################
write_software_versions() {
  local out="${WORKDIR}/Log/software_versions.txt"
  activate_env tools
  {
    echo "pipeline      $(cat "${PIPELINE_DIR}/VERSION")"
    echo "fastqc        $(version_of fastqc --version)"
    echo "fastp         $(version_of fastp --version)"
    echo "bowtie2       $(version_of bowtie2 --version)"
    echo "samtools      $(version_of samtools --version)"
    echo "deeptools     $(version_of bamCoverage --version)"
    echo "macs2         $(version_of macs2 --version)"
    echo "bedtools      $(version_of bedtools --version)"
    echo "featureCounts $(version_of featureCounts -v)"
  } > "${out}"
  activate_env R
  echo "R             $(version_of R --version)" >> "${out}"
}

#######################################
# FastQC of every raw fastq file (each lane) and a MultiQC report, to check
# the data before (or while) it is processed.
# Globals:
#   WORKDIR, SAMPLE_LIST, DATADIR, THREADS
#######################################
run_qc_part() {
  activate_env tools
  local out_dir="${WORKDIR}/QC/FastQC_raw"
  mkdir -p "${out_dir}" "${WORKDIR}/QC/MultiQC"

  local fastqs=() sample_id read lane lanes
  while read -r sample_id; do
    for read in 1 2; do
      IFS=',' read -ra lanes <<< "$(sample_info "${sample_id}" "fastq${read}")"
      for lane in "${lanes[@]}"; do fastqs+=("$(full_path "${lane}")"); done
    done
  done < <(all_samples)

  log "FastQC of ${#fastqs[@]} raw fastq files"
  fastqc -q -t "${THREADS}" -o "${out_dir}" "${fastqs[@]}"
  multiqc -f -q -o "${WORKDIR}/QC/MultiQC" -n multiqc_raw_fastq.html \
    "${out_dir}"
  touch "${WORKDIR}/Status/qc.done"
}

#######################################
# Consensus peaks, annotation and count matrices of each mark; QC summary.
# Globals:
#   CONFIG, WORKDIR, SAMPLE_LIST, IGG_MARK
#######################################
run_matrix_part() {
  local marks=() mark
  while read -r mark; do marks+=("${mark}"); done < <(all_marks)

  for mark in "${marks[@]}"; do
    log "---------- ${mark} ----------"
    bash "${CODE_DIR}/05_consensus.sh" "${CONFIG}" "${mark}"
    bash "${CODE_DIR}/06_annotate.sh" "${CONFIG}" "${mark}"
    bash "${CODE_DIR}/07_count_matrix.sh" "${CONFIG}" "${mark}"
    touch "${WORKDIR}/Status/${mark}.matrix.done"
  done

  # QC summary of all samples: MultiQC (raw and trimmed FastQC, fastp,
  # bowtie2, samtools), table and plots, software versions
  activate_env tools
  multiqc -f -q -o "${WORKDIR}/QC/MultiQC" -n multiqc_report.html \
    "${WORKDIR}/QC" "${WORKDIR}/BAM" || log "MultiQC failed, continue"
  write_software_versions
  activate_env R
  Rscript "${R_CODE_DIR}/qc_summary.R" --workdir="${WORKDIR}" \
    --sample-list="${SAMPLE_LIST}" --igg-mark="${IGG_MARK}"
  touch "${WORKDIR}/Status/qc_summary.done"
}

#######################################
# Differential peaks and summary figure of each mark; HTML report.
# Globals:
#   CONFIG, WORKDIR
#######################################
run_diff_part() {
  local marks=() mark
  while read -r mark; do marks+=("${mark}"); done < <(all_marks)

  for mark in "${marks[@]}"; do
    log "---------- ${mark} ----------"
    bash "${CODE_DIR}/08_diff.sh" "${CONFIG}" "${mark}"
    bash "${CODE_DIR}/09_figures.sh" "${CONFIG}" "${mark}"
    touch "${WORKDIR}/Status/${mark}.diff.done"
  done
  bash "${CODE_DIR}/10_report.sh" "${CONFIG}"
  touch "${WORKDIR}/Status/report.done"
}

main() {
  read_config "$1"
  case "$2" in
    qc) run_qc_part ;;
    matrix) run_matrix_part ;;
    diff) run_diff_part ;;
    *) die "Part must be qc, matrix or diff, not '$2'" ;;
  esac
  log "Finished ($2). Results in ${WORKDIR}"
}

main "$@"
