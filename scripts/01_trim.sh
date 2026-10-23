#!/bin/bash
#
# Step 1: adapter trimming (fastp) and quality check of the trimmed reads
# (FastQC) of one sample. Trimming matters for CUT&RUN: many fragments are
# shorter than the read length, so reads run into the adapter.
# (FastQC of the raw reads is done for all samples by the "qc" part.)
#
# Usage: 01_trim.sh <config.sh> <sample_id>
# Output: Fastq/<sample_id>_R1/2.trim.fastq.gz, QC/<sample_id>.fastp.html,
#         QC/FastQC_trimmed/

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p Fastq QC

  # Raw reads (several lanes are merged first)
  local fastq1 fastq2
  fastq1="$(merge_lanes "${sample_id}" 1)"
  fastq2="$(merge_lanes "${sample_id}" 2)"

  # Reads for the alignment
  local trim1="Fastq/${sample_id}_R1.trim.fastq.gz"
  local trim2="Fastq/${sample_id}_R2.trim.fastq.gz"

  # 1.1 Trim adapters and low quality bases
  #   --detect_adapter_for_pe : find adapters from the overlap of R1 and R2
  #   --length_required 20    : drop reads shorter than 20 bp after trimming
  #   -h / -j                 : html report / json summary for MultiQC
  if [[ "${RUN_TRIMMING}" == true ]]; then
    local fastp_args
    read -ra fastp_args <<< "${FASTP_ARGS}"
    fastp -i "${fastq1}" -I "${fastq2}" -o "${trim1}" -O "${trim2}" \
      -w "${THREADS}" -R "${sample_id}" \
      -h "QC/${sample_id}.fastp.html" -j "QC/${sample_id}.fastp.json" \
      "${fastp_args[@]}" 2> "QC/${sample_id}.fastp.log"
  else
    # No trimming: link the raw reads
    ln -sf "${fastq1}" "${trim1}"
    ln -sf "${fastq2}" "${trim2}"
  fi

  # 1.2 Quality check of the trimmed reads
  if [[ "${RUN_FASTQC}" == true && "${RUN_TRIMMING}" == true ]]; then
    mkdir -p QC/FastQC_trimmed
    fastqc -q -t 2 -o QC/FastQC_trimmed "${trim1}" "${trim2}"
  fi
}

main "$@"
