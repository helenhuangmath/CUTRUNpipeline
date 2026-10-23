#!/bin/bash
#
# Runs the steps of one sample. Finished steps are skipped.
#
# Usage: run_sample.sh <config.sh> <sample_id> align   (trim, align, bigwig)
#        run_sample.sh <config.sh> <sample_id> peaks   (peak calling)

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local sample_id="$2"
  local part="$3"

  case "${part}" in
    align)
      run_step "${sample_id}.trim" 01_trim.sh "${CONFIG}" "${sample_id}"
      run_step "${sample_id}.align" 02_align.sh "${CONFIG}" "${sample_id}"
      run_step "${sample_id}.bigwig" 03_bigwig.sh "${CONFIG}" "${sample_id}"
      ;;
    peaks)
      run_step "${sample_id}.peaks" 04_peaks.sh "${CONFIG}" "${sample_id}"
      ;;
    *)
      die "Part must be align or peaks, not '${part}'"
      ;;
  esac
}

main "$@"
