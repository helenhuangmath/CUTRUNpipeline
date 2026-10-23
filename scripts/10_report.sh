#!/bin/bash
#
# Step 10: interactive HTML report of all marks: search a gene to see its
# peaks, differential statistics, signal profile and normalized counts.
# The report is one self-contained file that opens offline in a browser.
#
# Usage: 10_report.sh <config.sh>
# Output: Report/CUTRUN_report.html

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  activate_env R
  Rscript "${R_CODE_DIR}/report.R" \
    --workdir="${WORKDIR}" \
    --sample-list="${SAMPLE_LIST}" \
    --igg-mark="${IGG_MARK}" \
    --species="${SPECIES}" \
    --project="${PROJECT_NAME}" \
    --version="$(cat "${PIPELINE_DIR}/VERSION")" \
    --out="${WORKDIR}/Report/CUTRUN_report.html"
}

main "$@"
