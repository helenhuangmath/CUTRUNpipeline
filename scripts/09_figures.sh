#!/bin/bash
#
# Step 9: signal around each union peak (deepTools) and a Nature-style
# summary figure of one mark.
#
# Usage: 09_figures.sh <config.sh> <mark>
# Output: Figures/<mark>/signal_<mark>.matrix.gz   signal of every sample,
#                                                  +-2 kb of each peak center
#         Figures/<mark>/Figure_<mark>_summary.pdf  multi-panel figure
#         Figures/<mark>/signal_profile_<mark>.pdf, signal_heatmap_<mark>.pdf

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local mark="$2"
  activate_env tools

  local out_dir="${WORKDIR}/Figures/${mark}"
  mkdir -p "${out_dir}"

  local bw_type=rpm
  if [[ "${BW_NORM}" == spikein ]]; then bw_type=spikein; fi
  local bigwigs=() labels=() sample
  while read -r sample; do
    bigwigs+=("${WORKDIR}/Wig/${sample}.${bw_type}.bw")
    labels+=("${sample}")
  done < <(samples_of "${mark}")

  # 9.1 Signal of each sample around the center of each union peak
  #   -b / -a 2000 : 2 kb before and after the center, in 100-bp bins
  #   --sortRegions keep : same order as the union peaks
  computeMatrix reference-point --referencePoint center -b 2000 -a 2000 \
    --binSize 100 --missingDataAsZero --sortRegions keep \
    -R "${WORKDIR}/Consensus/${mark}/peak_${mark}_union.bed" \
    -S "${bigwigs[@]}" --samplesLabel "${labels[@]}" -p "${THREADS}" \
    -o "${out_dir}/signal_${mark}.matrix.gz" \
    2> "${out_dir}/computeMatrix_${mark}.log"

  # 9.2 Summary figure
  activate_env R
  Rscript "${R_CODE_DIR}/figures.R" --mark="${mark}" --workdir="${WORKDIR}" \
    --sample-list="${SAMPLE_LIST}" --peak-callers="${PEAK_CALLERS}"
}

main "$@"
