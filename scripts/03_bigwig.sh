#!/bin/bash
#
# Step 3: normalized bigwig track of one sample, for IGV / UCSC.
#
# Usage: 03_bigwig.sh <config.sh> <sample_id>
# Output: Wig/<sample_id>.rpm.bw (CPM) or Wig/<sample_id>.spikein.bw

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p Wig

  # 3.1 Normalization: CPM (counts per million fragments) or spike-in
  local norm_args=(--normalizeUsing CPM)
  local bw="Wig/${sample_id}.rpm.bw"
  if [[ "${BW_NORM}" == spikein ]]; then
    local spikein_pairs scale
    spikein_pairs="$(awk 'NR == 2 { print $7 }' \
      "BAM/${sample_id}.align_stats.txt")"
    scale="$(awk -v c="${SPIKEIN_SCALE}" -v n="${spikein_pairs}" \
      'BEGIN { print c / n }')"
    norm_args=(--normalizeUsing None --scaleFactor "${scale}")
    bw="Wig/${sample_id}.spikein.bw"
  fi

  # 3.2 bamCoverage
  #   --extendReads --samFlagInclude 64 : each fragment counted once (read 1),
  #                                       over the whole fragment
  #   -bl : skip blacklist regions
  local bw_args
  read -ra bw_args <<< "${BW_ARGS}"
  bamCoverage -b "BAM/${sample_id}.map.dedup.bam" -o "${bw}" \
    --extendReads --samFlagInclude 64 -bl "${BLACKLIST}" -p "${THREADS}" \
    "${norm_args[@]}" "${bw_args[@]}" 2> "Wig/${sample_id}.bamCoverage.log"
}

main "$@"
