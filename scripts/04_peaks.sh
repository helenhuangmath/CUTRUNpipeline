#!/bin/bash
#
# Step 4: peak calling of one sample with MACS2 and/or FlexPeak
# (PEAK_CALLERS), using the settings of its mark (config/peak_settings.txt).
# Weak peaks and blacklist peaks are removed.
#
# Usage: 04_peaks.sh <config.sh> <sample_id>
# Output: Peak/MACS2/, Peak/FlexPeak/    raw output of each peak caller
#         Peak/peak_<id>.macs2.bed       cleaned peaks of each caller; columns:
#         Peak/peak_<id>.flexpeak.bed    chr, start, end, name, score, strand,
#                                        fold enrichment, -log10(p), -log10(q)
#         Peak/peak_<id>.final.bed       peaks used for the analysis
#                                        (ANALYSIS_PEAKS)
#         Peak/peak_<id>.callers.txt     MACS2 vs FlexPeak: numbers, overlap

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

#######################################
# Calls peaks with MACS2.
# Arguments:
#   Sample ID, mark, IgG control BAM (or "").
# Outputs:
#   Path of the MACS2 peak file.
#######################################
call_macs2() {
  local sample_id="$1"
  local macs2_args
  read -ra macs2_args <<< "$(peak_setting "$2" macs2_options)"
  if [[ -n "$3" ]]; then macs2_args+=(-c "$3"); fi
  mkdir -p Peak/MACS2

  #   -f BAMPE       : use the real fragment of each read pair
  #   --keep-dup all : duplicates were already handled in step 2
  log "${sample_id}: macs2 ${macs2_args[*]}"
  macs2 callpeak -t "BAM/${sample_id}.map.dedup.bam" \
    -f "${MACS2_FORMAT}" -g "${REFSIZE}" --keep-dup all "${macs2_args[@]}" \
    --outdir Peak/MACS2 -n "peak_${sample_id}" \
    2> "Peak/MACS2/peak_${sample_id}.macs2.log"

  local peak_type=narrowPeak
  if [[ " ${macs2_args[*]} " == *" --broad "* ]]; then peak_type=broadPeak; fi
  echo "Peak/MACS2/peak_${sample_id}_peaks.${peak_type}"
}

#######################################
# Calls peaks with FlexPeak (narrow peaks and broad domains in one run; the
# broadPeak file has one line per region).
# Arguments:
#   Sample ID, mark, IgG control BAM (or "").
# Outputs:
#   Path of the FlexPeak peak file.
#######################################
call_flexpeak() {
  local sample_id="$1"
  local flexpeak_args
  read -ra flexpeak_args <<< "$(peak_setting "$2" flexpeak_options)"
  if [[ -n "$3" ]]; then flexpeak_args+=(-c "$3"); fi
  mkdir -p Peak/FlexPeak

  #   --mapq      : reads were already filtered in step 2
  #   --cache-dir : coverage cache, kept out of the BAM folder
  log "${sample_id}: flexpeak ${flexpeak_args[*]}"
  flexpeak callpeak -t "BAM/${sample_id}.map.dedup.bam" "${flexpeak_args[@]}" \
    --mapq "${MIN_MAPQ}" --blacklist "${BLACKLIST}" \
    --cache-dir Peak/FlexPeak/cache --outdir Peak/FlexPeak \
    -n "peak_${sample_id}" 2> "Peak/FlexPeak/peak_${sample_id}.flexpeak.log"
  echo "Peak/FlexPeak/peak_${sample_id}_peaks.broadPeak"
}

#######################################
# Removes peaks with low fold enrichment and peaks in blacklist regions.
# Arguments:
#   Peak file (MACS2 / FlexPeak format), minimum fold enrichment, output BED.
#######################################
clean_peaks() {
  awk -v min_fold="$2" 'BEGIN { OFS = "\t" }
      !/^(#|track|browser)/ && $7 >= min_fold {
        print $1, $2, $3, $4, $5, ".", $7, $8, $9
      }' "$1" \
    | bedtools intersect -v -a - -b "${BLACKLIST}" \
    | sort -k1,1 -k2,2n > "$3"
}

#######################################
# Compares the peaks of MACS2 and FlexPeak: numbers and base-pair Jaccard.
# Arguments:
#   Sample ID.
#######################################
compare_callers() {
  local macs2="Peak/peak_$1.macs2.bed"
  local flexpeak="Peak/peak_$1.flexpeak.bed"
  local jaccard
  jaccard="$(bedtools jaccard -a "${macs2}" -b "${flexpeak}" \
    | awk 'NR == 2 { print $3 }')"
  {
    printf 'sample\tmacs2_peaks\tflexpeak_peaks\tjaccard\n'
    printf '%s\t%s\t%s\t%s\n' "$1" "$(wc -l < "${macs2}" | tr -d ' ')" \
      "$(wc -l < "${flexpeak}" | tr -d ' ')" "${jaccard}"
  } > "Peak/peak_$1.callers.txt"
  log "$1: MACS2 vs FlexPeak base-pair Jaccard ${jaccard}"
}

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p Peak

  local mark control control_bam="" min_fold
  mark="$(sample_info "${sample_id}" mark)"
  control="$(sample_info "${sample_id}" control)"
  min_fold="$(peak_setting "${mark}" min_fold)"
  if [[ "${IGG_MODE}" == control && "${control}" != "-" ]]; then
    control_bam="BAM/${control}.map.dedup.bam"
  fi

  # 4.1 Call and clean peaks with each peak caller
  local prefix="Peak/peak_${sample_id}"
  local callers caller raw_peaks
  read -ra callers <<< "${PEAK_CALLERS}"
  for caller in "${callers[@]}"; do
    raw_peaks="$("call_${caller}" "${sample_id}" "${mark}" "${control_bam}")"
    clean_peaks "${raw_peaks}" "${min_fold}" "${prefix}.${caller}.bed"
    log "${sample_id} ${caller}: $(awk '!/^#/' "${raw_peaks}" | wc -l) peaks," \
      "$(wc -l < "${prefix}.${caller}.bed") after cleaning"
  done

  # 4.2 Peaks used for the analysis: one caller, or both merged
  if [[ "${ANALYSIS_PEAKS}" == both ]]; then
    cat "${prefix}.macs2.bed" "${prefix}.flexpeak.bed" \
      | sort -k1,1 -k2,2n \
      | bedtools merge -i - > "${prefix}.final.bed"
  else
    cp "${prefix}.${ANALYSIS_PEAKS}.bed" "${prefix}.final.bed"
  fi

  if [[ " ${PEAK_CALLERS} " == *" macs2 "* \
    && " ${PEAK_CALLERS} " == *" flexpeak "* ]]; then
    compare_callers "${sample_id}"
  fi
}

main "$@"
