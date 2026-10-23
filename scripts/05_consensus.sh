#!/bin/bash
#
# Step 5: consensus peaks of each condition and union peaks of one mark.
#
# Usage: 05_consensus.sh <config.sh> <mark>
# Output: Consensus/<mark>/peak_consensus_<condition>.bed
#         Consensus/<mark>/peak_<mark>_union.bed (Chr Start End peakN Width .)
#         Consensus/<mark>/peak_<mark>_union.saf (for featureCounts)
#         Consensus/<mark>/plot_*.pdf (Venn diagrams, peak width)

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

#######################################
# Writes the IgG peak regions of all IgG samples, extended by IGG_EXTEND bp.
# Globals:
#   WORKDIR, IGG_MARK, IGG_EXTEND, GNM_SIZES
# Arguments:
#   Output BED file.
#######################################
make_igg_regions() {
  local sample
  while read -r sample; do
    cut -f 1-3 "${WORKDIR}/Peak/peak_${sample}.final.bed"
  done < <(samples_of "${IGG_MARK}") \
    | sort -k1,1 -k2,2n \
    | bedtools slop -b "${IGG_EXTEND}" -g "${GNM_SIZES}" -i - \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "$1"
}

main() {
  read_config "$1"
  local mark="$2"
  activate_env tools

  local out_dir="${WORKDIR}/Consensus/${mark}"
  mkdir -p "${out_dir}"
  cd "${out_dir}" || exit 1
  rm -f ./*.bed ./*.saf ./*.txt ./*.pdf  # results of an earlier run

  local consensus_gap union_gap min_width
  consensus_gap="$(peak_setting "${mark}" consensus_gap)"
  union_gap="$(peak_setting "${mark}" union_gap)"
  min_width="$(peak_setting "${mark}" min_width)"
  local conditions=() cond
  while read -r cond; do
    conditions+=("${cond}")
  done < <(conditions_of "${mark}")

  # 5.1 Consensus peaks: in >= MIN_REPS replicates of a condition. Replicate
  # peaks closer than consensus_gap count as one. Column 4 of
  # rep_support.bed = number of replicates with the peak, column 5 = which.
  local sample reps n_reps min_reps
  for cond in "${conditions[@]}"; do
    reps=()
    while read -r sample; do
      reps+=("${sample}")
    done < <(samples_of "${mark}" "${cond}")
    n_reps="${#reps[@]}"
    min_reps=$(( n_reps < MIN_REPS ? n_reps : MIN_REPS ))

    for sample in "${reps[@]}"; do
      awk -v sample="${sample}" 'BEGIN { OFS = "\t" }
          { print $1, $2, $3, sample }' \
        "${WORKDIR}/Peak/peak_${sample}.final.bed"
    done \
      | sort -k1,1 -k2,2n \
      | bedtools merge -d "${consensus_gap}" -c 4,4 \
        -o count_distinct,distinct -i - > "peak_${cond}.rep_support.bed"

    awk -v min_reps="${min_reps}" 'BEGIN { OFS = "\t" }
        $4 >= min_reps { print $1, $2, $3 }' \
      "peak_${cond}.rep_support.bed" > "peak_consensus_${cond}.bed"
    log "${mark} ${cond}: $(wc -l < "peak_consensus_${cond}.bed")" \
      "consensus peaks (in >= ${min_reps} of ${n_reps} replicates)"
  done

  # 5.2 Union peaks: merge consensus peaks of all conditions, remove IgG
  # regions and narrow peaks, add peak IDs
  for cond in "${conditions[@]}"; do
    cat "peak_consensus_${cond}.bed"
  done \
    | sort -k1,1 -k2,2n \
    | bedtools merge -d "${union_gap}" -i - > union.tmp.bed

  if [[ "${IGG_MODE}" == subtract && -n "$(samples_of "${IGG_MARK}")" ]]; then
    make_igg_regions IgG_regions.bed
    bedtools subtract -a union.tmp.bed -b IgG_regions.bed > union.tmp2.bed
    mv union.tmp2.bed union.tmp.bed
    log "${mark}: removed $(wc -l < IgG_regions.bed) IgG regions"
  fi

  local union="peak_${mark}_union.bed"
  awk -v min_width="${min_width}" -v mark="${mark}" 'BEGIN { OFS = "\t" }
      $3 - $2 >= min_width {
        n++
        print $1, $2, $3, mark "_peak" n, $3 - $2, "."
      }' union.tmp.bed > "${union}"
  rm union.tmp.bed
  [[ -s "${union}" ]] || die "${mark}: no union peaks; check peak calling"
  log "${mark}: $(wc -l < "${union}") union peaks"

  # SAF format for featureCounts: Start is 1-based
  awk 'BEGIN { OFS = "\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
       { print $4, $1, $2 + 1, $3, "." }' "${union}" > "peak_${mark}_union.saf"

  # Union peaks found in each condition (for the Venn diagram)
  for cond in "${conditions[@]}"; do
    bedtools intersect -u -a "${union}" -b "peak_consensus_${cond}.bed" \
      | cut -f 4 > "peak_${cond}_in_union.txt"
  done

  # 5.3 Plots: replicate Venn, condition Venn, peak widths
  activate_env R
  Rscript "${R_CODE_DIR}/peak_summary.R" --dir="${out_dir}" --mark="${mark}" \
    --conditions="$(IFS=,; echo "${conditions[*]}")" --min-reps="${MIN_REPS}"
}

main "$@"
