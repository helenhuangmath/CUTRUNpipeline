#!/bin/bash
#
# Step 7: read count matrix of all samples of one mark in its union peaks.
# Counts fragments of every sample (featureCounts), then writes the raw and
# the normalized count matrix, size factors, FRiP and a sample PCA.
#
# Usage: 07_count_matrix.sh <config.sh> <mark>
# Output: ReadCount/<mark>/RawRC_<mark>.txt   raw counts, peaks x samples
#         ReadCount/<mark>/NormRC_<mark>.txt  DESeq2-normalized counts
#         ReadCount/<mark>/sizeFactors_<mark>.txt, plot_PCA_<mark>.pdf, ...

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local mark="$2"
  activate_env tools

  local out_dir="${WORKDIR}/ReadCount/${mark}"
  mkdir -p "${out_dir}"

  local bams=()
  local sample
  while read -r sample; do
    bams+=("${WORKDIR}/BAM/${sample}.map.dedup.bam")
  done < <(samples_of "${mark}")

  # 7.1 Count fragments (read pairs) of each sample in each union peak
  #   -p --countReadPairs : count fragments, not reads (subread >= 2.0.2)
  featureCounts -p --countReadPairs -F SAF -T "${THREADS}" \
    -a "${WORKDIR}/Consensus/${mark}/peak_${mark}_union.saf" \
    -o "${out_dir}/featureCounts_${mark}.txt" "${bams[@]}" \
    2> "${out_dir}/featureCounts_${mark}.log"

  # 7.2 Raw and normalized count matrices, FRiP, PCA
  activate_env R
  Rscript "${R_CODE_DIR}/count_matrix.R" \
    --mark="${mark}" \
    --counts="${out_dir}/featureCounts_${mark}.txt" \
    --anno="${WORKDIR}/Annotation/${mark}/peak_${mark}_union.anno.txt" \
    --sample-list="${SAMPLE_LIST}" \
    --outdir="${out_dir}"
}

main "$@"
