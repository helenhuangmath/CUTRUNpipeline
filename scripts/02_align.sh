#!/bin/bash
#
# Step 2: alignment of one sample: bowtie2 -> filter -> remove duplicates.
#
# Usage: 02_align.sh <config.sh> <sample_id>
# Output: BAM/<sample_id>.map.dedup.bam (+ .bai), .bowtie2.log,
#         .align_stats.txt, .fragment_size.txt

source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

main() {
  read_config "$1"
  local sample_id="$2"
  activate_env tools
  cd "${WORKDIR}" || exit 1
  mkdir -p BAM

  local trim1="Fastq/${sample_id}_R1.trim.fastq.gz"
  local trim2="Fastq/${sample_id}_R2.trim.fastq.gz"
  local sort_bam="BAM/${sample_id}.map.sort.bam"
  local dedup_bam="BAM/${sample_id}.map.dedup.bam"
  local bowtie2_args
  read -ra bowtie2_args <<< "${BOWTIE2_ARGS}"

  # 2.1 Align (Henikoff lab settings) and keep good read pairs
  #   --local --very-sensitive-local : local alignment, soft-clips adapters
  #   --no-mixed --no-discordant     : both reads of a pair align together
  #   -I 10 -X 700                   : fragments of 10 to 700 bp
  #   -F 1804 : remove unmapped, mate unmapped, secondary, QC-fail, duplicate
  #   -f 2    : keep properly paired reads
  #   -q      : minimum mapping quality, removes multi-mapping reads
  #   -L      : keep chr1 to chrY (removes chrM, random and chrUn)
  #   fixmate : flags reads whose mate was removed; the 2nd filter drops them
  bowtie2 "${bowtie2_args[@]}" -p "${THREADS}" -x "${BOWTIE_IDX}" \
    -1 "${trim1}" -2 "${trim2}" 2> "BAM/${sample_id}.bowtie2.log" \
    | samtools view -u -F 1804 -f 2 -q "${MIN_MAPQ}" -L "${CHR1TOY_BED}" - \
    | samtools fixmate -m - - \
    | samtools view -u -F 1804 -f 2 - \
    | samtools sort -@ "${THREADS}" -m 1G -T "BAM/${sample_id}.tmp" \
      -o "${sort_bam}" -

  # 2.2 Remove PCR duplicates (REMOVE_DUPLICATES=false: only mark them)
  local dedup_args=(-s -@ "${THREADS}")
  if [[ "${REMOVE_DUPLICATES}" == true ]]; then dedup_args+=(-r); fi
  samtools markdup "${dedup_args[@]}" "${sort_bam}" "${dedup_bam}" \
    2> "BAM/${sample_id}_dedup.log"
  samtools index "${dedup_bam}"
  samtools flagstat "${dedup_bam}" > "${dedup_bam}_flagstat"

  # 2.3 Spike-in (optional): number of fragments aligned to e.g. E. coli
  local spikein_pairs="NA"
  if [[ -n "${SPIKEIN_IDX}" ]]; then
    spikein_pairs="$(bowtie2 "${bowtie2_args[@]}" --no-overlap --no-dovetail \
      -p "${THREADS}" -x "${SPIKEIN_IDX}" -1 "${trim1}" -2 "${trim2}" \
      2> "BAM/${sample_id}.spikein.bowtie2.log" \
      | samtools view -c -f 66 -F 1804 -q "${MIN_MAPQ}" -)"
  fi

  # 2.4 Fragment sizes (insert size of read 1): columns size, count
  samtools view -f 66 "${dedup_bam}" \
    | awk '{ s = ($9 < 0 ? -$9 : $9); n[s]++ }
           END { for (s in n) print s "\t" n[s] }' \
    | sort -k1,1n > "BAM/${sample_id}.fragment_size.txt"

  # 2.5 Read pairs at each step (read 1 only, -f 64)
  #   input_pairs    : read pairs given to bowtie2
  #   filtered_pairs : read pairs after 2.1
  #   unique_pairs   : read pairs that are not PCR duplicates
  local input_pairs align_rate filtered_pairs unique_pairs dup_rate
  input_pairs="$(awk '/reads; of these:/ { print $1 }' \
    "BAM/${sample_id}.bowtie2.log")"
  align_rate="$(awk '/overall alignment rate/ { sub("%", "", $1); print $1 }' \
    "BAM/${sample_id}.bowtie2.log")"
  filtered_pairs="$(samtools view -c -f 64 "${sort_bam}")"
  unique_pairs="$(samtools view -c -f 64 -F 1024 "${dedup_bam}")"
  dup_rate="$(awk -v u="${unique_pairs}" -v f="${filtered_pairs}" \
    'BEGIN { printf "%.4f", 1 - u / f }')"

  {
    printf 'sample\tinput_pairs\talign_rate\tfiltered_pairs\tunique_pairs'
    printf '\tdup_rate\tspikein_pairs\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${sample_id}" "${input_pairs}" \
      "${align_rate}" "${filtered_pairs}" "${unique_pairs}" "${dup_rate}" \
      "${spikein_pairs}"
  } > "BAM/${sample_id}.align_stats.txt"
  log "${sample_id}: ${input_pairs} read pairs, ${align_rate}% aligned," \
    "${unique_pairs} unique pairs (duplicate rate ${dup_rate})"

  # Remove intermediate files to save disk space
  rm -f "${sort_bam}" "${trim1}" "${trim2}" \
    "Fastq/${sample_id}_R1.merged.fastq.gz" \
    "Fastq/${sample_id}_R2.merged.fastq.gz"
}

main "$@"
