#!/bin/bash
#
# Makes a small simulated CUT&RUN data set to test the whole pipeline.
#
# Usage: make_test_data.sh <output_dir>
#
# Genome : chr1, chr2 (500 kb each, random sequence) and chrM
# Peaks  : 100 H3K27ac sites in Naive, Arm and Cl13 (2 replicates each)
#            sites 1-10      also in IgG       -> removed as IgG regions
#            every 7th site  4x lower in Naive -> Up vs Naive
#            every 10th site 4x higher in Cl13 -> Up in Cl13
# Reads  : paired-end 50 bp; fragments of ~150 bp (1 nucleosome) and ~300 bp
#          (2 nucleosomes); background reads; reads on chrM (to be removed)
# Lanes  : Naive_H3K27ac_R1 is split into 2 lanes (tests lane merging)

set -eo pipefail

#######################################
# Writes the genome, blacklist and genes. build_reference.sh uses these files
# instead of downloading them.
# Arguments:
#   Genome folder.
#######################################
make_reference() {
  local gnm_idx="$1"
  mkdir -p "${gnm_idx}"

  awk 'BEGIN {
      srand(1); split("A C G T", base, " ")
      n = split("chr1 500000 chr2 500000 chrM 16000", genome, " ")
      for (c = 1; c < n; c += 2) {
        print ">" genome[c]
        for (i = 0; i < genome[c + 1] / 1000; i++) {
          line = ""
          for (j = 0; j < 1000; j++) line = line base[int(rand() * 4) + 1]
          print line
        }
      }
    }' | gzip > "${gnm_idx}/mm10.fa.gz"

  # One blacklist region, between two peak sites
  printf 'chr1\t207500\t208500\n' | gzip > "${gnm_idx}/mm10-blacklist.v2.bed.gz"

  # A gene every 50 kb, with 2 exons
  awk 'BEGIN {
      OFS = "\t"
      for (c = 1; c <= 2; c++) for (s = 20000; s < 480000; s += 50000) {
        g++; strand = (g % 2 ? "+" : "-"); e = s + 20000
        id = "gene_id \"ENSMUSG" g ".1\"; gene_name \"Gene" g "\"; " \
          "gene_type \"protein_coding\";"
        tx = id " transcript_id \"ENSMUST" g ".1\";"
        print "chr" c, "TEST", "gene", s, e, ".", strand, ".", id
        print "chr" c, "TEST", "transcript", s, e, ".", strand, ".", tx
        print "chr" c, "TEST", "exon", s, s + 2000, ".", strand, ".", tx
        print "chr" c, "TEST", "exon", e - 3000, e, ".", strand, ".", tx
      }
    }' | gzip > "${gnm_idx}/genes.gtf.gz"
}

#######################################
# Writes the sample list and the comparison list.
# Arguments:
#   Output folder.
#######################################
make_lists() {
  local naive1="Naive_H3K27ac_R1_L1_R1.fastq.gz,Naive_H3K27ac_R1_L2_R1.fastq.gz"
  local naive2="Naive_H3K27ac_R1_L1_R2.fastq.gz,Naive_H3K27ac_R1_L2_R2.fastq.gz"
  {
    printf 'sample\tfastq1\tfastq2\tmark\tcondition\treplicate\tcontrol\t'
    printf 'batch\n'
    printf 'Naive_H3K27ac_R1\t%s\t%s\tH3K27ac\tNaive\t1\tIgG_R1\tb1\n' \
      "${naive1}" "${naive2}"
    local row
    for row in Naive_H3K27ac_R2:Naive:2:IgG_R2:b2 \
      Arm_H3K27ac_R1:Arm:1:IgG_R1:b1 Arm_H3K27ac_R2:Arm:2:IgG_R2:b2 \
      Cl13_H3K27ac_R1:Cl13:1:IgG_R1:b1 Cl13_H3K27ac_R2:Cl13:2:IgG_R2:b2 \
      IgG_R1:IgG:1:-:b1 IgG_R2:IgG:2:-:b2; do
      IFS=: read -r sample cond rep control batch <<< "${row}"
      local mark=H3K27ac
      if [[ "${cond}" == IgG ]]; then mark=IgG; fi
      printf '%s\t%s_R1.fastq.gz\t%s_R2.fastq.gz\t%s\t%s\t%s\t%s\t%s\n' \
        "${sample}" "${sample}" "${sample}" "${mark}" "${cond}" "${rep}" \
        "${control}" "${batch}"
    done
  } > "$1/sample_list.txt"

  printf 'mark\ttreatment\treference\n*\tArm\tNaive\n*\tCl13\tNaive\n' \
    > "$1/comparison_list.txt"
  printf '*\tCl13\tArm\n' >> "$1/comparison_list.txt"
}

#######################################
# Writes the paired-end reads of one sample.
# Arguments:
#   Genome fasta.gz, sample ID, random seed, output prefix.
#######################################
make_reads() {
  gzip -dc "$1" \
    | awk -v sample="$2" -v seed="$3" -v r1="$4_R1.fastq" -v r2="$4_R2.fastq" '
      function revcomp(s,   r, i) {
        r = ""
        for (i = length(s); i > 0; i--) r = r comp[substr(s, i, 1)]
        return r
      }
      function fragment(chr, start, len,   seq) {
        seq = substr(genome[chr], start, len); n++
        print "@" sample "_" n "/1\n" substr(seq, 1, 50) "\n+\n" qual > r1
        print "@" sample "_" n "/2\n" revcomp(substr(seq, len - 49, 50)) \
          "\n+\n" qual > r2
      }
      function fragment_length() {
        return rand() < 0.85 ? 120 + int(rand() * 60) : 280 + int(rand() * 60)
      }
      /^>/ { chr = substr($1, 2); next }
      { genome[chr] = genome[chr] $0 }
      END {
        srand(seed)
        comp["A"] = "T"; comp["C"] = "G"; comp["G"] = "C"; comp["T"] = "A"
        qual = sprintf("%50s", ""); gsub(/ /, "I", qual)

        n_sites = (sample ~ /IgG/ ? 10 : 100)
        for (k = 1; k <= n_sites; k++) {
          depth = 40
          if (sample ~ /Naive/ && k % 7 == 0) depth = 10
          if (sample ~ /Cl13/ && k % 10 == 0) depth = 160
          chr = (k % 2 ? "chr1" : "chr2"); site = 5000 + int((k - 1) / 2) * 9000
          for (j = 0; j < depth; j++) {
            fragment(chr, site + int(rand() * 800), fragment_length())
          }
        }
        for (j = 0; j < 3000; j++) {
          start = 1 + int(rand() * 498000)
          fragment("chr" (1 + j % 2), start, fragment_length())
        }
        for (j = 0; j < 300; j++) {
          fragment("chrM", 1 + int(rand() * 15500), fragment_length())
        }
      }'
  gzip -f "$4_R1.fastq" "$4_R2.fastq"
}

#######################################
# Splits the reads of one sample into 2 lanes.
# Arguments:
#   Raw folder, sample ID.
#######################################
split_lanes() {
  local read
  for read in R1 R2; do
    gzip -dc "$1/$2_${read}.fastq.gz" > "$1/lanes.tmp"
    head -n 8000 "$1/lanes.tmp" | gzip > "$1/$2_L1_${read}.fastq.gz"
    tail -n +8001 "$1/lanes.tmp" | gzip > "$1/$2_L2_${read}.fastq.gz"
    rm "$1/lanes.tmp" "$1/$2_${read}.fastq.gz"
  done
}

main() {
  local out_dir="$1"
  if [[ -z "${out_dir}" ]]; then
    echo "Usage: make_test_data.sh <output_dir>" >&2
    exit 1
  fi
  local gnm_idx="${out_dir}/indexes/mm10"
  mkdir -p "${out_dir}/Raw"

  make_reference "${gnm_idx}"
  make_lists "${out_dir}"

  local seed=0 sample
  for sample in Naive_H3K27ac_R1 Naive_H3K27ac_R2 Arm_H3K27ac_R1 \
    Arm_H3K27ac_R2 Cl13_H3K27ac_R1 Cl13_H3K27ac_R2 IgG_R1 IgG_R2; do
    seed=$(( seed + 1 ))
    make_reads "${gnm_idx}/mm10.fa.gz" "${sample}" "${seed}" \
      "${out_dir}/Raw/${sample}"
  done
  split_lanes "${out_dir}/Raw" Naive_H3K27ac_R1
  echo "Test data ready in ${out_dir}"
}

main "$@"
