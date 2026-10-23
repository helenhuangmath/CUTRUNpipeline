#!/bin/bash
#
# Prepares the reference files of one genome (run once per genome).
#
# Usage: build_reference.sh <species> <libraries_dir> [threads]
#   species       : mm10 | hg38 | hg19 | ecoli (spike-in, bowtie2 index only)
#   libraries_dir : folder of all genomes; files go to <libraries_dir>/<species>
# The bowtie2 index takes 1-3 hours, so submit it as a job, e.g.
#   bsub -n 8 -M 32000 -o ref.log bash build_reference.sh mm10 ~/indexes 8
# Files already in the folder are not downloaded again.

source "$(dirname "${BASH_SOURCE[0]}")/../scripts/common.sh"

readonly UCSC="https://hgdownload.soe.ucsc.edu/goldenPath"
readonly BLACKLISTS="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/\
master/lists"
readonly GENCODE="https://ftp.ebi.ac.uk/pub/databases/gencode"
readonly GENCODE_MOUSE="${GENCODE}/Gencode_mouse/release_M25"
readonly GENCODE_HUMAN="${GENCODE}/Gencode_human/release_44"
readonly NCBI_ECOLI="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/\
GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"

#######################################
# Downloads a file unless it is already there.
# Arguments:
#   URL, output file.
#######################################
download() {
  if [[ ! -s "$2" ]]; then
    log "Downloading $1"
    curl -fL --retry 3 -o "$2" "$1"
  fi
}

main() {
  local species="$1"
  local libraries_dir="$2"
  local threads="${3:-8}"
  if [[ -z "${species}" || -z "${libraries_dir}" ]]; then
    sed -n '3,11p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
  fi

  local genome_url blacklist_url gtf_url refsize
  case "${species}" in
    mm10)
      genome_url="${UCSC}/mm10/bigZips/mm10.fa.gz"
      blacklist_url="${BLACKLISTS}/mm10-blacklist.v2.bed.gz"
      gtf_url="${GENCODE_MOUSE}/gencode.vM25.annotation.gtf.gz"
      refsize=1.87e9
      ;;
    hg38)
      genome_url="${UCSC}/hg38/bigZips/hg38.fa.gz"
      blacklist_url="${BLACKLISTS}/hg38-blacklist.v2.bed.gz"
      gtf_url="${GENCODE_HUMAN}/gencode.v44.annotation.gtf.gz"
      refsize=2.7e9
      ;;
    hg19)
      genome_url="${UCSC}/hg19/bigZips/hg19.fa.gz"
      blacklist_url="${BLACKLISTS}/hg19-blacklist.v2.bed.gz"
      gtf_url="${GENCODE_HUMAN}/GRCh37_mapping"
      gtf_url+="/gencode.v44lift37.annotation.gtf.gz"
      refsize=2.7e9
      ;;
    ecoli)
      genome_url="${NCBI_ECOLI}"
      ;;
    *)
      die "species must be mm10, hg38, hg19 or ecoli"
      ;;
  esac

  # Conda environments of the pipeline (set CONDA_ENV= to use tools on PATH)
  CONDA_ENV="${CONDA_ENV-cutrun}"
  CONDA_ENV_R="${CONDA_ENV_R-cutrun_r}"
  activate_env tools

  local gnm_idx="${libraries_dir}/${species}"
  mkdir -p "${gnm_idx}/bowtie2"
  cd "${gnm_idx}" || exit 1

  # 1. Genome, chromosome sizes, bowtie2 index
  download "${genome_url}" "${species}.fa.gz"
  if [[ ! -s "${species}.fa" ]]; then
    gzip -dc "${species}.fa.gz" > "${species}.fa"
  fi
  samtools faidx "${species}.fa"
  cut -f 1,2 "${species}.fa.fai" > "${species}.genomesize"
  if [[ ! -f "bowtie2/${species}.1.bt2" ]]; then
    log "Building the bowtie2 index"
    bowtie2-build --threads "${threads}" "${species}.fa" "bowtie2/${species}" \
      > bowtie2/build.log 2>&1
  fi
  if [[ "${species}" == ecoli ]]; then
    log "Done. In config.sh set: SPIKEIN_IDX=\"${gnm_idx}/bowtie2/ecoli\""
    return 0
  fi

  # 2. chr1 to chrY (mouse chr1-19, X, Y; human chr1-22, X, Y)
  awk 'BEGIN { OFS = "\t" } $1 ~ /^chr([0-9]+|X|Y)$/ { print $1, 0, $2 }' \
    "${species}.genomesize" | sort -k1,1 -k2,2n > "${species}.chr1toY.bed"

  # 3. Blacklist (ENCODE v2): remove zero-length regions, sort, merge
  download "${blacklist_url}" "${species}-blacklist.v2.bed.gz"
  gzip -dc "${species}-blacklist.v2.bed.gz" \
    | awk 'BEGIN { OFS = "\t" } $3 > $2 { print $1, $2, $3 }' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "${species}-blacklist.clean.bed"

  # 4. Genes (GENCODE): TxDb for ChIPseeker and gene table
  download "${gtf_url}" genes.gtf.gz
  activate_env R
  Rscript "${R_CODE_DIR}/build_txdb.R" --gtf=genes.gtf.gz \
    --out-txdb=genes.txdb.sqlite --out-genes=genes.txt

  log "Reference ready in ${gnm_idx}. In config.sh set:" \
    "SPECIES=\"${species}\", LIBRARIES=\"${libraries_dir}\"," \
    "REFSIZE=\"${refsize}\""
}

main "$@"
