#!/bin/bash
#
# Shared constants and functions of the CUT&RUN pipeline.
# This file is a library: every pipeline script sources it. Do not run it.

set -Eeo pipefail
trap 'err "failed in $(basename "$0") line ${LINENO}: ${BASH_COMMAND}"' ERR

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly PIPELINE_DIR
readonly CODE_DIR="${PIPELINE_DIR}/scripts"
# shellcheck disable=SC2034  # used by the scripts that source this file
readonly R_CODE_DIR="${PIPELINE_DIR}/scripts/R"

#######################################
# Prints a message with a time stamp to STDERR.
# Arguments:
#   Message.
#######################################
log() {
  echo "[$(date '+%F %T')] $*" >&2
}

#######################################
# Prints an error message with a time stamp to STDERR.
# Arguments:
#   Message.
#######################################
err() {
  echo "[$(date '+%F %T')] ERROR: $*" >&2
}

#######################################
# Prints an error message and exits with status 1.
# Arguments:
#   Message.
#######################################
die() {
  err "$@"
  exit 1
}

#######################################
# Reads config.sh and sets default values for settings not in it.
# Globals:
#   CONFIG and all settings of config.sh (written).
# Arguments:
#   Path of config.sh.
#######################################
read_config() {
  [[ -f "$1" ]] || die "Config file not found: $1"
  CONFIG="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  # shellcheck source=../config/config.sh
  source "${CONFIG}"

  # Defaults
  : "${PEAK_SETTINGS:=${PIPELINE_DIR}/config/peak_settings.txt}"
  : "${SCHEDULER:=lsf}"
  : "${THREADS:=8}"
  : "${MEMORY:=32000}"
  : "${RUN_FASTQC:=true}"
  : "${RUN_TRIMMING:=true}"
  : "${BOWTIE2_ARGS:=--local --very-sensitive-local --no-unal --no-mixed \
--no-discordant --phred33 -I 10 -X 700}"
  : "${MIN_MAPQ:=10}"
  : "${REMOVE_DUPLICATES:=true}"
  : "${SPIKEIN_SCALE:=10000}"
  : "${BW_NORM:=CPM}"
  : "${PEAK_CALLERS:=macs2}"
  : "${ANALYSIS_PEAKS:=macs2}"
  : "${MACS2_FORMAT:=BAMPE}"
  : "${IGG_MARK:=IgG}"
  : "${IGG_MODE:=subtract}"
  : "${IGG_EXTEND:=200}"
  : "${MIN_REPS:=2}"
  : "${TSS_UP:=2000}"
  : "${TSS_DOWN:=500}"
  : "${FLANK:=10000}"
  : "${DESIGN:=~ condition}"
  : "${MIN_COUNT:=15}"
  : "${PADJ:=0.05}"
  : "${LFC:=0.585}"
  : "${FORCE:=false}"

  mkdir -p "${WORKDIR}/Log" "${WORKDIR}/Status"
}

#######################################
# Activates the software of the pipeline: a conda environment or modules.
# Globals:
#   CONDA_ENV, CONDA_ENV_R, MODULES, MODULES_R
# Arguments:
#   "tools" or "R".
#######################################
activate_env() {
  local env_name="${CONDA_ENV}"
  local modules="${MODULES}"
  if [[ "$1" == R ]]; then
    env_name="${CONDA_ENV_R}"
    modules="${MODULES_R}"
  fi

  if [[ -n "${modules}" ]]; then
    type module &> /dev/null || source /etc/profile.d/modules.sh
    local module_list
    read -ra module_list <<< "${modules}"
    module load "${module_list[@]}"
  fi

  if [[ -n "${env_name}" ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_name}"
  fi
}

#######################################
# Prints a table without comment lines, empty lines and Windows line ends.
# Arguments:
#   Table file.
# Outputs:
#   The table to STDOUT.
#######################################
read_table() {
  sed -e 's/\r$//' -e '/^#/d' -e '/^[[:space:]]*$/d' "$1"
}

#######################################
# Prints all sample IDs of the sample list.
# Globals:
#   SAMPLE_LIST
# Outputs:
#   One sample ID per line.
#######################################
all_samples() {
  read_table "${SAMPLE_LIST}" | awk -F'\t' 'NR > 1 { print $1 }'
}

#######################################
# Prints one value of a sample from the sample list.
# Globals:
#   SAMPLE_LIST
# Arguments:
#   Sample ID, column name. Example: sample_info Naive_R1 condition
# Outputs:
#   The value, e.g. "Naive".
#######################################
sample_info() {
  read_table "${SAMPLE_LIST}" \
    | awk -F'\t' -v id="$1" -v column="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $1 == id { print $col[column] }'
}

#######################################
# Prints the samples of one mark, and optionally of one condition.
# Globals:
#   SAMPLE_LIST
# Arguments:
#   Mark, condition (optional).
# Outputs:
#   One sample ID per line.
#######################################
samples_of() {
  read_table "${SAMPLE_LIST}" \
    | awk -F'\t' -v mark="$1" -v cond="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["mark"] == mark && (cond == "" || $col["condition"] == cond) {
          print $1
        }'
}

#######################################
# Prints the conditions of one mark, in the order of the sample list.
# Globals:
#   SAMPLE_LIST
# Arguments:
#   Mark.
# Outputs:
#   One condition per line.
#######################################
conditions_of() {
  read_table "${SAMPLE_LIST}" \
    | awk -F'\t' -v mark="$1" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["mark"] == mark && !seen[$col["condition"]]++ {
          print $col["condition"]
        }'
}

#######################################
# Prints all marks except IgG.
# Globals:
#   SAMPLE_LIST, IGG_MARK
# Outputs:
#   One mark per line.
#######################################
all_marks() {
  read_table "${SAMPLE_LIST}" \
    | awk -F'\t' -v igg="${IGG_MARK}" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $col["mark"] != igg && !seen[$col["mark"]]++ { print $col["mark"] }'
}

#######################################
# Prints a peak setting of a mark; marks not listed use the "default" row.
# Globals:
#   PEAK_SETTINGS
# Arguments:
#   Mark, column name. Example: peak_setting H3K27ac macs2_options
# Outputs:
#   The value.
#######################################
peak_setting() {
  read_table "${PEAK_SETTINGS}" \
    | awk -F'\t' -v mark="$1" -v column="$2" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        $1 == mark { value = $col[column]; found = 1 }
        $1 == "default" { default_value = $col[column] }
        END { print (found ? value : default_value) }'
}

#######################################
# Prints the full path of a fastq file (relative paths are in DATADIR).
# Globals:
#   DATADIR
# Arguments:
#   File name or path.
#######################################
full_path() {
  if [[ "$1" == /* ]]; then
    echo "$1"
  else
    echo "${DATADIR}/$1"
  fi
}

#######################################
# Prints the raw fastq file of read 1 or 2 of a sample. Several lanes
# (comma-separated in the sample list) are first merged into
# Fastq/<sample>_R<read>.merged.fastq.gz.
# Globals:
#   WORKDIR, SAMPLE_LIST, DATADIR
# Arguments:
#   Sample ID, read (1 or 2).
# Outputs:
#   Path of the fastq file.
#######################################
merge_lanes() {
  local files
  files="$(sample_info "$1" "fastq$2")"
  if [[ "${files}" != *,* ]]; then
    full_path "${files}"
    return 0
  fi

  local merged="${WORKDIR}/Fastq/$1_R$2.merged.fastq.gz"
  local lanes lane
  IFS=',' read -ra lanes <<< "${files}"
  : > "${merged}"
  for lane in "${lanes[@]}"; do
    # gzip files can be joined by simple concatenation
    cat "$(full_path "${lane}")" >> "${merged}"
  done
  echo "${merged}"
}

#######################################
# Runs a pipeline step once: a finished step leaves Status/<step>.done and is
# skipped the next time (unless FORCE=true).
# Globals:
#   WORKDIR, FORCE
# Arguments:
#   Step name, script in scripts/, arguments of the script.
#   Example: run_step Naive_R1.align 02_align.sh config.sh Naive_R1
#######################################
run_step() {
  local step="$1"
  shift
  if [[ -f "${WORKDIR}/Status/${step}.done" && "${FORCE}" != true ]]; then
    log "${step}: already done, skip"
    return 0
  fi
  log "${step}: start"
  bash "${CODE_DIR}/$1" "${@:2}"
  touch "${WORKDIR}/Status/${step}.done"
  log "${step}: done"
}

#######################################
# Submits a pipeline script as a cluster job, or runs it here when
# SCHEDULER=local.
# Globals:
#   WORKDIR, SCHEDULER, THREADS, MEMORY, QUEUE, EXTRA_SUBMIT_ARGS
# Arguments:
#   Job name, job IDs to wait for (space-separated, may be empty),
#   script in scripts/, arguments of the script.
# Outputs:
#   The job ID ("local" when run here).
#######################################
submit_job() {
  local job_name="$1"
  local wait_ids
  read -ra wait_ids <<< "$2"
  local cmd=(bash "${CODE_DIR}/$3" "${@:4}")
  local out_log="${WORKDIR}/Log/${job_name}.out"
  local err_log="${WORKDIR}/Log/${job_name}.err"
  local extra_args
  read -ra extra_args <<< "${EXTRA_SUBMIT_ARGS}"

  case "${SCHEDULER}" in
    lsf)
      local args=(-J "${job_name}" -n "${THREADS}" -M "${MEMORY}"
        -R "rusage[mem=${MEMORY}] span[hosts=1]"
        -o "${out_log}" -e "${err_log}")
      if [[ -n "${QUEUE}" ]]; then args+=(-q "${QUEUE}"); fi
      if (( ${#wait_ids[@]} > 0 )); then
        local dependency
        dependency="$(printf 'done(%s) && ' "${wait_ids[@]}")"
        args+=(-w "${dependency% && }")
      fi
      bsub "${args[@]}" "${extra_args[@]}" "${cmd[@]}" \
        | sed -n 's/^Job <\([0-9]*\)>.*/\1/p'
      ;;
    slurm)
      local args=(--parsable -J "${job_name}" -c "${THREADS}"
        --mem="${MEMORY}M" -o "${out_log}" -e "${err_log}")
      if [[ -n "${QUEUE}" ]]; then args+=(-p "${QUEUE}"); fi
      if (( ${#wait_ids[@]} > 0 )); then
        args+=(--dependency="afterok:$(IFS=:; echo "${wait_ids[*]}")"
          --kill-on-invalid-dep=yes)
      fi
      sbatch "${args[@]}" "${extra_args[@]}" \
        --wrap "$(printf '%q ' "${cmd[@]}")" \
        | cut -d ';' -f 1
      ;;
    local)
      log "Running ${job_name} (log: ${err_log})"
      "${cmd[@]}" < /dev/null > "${out_log}" 2> "${err_log}" \
        || die "${job_name} failed, see ${err_log}"
      echo local
      ;;
    *)
      die "SCHEDULER must be lsf, slurm or local"
      ;;
  esac
}
