#!/bin/bash
#
# Runs the pipeline in the Docker image (one computer, no cluster). The
# folders of config.sh are mounted at the same paths, so the config is used
# as it is; all jobs run in the container (SCHEDULER=local).
#
# Usage: run_docker.sh <config.sh> [qc|align|peaks|matrix|diff]
#        run_docker.sh --status <config.sh>
# Image: CUTRUN_IMAGE (default cutrunpipeline), built with
#        docker build --platform linux/amd64 -t cutrunpipeline .

set -eo pipefail

readonly IMAGE="${CUTRUN_IMAGE:-cutrunpipeline}"
readonly IMAGE_DIR="/opt/CUTRUNpipeline"
HOST_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly HOST_PIPELINE_DIR

mounts=()

#######################################
# Adds a folder to the mounts, at the same path (once).
# Globals:
#   mounts
# Arguments:
#   Folder.
#######################################
add_mount() {
  local dir="$1"
  if [[ -z "${dir}" ]]; then return 0; fi
  [[ -d "${dir}" ]] || { echo "Folder not found: ${dir}" >&2; exit 1; }
  dir="$(cd "${dir}" && pwd)"
  if [[ " ${mounts[*]} " == *" -v ${dir}:${dir} "* ]]; then return 0; fi
  mounts+=(-v "${dir}:${dir}")
}

main() {
  local script=run_pipeline.sh
  if [[ "$1" == --status ]]; then
    script=check_status.sh
    shift
  fi
  if [[ ! -f "$1" ]]; then
    echo "Usage: run_docker.sh [--status] <config.sh> [part]" >&2
    exit 1
  fi
  local config
  config="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"

  # Folders named in the config (read here, on this computer)
  local folders
  folders="$(PIPELINE_DIR="${HOST_PIPELINE_DIR}" bash -c '
      source "$1"
      mkdir -p "${WORKDIR}"
      printf "%s\n" "${WORKDIR}" "${DATADIR}" "${LIBRARIES}" \
        "$(dirname "${SAMPLE_LIST}")" "$(dirname "${PEAK_SETTINGS}")"
      if [[ -n "${COMPARISON_LIST}" ]]; then dirname "${COMPARISON_LIST}"; fi
      if [[ -n "${SPIKEIN_IDX}" ]]; then dirname "${SPIKEIN_IDX}"; fi' \
    _ "${config}")"
  add_mount "$(dirname "${config}")"
  local dir
  while IFS= read -r dir; do add_mount "${dir}"; done <<< "${folders}"

  # Files are owned by you, not by root
  docker run --rm --platform linux/amd64 \
    --user "$(id -u):$(id -g)" -e HOME=/tmp \
    "${mounts[@]}" "${IMAGE}" \
    bash "${IMAGE_DIR}/${script}" "${config}" "${@:2}"
}

main "$@"
