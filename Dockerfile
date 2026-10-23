# All software of the CUT&RUN pipeline in one image, for running it on one
# computer (laptop, workstation, cloud VM). On an HPC cluster use the conda
# environments instead (README section 1).
#
# Build:  docker build --platform linux/amd64 -t cutrunpipeline .
# Run:    bash docker/run_docker.sh config.sh       (README section 1)
#
# The image has the two conda environments of envs/*.yml (cutrun, cutrun_r),
# so the default CONDA_ENV / CONDA_ENV_R settings of config.sh work unchanged.
# linux/amd64: some bioconda tools have no ARM build; Apple Silicon runs the
# image through emulation.

FROM condaforge/miniforge3:26.7.2-0

COPY envs/ /tmp/envs/
RUN conda env create -f /tmp/envs/cutrun.yml \
    && conda env create -f /tmp/envs/cutrun_r.yml \
    && conda clean --all --yes \
    && rm -rf /tmp/envs

# Optional second peak caller, e.g.
#   docker build --build-arg FLEXPEAK="flexpeak[bam] @ git+<URL>" ...
ARG FLEXPEAK=""
RUN if [ -n "${FLEXPEAK}" ]; then \
      conda run -n cutrun pip install --no-cache-dir "${FLEXPEAK}"; \
    fi

COPY . /opt/CUTRUNpipeline
WORKDIR /opt/CUTRUNpipeline

# Jobs run inside the container, not on a cluster scheduler
ENV CUTRUN_SCHEDULER=local
# Environments used by tests/run_tests.sh, the default command
ENV TEST_CONDA_ENV=cutrun TEST_CONDA_ENV_R=cutrun_r

CMD ["bash", "/opt/CUTRUNpipeline/tests/run_tests.sh", "/tmp/cutrun_test"]
