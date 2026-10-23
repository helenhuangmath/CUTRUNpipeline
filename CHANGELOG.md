# Changelog

All notable changes are listed here. Version numbers follow [Semantic Versioning](https://semver.org):
MAJOR (results or inputs change in a way that is not backward compatible), MINOR (new features),
PATCH (bug fixes).

## [Unreleased]

### Added
- Docker image with all software (`Dockerfile`) and `docker/run_docker.sh` to run the pipeline on
  one computer with the same `config.sh`; tested in CI.

### Fixed
- A failed command inside `$(...)` did not stop the script (bash `inherit_errexit`): e.g. after
  MACS2 failed, the peak step went on with a missing peak file and hid the real error.
- When MACS2 or FlexPeak fails, the end of its log is shown in the job log.
- CI: ShellCheck settings (`.shellcheckrc`) were missing, so the lint job failed on notes about
  sourced files and settings from `config.sh`.
- CI: MACS2 logs are kept when the tests fail; actions updated to Node 24; faster lintr install.

## [0.1.0] - 2026-10-23

First version, rebuilt from the lab's CUT&RUN scripts (2020-2026) into one pipeline.

### Added
- One command for the whole analysis (`run_pipeline.sh`), driven by `config.sh`, `sample_list.txt`
  and `comparison_list.txt`; jobs for LSF, SLURM or a single computer, with job dependencies.
- Input checks before any job is submitted (`scripts/validate.sh`).
- Trimming with fastp; FastQC and MultiQC reports.
- Rerun after a failure: finished steps are skipped (`Status/*.done`); `check_status.sh`.
- MACS2 settings of each histone mark in one table (`config/peak_settings.txt`).
- Consensus peaks for any number of replicates; IgG regions removed from union peaks.
- Generic DESeq2 analysis for any conditions and comparisons, with batch correction.
- QC summary of all samples (read pairs, alignment, duplicates, peaks, FRiP, fragment sizes).
- Reproducibility records: settings of every run, software versions, R sessionInfo.
- Reference builder for mm10, hg38, hg19 and E. coli spike-in (`reference/build_reference.sh`).
- Conda environments (`envs/`), automated tests on simulated data (`tests/`), GitHub Actions CI.
- FlexPeak as a second peak caller next to MACS2 (`PEAK_CALLERS`, `ANALYSIS_PEAKS`), with a
  per-sample comparison of the two callers.
- `qc` part: FastQC and MultiQC of the raw fastq files; FastQC of the trimmed reads.
- Separate count matrix step (`07_count_matrix.sh`): raw, DESeq2-normalized and VST matrices of
  all samples; the DESeq2 step reuses its size factors.
- Nature-style plots (7 pt, journal column widths, colour-blind-safe colours) and a multi-panel
  summary figure of each mark with signal profiles and heatmaps (`09_figures.sh`).
- Searchable HTML report of all marks (`10_report.sh`): find a gene and see its peaks, statistics,
  signal profile and counts.
- Code style: Google Shell and R style guides, checked with ShellCheck and lintr in CI.
- MIT license; contributor guide, code of conduct, contributors list, issue and pull request templates.
