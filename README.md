# CUTRUNpipeline

A CUT&RUN analysis pipeline for histone marks and chromatin proteins, from paired-end FASTQ files to peaks,
gene annotation, count matrices, differential peaks, publication-style figures and a searchable HTML report.
Written in plain bash and R; runs on an HPC cluster (LSF or SLURM) or on one computer.

Fill in **a sample list** and **a config file**, and run **one command**.

```
FASTQ ─▶ FastQC + MultiQC ─▶ fastp ─▶ bowtie2 ─▶ filter + remove duplicates ─▶ BAM ─▶ bigwig
                                                                                 │
                                   peaks of each sample: MACS2 and/or FlexPeak (settings per mark)
                                                                                 │
                        consensus peaks of each condition ─▶ union peaks of each mark (IgG removed)
                                                                                 │
                 annotation ─▶ raw / normalized / VST count matrices of all samples ─▶ QC summary
                                                                                 │
              DESeq2 differential peaks ─▶ Nature-style figures ─▶ HTML report (search a gene)
```

Each step, what it does and what to check: [PIPELINE_STEPS.txt](PIPELINE_STEPS.txt).

---

## 1. Install (once)

You need conda (Miniconda or Miniforge). Create the two software environments:

```bash
git clone <repository URL> ~/CUTRUNpipeline
cd ~/CUTRUNpipeline
conda env create -f envs/cutrun.yml      # tools: fastp, bowtie2, samtools, MACS2, deepTools, ...
conda env create -f envs/cutrun_r.yml    # R: DESeq2, ChIPseeker, ggplot2, patchwork, ...

# FlexPeak (optional second peak caller), into the tools environment:
conda activate cutrun
pip install "/path/to/FlexPeak[bam,fast]"
```

Software versions are fixed in `envs/*.yml`, so every user gets the same versions. If your cluster uses
modules instead, leave `CONDA_ENV` empty in the config and list the modules in `MODULES`.

## 2. Reference genome (once per genome)

```bash
bsub -n 8 -M 32000 -o ref_mm10.log bash ~/CUTRUNpipeline/reference/build_reference.sh mm10 ~/indexes 8
```

Supported: `mm10`, `hg38`, `hg19`, and `ecoli` (spike-in). It downloads the genome (UCSC), the blacklist
(ENCODE v2) and genes (GENCODE), and makes the bowtie2 index (1-3 hours), `mm10.genomesize`,
`mm10.chr1toY.bed`, `mm10-blacklist.clean.bed` and the gene annotation in `~/indexes/mm10/`.

## 3. Set up a project

```bash
mkdir -p ~/MyProject && cd ~/MyProject
cp ~/CUTRUNpipeline/config/{config.sh,sample_list.txt,comparison_list.txt} .
```

**sample_list.txt**: one row per sample, tab-separated (in Excel: save as *Tab delimited text*).

```
sample            fastq1                        fastq2                        mark     condition  replicate  control  batch
Naive_H3K27ac_R1  Naive_H3K27ac_R1_R1.fastq.gz  Naive_H3K27ac_R1_R2.fastq.gz  H3K27ac  Naive      1          IgG_R1   b1
Arm_H3K27ac_R1    Arm_H3K27ac_R1_R1.fastq.gz    Arm_H3K27ac_R1_R2.fastq.gz    H3K27ac  Arm        1          IgG_R1   b1
IgG_R1            IgG_R1_R1.fastq.gz            IgG_R1_R2.fastq.gz            IgG      IgG        1          -        b1
```

- `fastq1/2` are inside `DATADIR`; several lanes: `L1_R1.fastq.gz,L2_R1.fastq.gz`
- `mark` selects the peak settings in [config/peak_settings.txt](config/peak_settings.txt); IgG samples: `IgG`
- several marks can be in one sample list; each mark is analysed separately

**comparison_list.txt**: comparisons for DESeq2 (`*` = every mark; log2FC > 0 = higher in treatment).

```
mark  treatment  reference
*     Arm        Naive
*     Cl13       Naive
```

**config.sh**: fill in sections 1-4 (project folders, reference, software, cluster). The other settings
have defaults; all are explained in the file. Peak calling:

| setting | meaning |
|---|---|
| `PEAK_CALLERS="macs2 flexpeak"` | run MACS2 and FlexPeak on every sample (or only one of them) |
| `ANALYSIS_PEAKS="macs2"` | peaks used for consensus, counts and DESeq2: `macs2`, `flexpeak` or `both` (merged) |

The options of each peak caller for each mark (`macs2_options`, `flexpeak_options`) are in
[config/peak_settings.txt](config/peak_settings.txt).

## 4. Run

```bash
bash ~/CUTRUNpipeline/run_pipeline.sh config.sh
```

First the inputs are checked (nothing is submitted if something is wrong). Then the jobs are submitted;
each starts when the jobs it needs are finished:

| part | jobs | what |
|---|---|---|
| `qc` | 1 | FastQC and MultiQC of all raw fastq files (runs next to the alignment) |
| `align` | 1 per sample | fastp trimming, FastQC of trimmed reads, bowtie2, filtering, duplicates, bigwig |
| `peaks` | 1 per sample | MACS2 and/or FlexPeak, comparison of the two |
| `matrix` | 1 | consensus and union peaks, annotation, count matrices, QC summary |
| `diff` | 1 | DESeq2, summary figures, HTML report |

Run one part only with `run_pipeline.sh config.sh qc|align|peaks|matrix|diff`, e.g. `qc` first to check
the data, or `diff` again after changing the comparisons or `DESIGN`.

## 5. Check progress and rerun

```bash
bash ~/CUTRUNpipeline/check_status.sh config.sh
```

- **A job failed**: fix the problem and run `run_pipeline.sh config.sh` again; finished steps are skipped.
  On LSF, first `bkill` the jobs waiting for the failed one.
- **Changed peak settings**: delete `Data/Status/*.peaks.done` and run again.
- **Redo everything**: `FORCE=true` in config.sh.

## 6. Results (in WORKDIR)

| folder | files |
|---|---|
| `Report/` | **`CUTRUN_report.html`**: search a gene, peak or locus; see its peaks, differential statistics, signal profile and normalized counts; links to UCSC and IGV. One file, opens offline |
| `Figures/<mark>/` | **`Figure_<mark>_summary.pdf`**: multi-panel figure (183 mm, 7 pt): peaks per caller, FRiP, genomic features, PCA, differential peaks, volcano plots, signal profiles, signal heatmap |
| `ReadCount/<mark>/` | **count matrices of all samples** (union peaks x samples): `RawRC_<mark>.txt` (raw), `NormRC_<mark>.txt` (DESeq2-normalized), `VST_<mark>.txt` (log2 scale), `NormRC_<mark>.CondMean.txt`, `sizeFactors_<mark>.txt` (size factors, FRiP), PCA, correlation |
| `DiffPeak/<mark>/` | per comparison: `deseq2_All_...`, `deseq2_Sig_...`, `_Up.bed`, `_Down.bed`, MA, volcano and heatmap plots; `deseq2_summary_<mark>.txt` |
| `QC/` | `MultiQC/multiqc_raw_fastq.html`, `MultiQC/multiqc_report.html` (all steps), `FastQC_raw/`, `FastQC_trimmed/`, fastp reports, `Summary/sample_QC_summary.txt` + plots |
| `BAM/` | `<sample>.map.dedup.bam` (filtered, duplicates removed), alignment stats, fragment sizes |
| `Wig/` | `<sample>.rpm.bw` tracks for IGV / UCSC |
| `Peak/` | `MACS2/`, `FlexPeak/` (raw output), `peak_<sample>.macs2.bed`, `.flexpeak.bed`, `.final.bed` (used for the analysis), `.callers.txt` (MACS2 vs FlexPeak) |
| `Consensus/<mark>/` | `peak_consensus_<condition>.bed`, `peak_<mark>_union.bed`, Venn diagrams, peak width |
| `Annotation/<mark>/` | `peak_<mark>_union.anno.txt` (nearest gene, promoter/intron/...), `.anno_allgenes.txt` |
| `Log/` | job logs, settings of every run (`run_<date>/`), `software_versions.txt` |

## 7. Tests and code quality

```bash
TEST_CONDA_ENV=cutrun TEST_CONDA_ENV_R=cutrun_r bash tests/run_tests.sh
```

Makes a small simulated data set with known peaks and known differences between conditions
([tests/make_test_data.sh](tests/make_test_data.sh)), runs the whole pipeline on it (a few minutes), and
checks 90 results: input checks catch common mistakes, chrM reads are removed, lanes are merged, both peak
callers find the simulated sites, the expected differential peaks are found, the figures and the HTML
report are made, and a rerun skips finished steps.

Code follows the [Google Shell Style Guide](https://google.github.io/styleguide/shellguide.html) and the
[Google R Style Guide](https://google.github.io/styleguide/Rguide.html), checked with ShellCheck
([.shellcheckrc](.shellcheckrc)) and lintr ([.lintr](.lintr)). GitHub Actions runs the linters and the
tests on every change ([.github/workflows/tests.yml](.github/workflows/tests.yml)).

## 8. Troubleshooting

| message | what to do |
|---|---|
| `The sample list header must be ...` | save the sample list as tab-separated text with this header |
| `fastq1 not found` | check `DATADIR` and the file names |
| `conda: command not found` | load conda (`module load miniconda`), or use `MODULES` |
| `flexpeak: command not found` | install FlexPeak into the tools environment (section 1), or `PEAK_CALLERS="macs2"` |
| very few peaks | check `QC/Summary` (fragments, FRiP); loosen the cutoff of the mark in `peak_settings.txt` |
| `not enough replicates` | each condition of a comparison needs at least 2 replicates |

## 9. Methods text

> Raw reads were checked with FastQC and MultiQC, trimmed with fastp, and aligned to the genome with
> bowtie2 (`--local --very-sensitive-local --no-mixed --no-discordant -I 10 -X 700`). Properly paired
> reads with mapping quality ≥ 10 on chr1-chrY were kept, and PCR duplicates were removed with samtools
> markdup. CPM-normalized coverage tracks were made with deepTools bamCoverage. Peaks were called for each
> sample with MACS2 (BAMPE mode) and FlexPeak, with settings for each mark, and peaks in ENCODE blacklist
> regions were removed. Consensus peaks were peaks found in at least two replicates of a condition;
> consensus peaks of all conditions were merged into union peaks, and regions with IgG peaks were removed.
> Fragments in union peaks were counted with featureCounts, normalized with DESeq2 size factors, and
> differential peaks were found with DESeq2 (FDR < 0.05, fold change > 1.5). Peaks were annotated with
> ChIPseeker (promoter: −2 kb to +500 bp of the TSS) using GENCODE genes.

Exact versions are in `Log/software_versions.txt` and `DiffPeak/<mark>/sessionInfo_<mark>.txt`.

## 10. Contributing, citation, license

- **Contributing**: questions, bug reports and changes are welcome. Open an issue (bug report or feature
  request template) or a pull request. Contributors are credited in
  [CONTRIBUTORS.md](CONTRIBUTORS.md); changes are listed in [CHANGELOG.md](CHANGELOG.md).
- **Citation**: please cite this pipeline ([CITATION.cff](CITATION.cff)) and the tools it uses: FastQC,
  MultiQC, fastp, Bowtie2, SAMtools, deepTools, MACS2, FlexPeak, BEDTools, Subread/featureCounts, DESeq2,
  ChIPseeker.
- **License**: [MIT](LICENSE). You may use, change and share the pipeline, also commercially, if you keep
  the copyright and license notice.
