#!/bin/bash
#
# Settings of one CUT&RUN project.
# Copy this file, sample_list.txt and comparison_list.txt into your project
# folder, edit them, then run:
#   bash ~/CUTRUNpipeline/run_pipeline.sh config.sh
# Use full paths. No spaces around "=".

# 1. Project
PROJECT_NAME="MyCutRun"                        # short name for cluster jobs
WORKDIR="/home/you/MyCutRun/Data"              # all results go here
DATADIR="/home/you/MyCutRun/Raw"               # raw fastq files
SAMPLE_LIST="/home/you/MyCutRun/sample_list.txt"
COMPARISON_LIST="/home/you/MyCutRun/comparison_list.txt"  # empty: no DESeq2

# MACS2 settings of each mark. To change them, copy the file and point here.
PEAK_SETTINGS="${PIPELINE_DIR}/config/peak_settings.txt"

# 2. Reference genome (made by reference/build_reference.sh)
SPECIES="mm10"
LIBRARIES="/home/you/indexes"                  # one folder per genome
GNM_IDX="${LIBRARIES}/${SPECIES}"
BOWTIE_IDX="${GNM_IDX}/bowtie2/${SPECIES}"     # without .1.bt2
GNM_SIZES="${GNM_IDX}/${SPECIES}.genomesize"
BLACKLIST="${GNM_IDX}/${SPECIES}-blacklist.clean.bed"
CHR1TOY_BED="${GNM_IDX}/${SPECIES}.chr1toY.bed"  # keep reads on chr1-chrY
TXDB="${GNM_IDX}/genes.txdb.sqlite"            # genes for ChIPseeker
GENE_TABLE="${GNM_IDX}/genes.txt"              # gene ID -> gene name
REFSIZE="1.87e9"                               # mouse 1.87e9, human 2.7e9
HOMER_GENOME=""                                # e.g. mm10 (empty: no Homer)

# 3. Software: conda environments from envs/*.yml (empty: tools on PATH)
CONDA_ENV="cutrun"
CONDA_ENV_R="cutrun_r"
MODULES=""                                     # or modules, e.g. "bowtie2"
MODULES_R=""                                   # e.g. "R/4.3"

# 4. Cluster
SCHEDULER="lsf"                                # lsf | slurm | local
QUEUE=""                                       # empty: default queue
THREADS=8                                      # CPUs per job
MEMORY=32000                                   # memory per job (MB)
EXTRA_SUBMIT_ARGS=""                           # e.g. "-W 24:00"

# 5. Trimming
RUN_FASTQC=true
RUN_TRIMMING=true                              # fastp adapter trimming
FASTP_ARGS="--detect_adapter_for_pe --length_required 20"

# 6. Alignment (Henikoff lab settings: local, proper pairs, 10-700 bp)
BOWTIE2_ARGS="--local --very-sensitive-local --no-unal --no-mixed \
--no-discordant --phred33 -I 10 -X 700"
MIN_MAPQ=10                                    # 0: keep multi-mapping reads
REMOVE_DUPLICATES=true                         # false: only mark them
SPIKEIN_IDX=""                                 # optional spike-in index
SPIKEIN_SCALE=10000                            # scale / spike-in fragments

# 7. Bigwig
BW_NORM="CPM"                                  # CPM | spikein
BW_ARGS="--binSize 10"                         # more bamCoverage options

# 8. Peaks
PEAK_CALLERS="macs2 flexpeak"                  # macs2 and/or flexpeak
ANALYSIS_PEAKS="macs2"                         # peaks for consensus/union:
                                               # macs2 | flexpeak | both
MACS2_FORMAT="BAMPE"                           # real fragment of each pair
IGG_MARK="IgG"                                 # "mark" of IgG samples
IGG_MODE="subtract"                            # subtract | control | none
IGG_EXTEND=200                                 # extend IgG regions (bp)
MIN_REPS=2                                     # consensus: in >= 2 reps

# 9. Annotation
TSS_UP=2000                                    # promoter: TSS - 2 kb
TSS_DOWN=500                                   #       to: TSS + 500 bp
FLANK=10000                                    # genes within 10 kb

# 10. Differential peaks (DESeq2)
DESIGN="~ condition"                           # or "~ batch + condition"
MIN_COUNT=15                                   # max raw count of a peak
PADJ=0.05                                      # FDR cutoff
LFC=0.585                                      # |log2FC| (0.585: 1.5 fold)

# true: redo all steps, also finished ones
FORCE=false
