# chip-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen)](https://www.nextflow.io/)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-276DC3?logo=r)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue?logo=python)](https://www.python.org/)
[![Conda](https://img.shields.io/badge/conda-enabled-44903d?logo=anaconda)](https://docs.conda.io)

**Nextflow-based pipeline for the end-to-end analysis of ChIP-seq datasets**

A modular, reproducible Nextflow (DSL2) pipeline for end-to-end ChIP-seq data processing: from raw paired-end FASTQ files to differential peak analysis, annotation, and visualizations.

---

## Repository structure

```
chip-nf/
├── main.nf                            # Main pipeline entry point
├── nextflow.config                    # Pipeline configuration and parameters
├── metadata.csv                       # Sample sheet (sample, fastq_1, fastq_2, control)
├── run_chipseq.sh                     # Convenience script to run the pipeline using SLURM
├── README.md
├── LICENSE
├── .gitignore
│
├── modules/                           # Modular Nextflow process definitions
│   ├── qc/
│   │   ├── fastqc.nf                  # FastQC on raw and trimmed reads
│   │   ├── fastp.nf                   # Adapter trimming and quality filtering
│   │   └── multiqc.nf                 # Aggregate QC report
│   ├── alignment/
│   │   └── bowtie2.nf                 # Bowtie2 alignment → sorted, deduplicated BAM
│   ├── peaks/
│   │   ├── macs2_narrow.nf            # MACS2 narrow peak calling (TFs)
│   │   ├── macs2_broad.nf             # MACS2 broad peak calling (histones)
│   │   ├── merge_peaks.nf             # Consensus peak set across replicates
│   │   └── deseq2_peaks.nf            # Differential peak analysis with DESeq2
│   ├── annotation/
│   │   ├── homer_setup.nf             # Setup the HOMER installation files
│   │   └── homer_annot.nf             # Peak annotation with HOMER
│   └── visualization/
│       ├── pca_peaks.nf               # PCA on consensus peak count matrix
│       ├── volcano_plot.nf            # Volcano plot of differential peaks
│       ├── plot_heatmap.nf            # Signal heatmap around peaks
│       └── genomic_distribution.nf      # Stacked bar: promoter vs intergenic distribution
│
├── bin/
│   └── pca_merged_peaks.R             # R helper script for PCA computation
│
├── Genomes/                           # Set up this directory with your own genome fasta and GTF
│   └── dm6
│ 
├── data/                              # Raw FASTQ input files (not tracked by git)
│   └── *.fastq.gz
│
└── expected_outputs/                  # Reference outputs for pipeline validation
```

---

## Features

- Quality control of raw and trimmed reads (FastQC + fastp)
- Paired-end alignment with Bowtie2, filtering by MAPQ
- MACS2 peak calling in **narrow** (TF) or **broad** (histone) mode
- Consensus peak set generation with BEDTools
- Peak quantification and PCA across all replicates (HOMER + R/ggplot2)
- Differential peak analysis between conditions (DESeq2)
- Volcano plot of significant differential peaks
- Genomic distribution stacked bar chart (HOMER annotation + R)
- deepTools heatmap of ChIP-seq signal at consensus peaks
- Aggregated MultiQC report (FastQC, fastp, Bowtie2, samtools, MACS2)
- Per-process conda environments
- `local` and `slurm` execution profiles

---

## Pipeline overview

```
FASTQ (paired-end)
    │
    ├─► FastQC (raw)
    │
    ├─► fastp (adapter trimming & quality filtering)
    │       └─► FastQC (trimmed)
    │
    └─► Bowtie2 alignment → samtools sort/index/flagstat
            │
            ├─► [ChIP BAMs] ──┐
            └─► [Input BAMs] ─┤
                              │
                        MACS2 peak calling
                        (narrow or broad)
                              │
                    ┌─────────┴──────────────┐
                    │                        │
              HOMER annotation         BEDTools merge
              (per sample)            (consensus peaks)
                    │                        │
              Genomic distribution    HOMER quantification
              stacked bar plot               │
                                    ┌────────┴──────┐
                                    │               │
                                  PCA             DESeq2
                                  plot              │
                                              Volcano plot
                                            deepTools heatmap
                                                    │
                                              MultiQC report
```

---

## Requirements

| Tool         | Version  | Notes                          |
|--------------|----------|--------------------------------|
| Nextflow     | ≥ 23.04  | `NXF_VER=23.10.0` recommended |
| conda/mamba  | any      | mamba strongly recommended     |
| Bowtie2      | 2.5.1    | installed via conda            |
| samtools     | 1.18     | installed via conda            |
| fastp        | 0.23.4   | installed via conda            |
| FastQC       | 0.12.1   | installed via conda            |
| MACS2        | 2.2.9.1  | installed via conda            |
| HOMER        | 4.11     | installed via conda            |
| BEDTools     | 2.31.1   | installed via conda            |
| DESeq2       | ≥ 1.40   | installed via conda (R)        |
| deepTools    | 3.5.4    | installed via conda            |
| MultiQC      | 1.21     | installed via conda            |
| R            | ≥ 4.3    | with ggplot2, ggrepel, dplyr   |

All software is managed automatically through per-process conda environments. No manual installation is required beyond Nextflow and conda/mamba.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/ccarloscr/chip-nf.git
cd chip-nf

# 2. (Optional but recommended) install mamba for faster env solving
conda install -n base -c conda-forge mamba

# 3. Download the reference genome in fasta format and build a Bowtie2 index
#    Example for dm6:
mkdir -p Genomes/dm6
bowtie2-build dm6.fa Genomes/dm6/dm6

# 4. Download the genome GTF annotation
cp dm6_annotation.gtf Genomes/dm6/
```

---

## Usage

### Local execution

```bash
nextflow run main.nf -profile standard
```

### SLURM cluster

```bash
nextflow run main.nf -profile slurm
or
sbatch run_chipnf.sh
```

### Resume a run after failure

```bash
nextflow run main.nf -profile slurm -resume
```

---

## Samplesheet format

The pipeline expects a comma-separated file (`metadata.csv`) with the following columns:

| Column    | Required | Description                                               |
|-----------|----------|-----------------------------------------------------------|
| `sample`  | yes      | Unique sample identifier (e.g. `WT_rep1`)                 |
| `fastq_1` | yes      | Path to R1 FASTQ file (gzipped)                          |
| `fastq_2` | yes      | Path to R2 FASTQ file (gzipped)                          |
| `control` | no       | Sample ID of the matched Input/control; leave empty for controls |

**Naming convention:** condition and replicate are parsed from the sample ID using `_` as separator. For example, `WT_rep1` → condition `WT`, replicate `rep1`.

```csv
sample,fastq_1,fastq_2,control
WT_rep1,data/WT_rep1_1.fastq.gz,data/WT_rep1_2.fastq.gz,Input_WT_rep1
WT_rep2,data/WT_rep2_1.fastq.gz,data/WT_rep2_2.fastq.gz,Input_WT_rep2
Input_WT_rep1,data/Input_WT_rep1_1.fastq.gz,data/Input_WT_rep1_2.fastq.gz,
Input_WT_rep2,data/Input_WT_rep2_1.fastq.gz,data/Input_WT_rep2_2.fastq.gz,
PH_rep1,data/PH_rep1_1.fastq.gz,data/PH_rep1_2.fastq.gz,Input_PH_rep1
PH_rep2,data/PH_rep2_1.fastq.gz,data/PH_rep2_2.fastq.gz,Input_PH_rep2
Input_PH_rep1,data/Input_PH_rep1_1.fastq.gz,data/Input_PH_rep1_2.fastq.gz,
Input_PH_rep2,data/Input_PH_rep2_1.fastq.gz,data/Input_PH_rep2_2.fastq.gz,
```

---

## Parameters

All parameters are defined in `nextflow.config` and can be overridden at the command line with `--param_name value`.

### Core

| Parameter    | Default                                 | Description                         |
|--------------|-----------------------------------------|-------------------------------------|
| `genome`     | `dm6`                                   | Genome assembly name                |
| `metadata`   | `${projectDir}/metadata.csv`            | Path to samplesheet                 |
| `output_dir` | `${projectDir}/Results`                 | Top-level output directory          |

### Trimming (fastp)

| Parameter                | Default | Description                                 |
|--------------------------|---------|---------------------------------------------|
| `fastp_cut_mean_quality` | `20`    | Sliding-window mean quality threshold       |
| `fastp_min_read_length`  | `30`    | Minimum read length after trimming (bp)     |

### Alignment (Bowtie2)

| Parameter         | Default                                      | Description                              |
|-------------------|----------------------------------------------|------------------------------------------|
| `bowtie2_index`   | `${projectDir}/Genomes/${genome}/${genome}/…` | Bowtie2 index prefix                     |
| `min_mapq`        | `30`                                         | Minimum MAPQ to retain aligned reads     |
| `max_insert_size` | `2000`                                       | Maximum fragment insert size (`-X` flag) |

### Peak calling (MACS2)

| Parameter      | Default   | Description                                         |
|----------------|-----------|-----------------------------------------------------|
| `peak_mode`    | `broad`   | `broad` for histones, `narrow` for transcription factors |
| `genome_size`  | `dm`      | Effective genome size passed to MACS2 (`-g`)        |
| `min_qval`     | `0.01`    | Minimum q-value threshold for peak calling          |
| `broad_cutoff` | `0.1`     | Broad peak score cutoff (`--broad-cutoff`)          |

### Annotation

| Parameter  | Default                                       | Description          |
|------------|-----------------------------------------------|----------------------|
| `gtf_file` | `${projectDir}/Genomes/${genome}/…_annotation.gtf` | Genome GTF file |

### Visualization

| Parameter      | Default | Description                              |
|----------------|---------|------------------------------------------|
| `volcano_top`  | `20`    | Number of top peaks to label on volcano plot |

---

## Output structure

```
Results/
├── QC/
│   ├── FastQC/               # Raw and trimmed FastQC reports (.zip, .html)
│   ├── Fastp/                # fastp HTML + JSON reports
│   └── MultiQC/              # Aggregated MultiQC report
├── trimmed/                  # Trimmed FASTQ files
├── Aligned/                  # Sorted BAM files + indices + Bowtie2 logs + flagstat
├── Peaks/                    # MACS2 peak files (.narrowPeak or .broadPeak, .xls)
├── Annotated/                # HOMER-annotated peak files + annotation stats
├── DESeq2/                   # DESeq2 differential results table (TSV)
└── Visualization/
    ├── PCA/                  # PCA plot (PNG) + HOMER count matrix
    ├── Volcano/              # Volcano plot (PNG)
    ├── GenomicDistribution/  # Stacked bar plot (PNG) + summary table (TSV)
    └── Heatmap/              # deepTools heatmap (PNG + PDF) + score matrix
```

---

## Notes

**Samplesheet conventions**
- The control column must be left empty for input samples
- Sample IDs are expected to follow the CONDITION_repN naming pattern (e.g. WT_rep1)

**Genome setup**
- The Bowtie2 index must be pre-built and placed at Genomes/<genome>/<genome>/<genome> relative to the project directory, matching the path hardcoded in `nextflow.config`
- The GTF annotation file is also expected under Genomes/<genome>/ — both need to be in place before running
- `genome_size` in the config uses MACS2 shorthand (dm, hs, mm…), not the full genome name

**Peak calling mode**
- `peak_mode = "broad"` is the default and intended for histone marks; switch to `"narrow" for transcription factors — this affects both MACS2 calling and which .xls files get passed to MultiQC

**Conda/Mamba**
- The pipeline uses per-module conda environments (not a single environment.yml), so Mamba must be installed and on PATH, as `conda.useMamba = true` is set globally
- First run will be slow while environments are built; subsequent runs with `-resume` reuse the cache at `conda_envs/`

**Resource defaults**
- Default is 4 CPUs / 8 GB per process; `process_high` label gets 12 CPUs / 24 GB
- Change `process.queue` in `nextflow.config` to match the name of your cluster queue
- `executor.queueSize = 4` limits concurrent jobs, raise it if your cluster allows more parallelism



---

## License

This project is licensed under the MIT License.
