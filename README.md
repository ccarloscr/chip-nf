# chip_nf

This workflow processes, maps and analises pre-filtered .fastq files from ChIP-seq experiments. The pipeline is optimized for histone marks, but the fragment length parameters can be adjusted for other proteins.

The main Nextflow script [`chipseq.nf`](chipseq.nf) orchestrates the pipeline by calling the scripts located in the [`Scripts/`](Scripts/) directory. Each script corresponds to a specific step in the workflow:

- [`Mapping.sh`](Scripts/Mapping.sh): Read mapping using HISAT2.
- [`Post-map-process.sh`](Scripts/Post-map-process.sh): Filtering, sorting and indexing of aligned files.
- [`Peak-calling.sh`](Scripts/Peak-calling.sh): Peak calling using MACS2.
- [`Peak-annotation.R`](Scripts/Peak-annotation.R): LiftOver from dm3 to dm6 (only if using dm3); Filtering of non-canonical chromosomes; Annotation of peak features.


## Installation

To install the pipeline clone the repository:
```bash
git clone https://github.com/ccarloscr/Chipseq-analysis.git
cd Chipseq-analysis
```

This workflow depends on multiple tools and libraries, which are installed in the Conda environment [environment.yml](environment.yml). Once the environment is created, it will be called automatically by the [`chipseq.nf`](chipseq.nf) script.

To create the required conda environment:
```bash
conda env create -f environment.yml -n chipseq_env
```

## Build the reference genome

The [`Mapping.sh`](Scripts/Mapping.sh) script uses HISAT2 for the alignment of reads. HISAT2 requires the reference genome to work. For this, you should download the fasta file of your reference genome and build it using the hisat2-build command.

The code below is used to download and build the dm3 or dm6 genomes of _Drosophila melanogaster_. For other genomes, search for the fasta file of your reference genome here https://hgdownload.soe.ucsc.edu/downloads.html.

#### dm3 genome
```bash
# Create the directory:
mkdir -p ~/Chipseq-analysis/Genomes/dm3
cd ~/Chipseq-analysis/Genomes/dm3

# Download the dm3 genome from UCSC
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz
gunzip dm3.fa.gz

# Activate the conda environment to get access to HISAT2
conda activate chipseq_env

# Build the dm3 genome
hisat2-build dm3.fa dm3_index
```

#### dm6 genome
```bash
# Create the directory:
mkdir -p ~/Chipseq-analysis/Genomes/dm6
cd ~/Chipseq-analysis/Genomes/dm6

# Download the dm6 genome from UCSC
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6.fa.gz

# Activate the conda environment to get access to HISAT2
conda activate chipseq_env

# Build the dm6 genome
hisat2-build dm6.fa dm6_index
```


## Configuration

#### Change the default parameters in chipseq.nf
The [`chipseq.nf`](chipseq.nf) script uses several parameters defined in the first lines (Default parameters). You should only change the following:

- params.home_dir: change to the absolute path were you cloned this repository.
- params.genome: reference genome used for mapping (default: dm3).
- params.max_mismatch: maximum number of mismatches allowed during the mapping step (default: 4).
- params.ext_size: average chip-seq fragment length or maximum peak size (default: 150).

#### Change the nextflow.config file
The main nextflow script [`chipseq.nf`](chipseq.nf) uses the [`nextflow.config`](nextflow.config) configuration. The provided options use SLURM to connect to the irbio01 HPC cluster. Change the [`nextflow.config`](nextflow.config) according to your preferences.

#### Change or ignore the run_chipseq.nf script
The [`run_chipseq.nf`](run_chipseq.nf) script is sent using sbatch to the cluster. Ignore if you prefer running the [`chipseq.nf`](chipseq.nf) script directly or configure the SBATCH options according to your preferences.

#### Change the metadata.csv file
The [`metadata.csv`](metadata.csv) file contains the necessary information to pair the input and experimental samples. In the first and second columns, you should place the sample names of the input and the corresponding experimental sample, respectively. Columns 3 to 4 are used to rename the output files according to the antibody and genotype of the samples. Column 5 designs the replicate number.


## License
This project is licensed under the MIT License.
