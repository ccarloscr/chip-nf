#!/bin/bash

#SBATCH --job-name=nextflow-head
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=nextflow_%j.log

export JAVA_HOME=~/java-17
export JAVA_CMD=~/java-17/bin/java
export PATH=~/java-17/bin:$PATH

~/nextflow_versions/24.04.3/nextflow run main.nf -profile slurm -resume
