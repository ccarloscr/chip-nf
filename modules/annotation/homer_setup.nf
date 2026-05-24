// HOMER genome setup process
//
// Problem this solves:
//   annotatePeaks.pl requires genome-specific data files (chromosome sizes,
//   gene models, etc.) that HOMER downloads via configureHomer.pl -install.
//   This is required in ANNOTATE_PEAKS and PCA_PEAKS. If added individually
//   in each module, it runs on every iteration --> waste of resources
//
// Solution:
//   Run configureHomer.pl  once and stage the resulting data directory into every
//   process that needs it.
//
// Output:
//   homer_data/  — a self-contained directory tree that HOMER tools find when
//                  HOMER_DATA_DIR (or $CONDA_PREFIX/share/homer) is set to it.

process HOMER_SETUP {
    label 'process_high'

    conda "bioconda::homer=4.11 conda-forge::perl"

    // The populated HOMER data directory, passed as a path to downstream processes.
    output:
    path "homer_data", emit: homer_data_dir

    script:
    """
    # Copy the conda-installed HOMER data directory into the work directory
    # so Nextflow can stage it as a regular output path. We copy rather than
    # symlink because the staged path must be self-contained across nodes.
    cp -rL \$CONDA_PREFIX/share/homer homer_data

    # Download genome-specific data into the local copy.
    # -install writes files under homer_data/data/genomes/${params.genome}/
    perl homer_data/configureHomer.pl -install ${params.genome} \\
        --homer-data homer_data
    """
}
