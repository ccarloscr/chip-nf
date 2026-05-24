// MultiQC process: aggregates QC reports from FastQC, fastp, Bowtie2,
// samtools flagstat, samtools markdup, and MACS2.
// Each tool's outputs are staged into named subdirectories so MultiQC's
// module auto-detection works correctly even when file names overlap.
process MULTIQC {
    label 'process_low'

    publishDir "${params.output_dir}/QC/MultiQC", mode: 'copy'

    conda "bioconda::multiqc=1.21"

    input:
    path fastqc_raw_reports,	   stageAs: "fastqc_raw/*"	// Raw FastQC .zip and .html
    path fastqc_trimmed_reports,   stageAs: "fastqc_trimmed/*"  // Trimmed FastQC .zip and .html
    path fastp_json,                stageAs: "fastp/*"            // Fastp json report"
    path bowtie2_logs,             stageAs: "bowtie2/*"         // Bowtie2 .log (stderr captured during alignment)
    path flagstat_logs,            stageAs: "samtools/*"        // samtools flagstat output files
    path markdup_logs,              stageAs: "samtools/*"        // samtools markdup output files
    path macs2_logs,               stageAs: "macs2/*"           // MACS2 _peaks.xls files

    output:
    path "multiqc_report.html",      emit: report
    path "multiqc_report_data/",     emit: data

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --dirs \\
        --dirs-depth 1
    """
}
