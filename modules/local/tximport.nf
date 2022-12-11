process TXIMPORT {
    label "process_medium"

    conda (params.enable_conda ? "bioconda::bioconductor-tximeta=1.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0' :
        'quay.io/biocontainers/bioconductor-tximeta:1.8.0--r40_0' }"

    input:
    path ("kallisto/*")
    path tx2gene

    output:
    path "gene_tpm.tsv"                 , emit: tpm_gene
    path "gene_counts.tsv"              , emit: counts_gene
    path "transcript_tpm.tsv"           , emit: tpm_transcript
    path "transcript_counts.tsv"        , emit: counts_transcript
    path "transcript_scaled_tpm.tsv"    , emit: scaled_tpm_transcript
    path "transcript_scaled_counts.tsv" , emit: scaled_counts_transcript
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in dalmolingroup/bulkrna/bin/
    """
    tximport.r \\
        kallisto \\
        $tx2gene

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
