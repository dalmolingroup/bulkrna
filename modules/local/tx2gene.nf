process TX2GENE {
    tag "$gtf"
    label "process_medium"

    conda (params.enable_conda ? "bioconda::bioconductor-tximeta=1.8.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.8.0--r40_0' :
        'quay.io/biocontainers/bioconductor-tximeta:1.8.0--r40_0' }"

    input:
    path gtf

    output:
    path "tx2gene.rds" , emit: tx2gene
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in dalmolingroup/bulkrna/bin/
    """
    gen_tx2gene.r \\
        $gtf \\
        tx2gene.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
