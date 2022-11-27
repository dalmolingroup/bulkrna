process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1' :
        'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1' }"

    input:
    tuple val(meta), path(reads)
    path reference
    val fragment_length
    val fragment_length_sd

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val(meta), path("*.log")      , emit: log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (meta.single_end) {
    """
    kallisto quant \
        -t $task.cpus \
        -i $reference \
        -o $meta.id \
        --single \
        --fragment-length $fragment_length \
        --sd $fragment_length_sd \
        $reads \
        2> ${meta.id}.kallisto.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """
    } else {
    """
    kallisto quant \
        -t $task.cpus \
        -i $reference \
        -o $meta.id \
        $reads \
        2> ${meta.id}.kallisto.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto 2>&1) | sed 's/^kallisto //; s/Usage.*\$//')
    END_VERSIONS
    """
    }
}