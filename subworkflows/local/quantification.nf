include { KALLISTO_INDEX } from '../../modules/nf-core/kallisto/index/main.nf'
include { KALLISTO_QUANT } from '../../modules/local/kallisto/quant/main.nf'
include { TX2GENE } from '../../modules/local/tx2gene.nf'
include { TXIMPORT } from '../../modules/local/tximport.nf'

workflow QUANTIFICATION {
    take:
        reads            // channel: [ val(meta), [ reads ] ]
        transcriptome
        gtf
        fragment_length
        fragment_length_sd

    main:
        ch_versions = Channel.empty()

        KALLISTO_INDEX ( transcriptome )
        KALLISTO_QUANT (
            reads,
            KALLISTO_INDEX.out.idx,
            fragment_length,
            fragment_length_sd
        )
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())

        TX2GENE (
            gtf
        )

        TXIMPORT (
            KALLISTO_QUANT.out.results.collect{it[1]},
            TX2GENE.out.tx2gene
        )

        ch_versions = ch_versions.mix(TXIMPORT.out.versions)

    emit:
        results = KALLISTO_QUANT.out.results
        logs = KALLISTO_QUANT.out.log
        versions = ch_versions
}
