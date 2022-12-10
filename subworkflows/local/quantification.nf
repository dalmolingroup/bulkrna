include { KALLISTO_INDEX } from '../../modules/nf-core/kallisto/index/main.nf'
include { KALLISTO_QUANT } from '../../modules/local/kallisto/quant/main.nf'
include { TX2GENE } from '../../modules/local/tx2gene.nf'

workflow QUANTIFICATION {
    take:
        reads            // channel: [ val(meta), [ reads ] ]
        transcriptome
        gtf
        fragment_length
        fragment_length_sd

    main:
        KALLISTO_INDEX ( transcriptome )
        KALLISTO_QUANT (
            reads,
            KALLISTO_INDEX.out.idx,
            fragment_length,
            fragment_length_sd
        )
        TX2GENE (
            gtf
        )

    emit:
        results = KALLISTO_QUANT.out.results
        logs = KALLISTO_QUANT.out.log
        tx2gene = TX2GENE.out.tx2gene
        kallisto_versions = KALLISTO_QUANT.out.versions
        tximport_versions = TX2GENE.out.versions
}
