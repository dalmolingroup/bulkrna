include { KALLISTO_INDEX } from '../../modules/nf-core/kallisto/index/main.nf'
include { KALLISTO_QUANT } from '../../modules/local/kallisto/quant/main.nf'

workflow QUANTIFICATION {
    take:
        reads            // channel: [ val(meta), [ reads ] ]
        transcriptome
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

    emit:
        results = KALLISTO_QUANT.out.results
        logs = KALLISTO_QUANT.out.log
        versions = KALLISTO_QUANT.out.versions
}