include { KALLISTO_INDEX } from '../../modules/nf-core/kallisto/index/main.nf'

workflow QUANTIFICATION {
    take:
        transcriptome

    main:
        KALLISTO_INDEX ( transcriptome )

    emit:
        versions = KALLISTO_INDEX.out.versions
}