include { MEGAHIT } from '../../../modules/nf-core/megahit'

workflow FASTQ_ASSEMBLY_BIN_TAXONOMY {

    take:
    reads
    db
    
    main:
    ch_reads = Channel.empty()
    ch_reads = reads.map { meta, reads -> [ meta, reads[0], reads[1] ]}

    MEGAHIT(
        ch_reads
    )
    emit:
    versions    = MEGAHIT.out.versions
    contigs     = MEGAHIT.out.contigs
}