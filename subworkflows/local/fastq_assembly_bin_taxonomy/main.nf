include { MEGAHIT                              } from '../../../modules/nf-core/megahit'
include { BWA_INDEX                            } from '../../../modules/nf-core/bwa/index'
include { BWA_MEM                              } from '../../../modules/nf-core/bwa/mem'
include { SAMTOOLS_INDEX                       } from '../../../modules/nf-core/samtools/index'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../modules/nf-core/metabat/jgisummarizebamcontigdepths/main.nf'
include { METABAT2_METABAT2                    } from '../../../modules/nf-core/metabat/metabat'
include { GTDBTK_CLASSIFYWF                    } from '../../../modules/nf-core/gtdbtk/classifywf'
include { GTDBTK_GTDBTONCBIMAJORITYVOTE        } from '../../../modules/nf-core/gtdbtk/gtdbtoncbimajorityvote/main.nf'


workflow FASTQ_ASSEMBLY_BIN_TAXONOMY {

    take:
    reads
    db
    
    main:
    ch_reads                = channel.empty()
    ch_db                   = channel.empty()
    ch_contigs              = channel.empty()
    ch_bwa_index            = channel.empty()
    ch_bins                 = channel.empty()
    ch_metabat_input        = channel.empty()
    ch_reads_contigs_index  = channel.empty()
    ch_versions             = channel.empty()
    //
    // Assembly with MEGAHIT
    //
    ch_reads = reads.map { meta, r -> [ meta, r[0], r[1] ]}
    MEGAHIT(
        ch_reads
    )
    ch_contigs  = MEGAHIT.out.contigs
    ch_versions = ch_versions.mix(MEGAHIT.out.versions)
    //
    // Index contigs with BWA
    //
    BWA_INDEX(
        ch_contigs
    )
    ch_bwa_index = BWA_INDEX.out.index
    ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
    //
    // Map reads to contigs with BWA-MEM
    //
    ch_reads_contigs_index = ch_reads.join(ch_contigs).join(ch_bwa_index)
    ch_reads_fixed   = ch_reads_contigs_index.map { meta, reads1, reads2, _contigs, _index -> [ meta, [reads1, reads2] ] }
    ch_index_fixed   = ch_reads_contigs_index.map { meta, _reads1, _reads2, _contigs, index -> [ meta, index ] }
    ch_fasta_fixed   = ch_reads_contigs_index.map { meta, _reads1, _reads2, contigs, _index -> [ meta, contigs ] }
    BWA_MEM(
        ch_reads_fixed,
        ch_index_fixed,
        ch_fasta_fixed,
        true
    )
    ch_bam      = BWA_MEM.out.bam
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    //
    // Index BAM with SAMtools
    //
    SAMTOOLS_INDEX(
        ch_bam
    )
    ch_index    = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    //
    // Generate depth file with jgi_summarize_bam_contig_depths
    //
    ch_align = ch_bam.join(ch_index)
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(
        ch_align
    )
    ch_depth      = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    ch_versions   = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.versions)
    //
    // Binning with MetaBAT2
    //
    ch_metabat_input = ch_contigs.join(ch_depth)
    METABAT2_METABAT2(
        ch_metabat_input
    )
    ch_bins      = METABAT2_METABAT2.out.fasta
    ch_versions = ch_versions.mix(METABAT2_METABAT2.out.versions)
    //
    // Taxonomy assignment with GTDB-Tk
    //
    ch_db = db.map {db_path -> [ 'gtdbtk', db_path ] }
    GTDBTK_CLASSIFYWF(
        ch_bins,
        ch_db,
        false
    )
    ch_gtdbtk_outdir    = GTDBTK_CLASSIFYWF.out.gtdb_outdir
    ch_versions         = ch_versions.mix(GTDBTK_CLASSIFYWF.out.versions)
    //
    // Converts the output classifications of GTDB-TK from GTDB taxonomy to NCBI taxonomy
    //
    ar53_metadata   = channel.value(file(params.gtdbtk_ar53_metadata)).map { file -> [ 'ar53', file ] }
    bac120_metadata = channel.value(file(params.gtdbtk_bac120_metadata)).map { file -> [ 'bac120', file ] }
    ch_input        = ch_gtdbtk_outdir.map { meta, outdir -> [ meta, outdir, [] ] }
    GTDBTK_GTDBTONCBIMAJORITYVOTE(
        ch_input,
        ar53_metadata,
        bac120_metadata
    )
    ch_versions = ch_versions.mix(GTDBTK_GTDBTONCBIMAJORITYVOTE.out.versions)

    emit:
    versions    = ch_versions
    contigs     = MEGAHIT.out.contigs
}