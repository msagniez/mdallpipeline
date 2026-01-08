// Import modules
include { SUBREAD_FEATURECOUNTS              } from '../../../modules/nf-core/subread/featurecounts/main.nf'     // count reads from bams on transcriptomic ref
include { ARRANGE_COUNTS                     } from '../../../modules/local/arrangecounts/main.nf'               // arrange featucounts count table to keep only geneID and raw counts
include { MDALL_COUNT_MATRIX                 } from '../../../modules/local/mdall/main.nf'                       // predict ALL subtype with MD-ALL

// Define the main workflow
workflow COUNTING {
    take:
    bam     //tuple with meta and path to bam files
    ref     //tuple with meta and path to transcriptomic annotation
    reference_dir // path to MD-ALL reference_data folder

    main:
    ch_versions = Channel.empty() // For collecting version info

    ch_count = bam
        .join(ref)
    SUBREAD_FEATURECOUNTS( ch_count )

    ARRANGE_COUNTS( SUBREAD_FEATURECOUNTS.out.counts )

    MDALL_COUNT_MATRIX( ARRANGE_COUNTS.out.arrcounts, reference_dir)

    // Collect versions from all modules
    ch_versions = SUBREAD_FEATURECOUNTS.out.versions
                  .mix(MDALL_COUNT_MATRIX.out.versions)

    emit:
    arrangedcounts = ARRANGE_COUNTS.out.arrcounts       // Final arranged count matrix
    predictions = MDALL_COUNT_MATRIX.out.classification // Final MD-ALL prediction matrix
    versions = ch_versions                              // All tool versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/