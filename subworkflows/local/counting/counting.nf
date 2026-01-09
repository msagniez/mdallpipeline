// Import modules
include { SUBREAD_FEATURECOUNTS              } from '../../../modules/nf-core/subread/featurecounts/main.nf'     // count reads from bams on transcriptomic ref using featurecounts
include { ARRANGE_COUNTS                     } from '../../../modules/local/arrangecounts/main.nf'               // arrange featucounts count table to keep only geneID and raw counts

// Define the main workflow
workflow COUNTING {
    take:
    bam     //tuple with meta and path to bam files
    ref     //tuple with meta and path to transcriptomic annotation

    main:
    ch_versions = Channel.empty() // For collecting version info

    ch_count = bam
        .join(ref)
    SUBREAD_FEATURECOUNTS( ch_count )

    ARRANGE_COUNTS( SUBREAD_FEATURECOUNTS.out.counts )

    // Collect versions from all modules
    ch_versions = SUBREAD_FEATURECOUNTS.out.versions

    emit:
    arrangedcounts = ARRANGE_COUNTS.out.arrcounts       // Final arranged count matrix
    versions = ch_versions                            // All tool versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/