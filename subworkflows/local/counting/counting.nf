// Import modules
include { SUBREAD_FEATURECOUNTS              } from '../../../modules/nf-core/subread/featurecounts/main.nf'     // count reads from bams on transcriptomic ref

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

    // Collect versions from all modules
    ch_versions = SUBREAD_FEATURECOUNTS.out.versions

    emit:
    bam      = SUBREAD_FEATURECOUNTS.out.counts       // Final sorted BAM with index
    versions = ch_versions                            // All tool versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/