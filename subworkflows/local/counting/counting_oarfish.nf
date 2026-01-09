// Import modules
include { OARFISH_QUANTIFY_FQ                } from '../../../modules/local/oarfish/main.nf'     // count reads from fastqs on transcriptomic ref using Oarfish
//include { ARRANGE_COUNTS                     } from '../../../modules/local/arrangecounts/main.nf'               // arrange featucounts count table to keep only geneID and raw counts

// Define the main workflow
workflow COUNTING_OARFISH {
    take:
    fastq     //tuple with meta and path to fastq files
    ref     //tuple with meta and path to transcriptomic annotation

    main:
    ch_versions = Channel.empty() // For collecting version info

    ch_count = fastq
        .join(ref)
    OARFISH_QUANTIFY_FQ( ch_count )

    //ARRANGE_COUNTS( SUBREAD_FEATURECOUNTS.out.counts )

    // Collect versions from all modules
    ch_versions = OARFISH_QUANTIFY_FQ.out.versions

    emit:
    oarfish_counts = OARFISH_QUANTIFY_FQ.out.quant       // Final arranged count matrix
    versions = ch_versions                               // All tool versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/