/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAPPING                } from '../subworkflows/local/mapping/mapping.nf'
include { COUNTING               } from '../subworkflows/local/counting/counting.nf'
include { COUNTING_OARFISH       } from '../subworkflows/local/counting/counting_oarfish.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mdallpipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MDALLPIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_mapping = Channel.empty()

    //
    // MODULE: Run MAPPING
    //prep channel
    ch_mapping = ch_samplesheet.map { meta, libtype, fastqFiles, gref, tref -> 
    tuple(meta, libtype, fastqFiles, gref)
    }
    //run MAPPING
    MAPPING ( ch_mapping )

    // MODULE: Run COUNTING with FeatureCounts
    //prep channel
    ch_bam = MAPPING.out.bam.map { meta, bam, bai ->
    tuple(meta,bam)
    }
    ch_ref = ch_samplesheet.map { meta, libtype, fastqFiles, gref, tref -> 
    tuple(meta, tref)
    }
    //run COUNTING with FeatureCounts
    COUNTING( ch_bam, ch_ref )

    // MODULE: Run COUNTING_OARFISH with Oarfish
    //prep channel
    ch_fastq = MAPPING.out.pychopped_fastq
    ch_ref = ch_samplesheet.map { meta, libtype, fastqFiles, gref, tref -> 
    tuple(meta, tref)
    }
    //run COUNTING_OARFISH with Oarfish
    COUNTING_OARFISH( ch_bam, ch_ref )
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'mdallpipeline_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    

    emit: 
    multiqc_report = Channel.empty() 
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
