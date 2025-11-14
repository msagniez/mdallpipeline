// Import modules
include { CAT_FASTQ              } from '../../../modules/local/catfastq/main.nf'
include { PYCHOPPER              } from '../../../modules/nf-core/pychopper/main.nf'      // orient and trim cDNA reads with pychopper
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_CDNA  } from '../../../modules/local/minimap2/main.nf'         // minimap2 alignment for pychopped cDNA reads
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_DRNA  } from '../../../modules/local/minimap2/main.nf'         // minimap2 alignment for dRNA reads
include { SAMTOOLS_TOBAM         } from '../../../modules/local/samtools/main.nf'         // Convert SAM to BAM
include { SAMTOOLS_SORT          } from '../../../modules/local/samtools/main.nf'         // Sort BAM
include { SAMTOOLS_INDEX         } from '../../../modules/local/samtools/main.nf'         // Index BAM

// Define the main workflow
workflow MAPPING {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty() // For collecting version info

    // Branch the samplesheet based on library type
    ch_samplesheet
        .branch { meta, libtype, fastqFiles, gref ->
            cdna: libtype == 'cDNA'
                return tuple(meta, libtype, fastqFiles, gref)
            drna: libtype == 'dRNA'
                return tuple(meta, libtype, fastqFiles, gref)
        }
    .set { branched_samples }

    // Concatenate FASTQ files for cDNA samples
    ch_cdna_for_concat = branched_samples.cdna
        .map { meta, libtype, fastqFiles, gref ->
            tuple(meta, fastqFiles)
        }

    CAT_FASTQ(ch_cdna_for_concat)

    // Run Pychopper on concatenated FASTQ
    PYCHOPPER(CAT_FASTQ.out.fastq)

    ch_mm2_cdna = PYCHOPPER.out.fastq
       .join(branched_samples.cdna.map { meta, libtype, fastqFiles, gref ->
            tuple(meta, gref)
        })

    MINIMAP2_ALIGN_CDNA(ch_mm2_cdna)

    ch_sam_cdna = MINIMAP2_ALIGN_CDNA.out.sam

    // Process dRNA samples directly
    ch_mm2_drna = branched_samples.drna
        .map { meta, libtype, fastqFiles, gref ->
            tuple(meta, fastqFiles, gref)
        }
    MINIMAP2_ALIGN_DRNA(ch_mm2_drna)

    ch_sam_drna = MINIMAP2_ALIGN_DRNA.out.sam

    // Merge all SAM outputs
    ch_all_sam = ch_sam_cdna
        .mix(ch_sam_drna)

    // Convert SAM to BAM
    SAMTOOLS_TOBAM(ch_all_sam)

    // Sort and index BAM
    SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sortedbam)

    // Collect versions from all modules
    ch_versions = PYCHOPPER.out.versions
        .mix(MINIMAP2_ALIGN_CDNA.out.versions)
        .mix(MINIMAP2_ALIGN_DRNA.out.versions)
        .mix(SAMTOOLS_TOBAM.out.versions)
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam      = SAMTOOLS_INDEX.out.bamfile_index       // Final sorted BAM with index
    versions = ch_versions                            // All tool versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/