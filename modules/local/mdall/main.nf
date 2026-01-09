#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MDALL_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'
    
    // Use module system instead of conda
    module 'R/4.2.0:DESeq2/1.38.0:SummarizedExperiment/1.28.0'
    
    input:
    tuple val(meta), path(count_file), path(vcf_file), path(fusioncatcher_file), path(cicero_file)
    path(reference_dir)  // Directory containing R objects and functions

    output:
    tuple val(meta), path("${meta.id}_classification.tsv"), emit: classification
    tuple val(meta), path("${meta.id}_details.rds"), emit: details
    path "versions.yml", emit: versions

    script:
    def fusioncatcher = fusioncatcher_file.name != 'NO_FILE' ? "\"${fusioncatcher_file}\"" : '""'
    def cicero = cicero_file.name != 'NO_FILE' ? "\"${cicero_file}\"" : '""'
    def min_read_cnt = task.ext.minReadCnt ?: 3
    def min_depth = task.ext.minDepth ?: 20
    def maf_min = task.ext.mafmin ?: 0.1
    def maf_max = task.ext.mafmax ?: 0.85
    def feature_n = task.ext.featureN_PG ?: "100,200,300,400,500,600,700,800,900,1000,1058"
    
    """
    #!/usr/bin/env Rscript
    
    # Set working directory
    setwd("${reference_dir}")
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(DESeq2)
        library(SummarizedExperiment)
        library(dplyr)
        library(umap)
        library(Rphenograph)
        library(e1071)
        library(stringr)
        library(vroom)
    })
    
    # Source MD-ALL functions
    source("${reference_dir}/mdall_core_functions.R")
    
    # Load reference objects
    cat("Loading reference data...\\n")
    load("${reference_dir}/obj_234_HTSeq.RData")
    load("${reference_dir}/obj_1821.RData")
    load("${reference_dir}/models_svm.RData")
    
    # Parse feature numbers
    featureN_PG <- as.numeric(unlist(strsplit("${feature_n}", ",")))
    
    # Run classification
    cat("\\nStarting MD-ALL classification for ${meta.id}\\n")
    cat("================================================\\n")
    
    result <- tryCatch({
        run_one_sample(
            sample_id = "${meta.id}",
            file_count = "${count_file}",
            file_vcf = "${vcf_file}",
            file_fusioncatcher = ${fusioncatcher},
            file_cicero = ${cicero},
            featureN_PG = featureN_PG,
            minReadCnt = ${min_read_cnt},
            minDepth = ${min_depth},
            mafmin = ${maf_min},
            mafmax = ${maf_max}
        )
    }, error = function(e) {
        cat("ERROR in classification:\\n")
        cat(conditionMessage(e), "\\n")
        data.frame(
            sample_id = "${meta.id}",
            error = conditionMessage(e),
            stringsAsFactors = FALSE
        )
    })
    
    # Write output
    write.table(result, 
                "${meta.id}_classification.tsv",
                sep = "\\t",
                row.names = FALSE,
                quote = FALSE)
    
    # Save detailed results if successful
    if (!"error" %in% colnames(result)) {
        saveRDS(result, "${meta.id}_details.rds")
    }
    
    # Write versions
    version_info <- c(
        '"${task.process}":',
        paste0('    r-base: "', R.version.string, '"'),
        paste0('    deseq2: "', packageVersion("DESeq2"), '"'),
        paste0('    rphenograph: "', packageVersion("Rphenograph"), '"'),
        paste0('    summarizedexperiment: "', packageVersion("SummarizedExperiment"), '"')
    )
    writeLines(version_info, "versions.yml")
    
    cat("\\n================================================\\n")
    cat("Classification complete for ${meta.id}\\n")
    """
}

process MDALL_BATCH_CLASSIFY {
    tag "batch_${meta.id}"
    label 'process_high'
    
    module 'R/4.2.0:DESeq2/1.38.0:SummarizedExperiment/1.28.0'

    input:
    tuple val(meta), path(sample_listing)
    path(reference_dir)

    output:
    tuple val(meta), path("${meta.id}_batch_classification.tsv"), emit: classification
    tuple val(meta), path("${meta.id}_batch_details.rds"), emit: details
    path "versions.yml", emit: versions

    script:
    def min_read_cnt = task.ext.minReadCnt ?: 3
    def min_depth = task.ext.minDepth ?: 20
    def maf_min = task.ext.mafmin ?: 0.1
    def maf_max = task.ext.mafmax ?: 0.85
    def feature_n = task.ext.featureN_PG ?: "100,200,300,400,500,600,700,800,900,1000,1058"
    
    """
    #!/usr/bin/env Rscript
    
    setwd("${reference_dir}")
    
    suppressPackageStartupMessages({
        library(DESeq2)
        library(SummarizedExperiment)
        library(dplyr)
        library(umap)
        library(Rphenograph)
        library(e1071)
        library(stringr)
        library(vroom)
    })
    
    source("${reference_dir}/mdall_core_functions.R")
    
    cat("Loading reference data...\\n")
    load("${reference_dir}/obj_234_HTSeq.RData")
    load("${reference_dir}/obj_1821.RData")
    load("${reference_dir}/models_svm.RData")
    
    featureN_PG <- as.numeric(unlist(strsplit("${feature_n}", ",")))
    
    cat("\\nProcessing batch samples from ${sample_listing}\\n")
    
    result <- tryCatch({
        run_multiple_samples(
            file_listing = "${sample_listing}",
            featureN_PG = featureN_PG,
            minReadCnt = ${min_read_cnt},
            minDepth = ${min_depth},
            mafmin = ${maf_min},
            mafmax = ${maf_max}
        )
    }, error = function(e) {
        cat("ERROR in batch processing:\\n")
        cat(conditionMessage(e), "\\n")
        list(df_sums = data.frame(error = conditionMessage(e)))
    })
    
    write.table(result\$df_sums,
                "${meta.id}_batch_classification.tsv",
                sep = "\\t",
                row.names = FALSE,
                quote = FALSE)
    
    if (!"error" %in% colnames(result\$df_sums)) {
        saveRDS(result, "${meta.id}_batch_details.rds")
    }
    
    version_info <- c(
        '"${task.process}":',
        paste0('    r-base: "', R.version.string, '"'),
        paste0('    deseq2: "', packageVersion("DESeq2"), '"'),
        paste0('    rphenograph: "', packageVersion("Rphenograph"), '"')
    )
    writeLines(version_info, "versions.yml")
    """
}

process MDALL_COUNT_MATRIX {
    tag "$meta.id"
    label 'process_low_time'
    label 'process_low_memory'
    label 'process_low_medium_cpu'
    
    //module 'R/4.2.0:DESeq2/1.38.0:SummarizedExperiment/1.28.0' //error command not found (module)

    input:
    tuple val(meta), path(count_matrix)
    path(reference_dir)

    output:
    tuple val(meta), path("${meta.id}_matrix_classification.tsv"), emit: classification
    tuple val(meta), path("${meta.id}_matrix_details.rds"), emit: details
    path "versions.yml", emit: versions

    script:
    def feature_n = task.ext.featureN_PG ?: "1000,1058"
    
    """
    #!/usr/bin/env Rscript
    
    setwd("${reference_dir}")
    
    suppressPackageStartupMessages({
        library(DESeq2)
        library(SummarizedExperiment)
        library(dplyr)
        library(umap)
        library(Rphenograph)
        library(e1071)
        library(vroom)
    })
    
    source("${reference_dir}/mdall_core_functions.R")
    
    load("${reference_dir}/obj_234_HTSeq.RData")
    load("${reference_dir}/obj_1821.RData")
    load("${reference_dir}/models_svm.RData")
    
    featureN_PG <- as.numeric(unlist(strsplit("${feature_n}", ",")))
    
    result <- tryCatch({
        run_countMatrix(
            file_countMatrix = "${count_matrix}",
            featureN_PG = featureN_PG
        )
    }, error = function(e) {
        list(df_sums = data.frame(error = conditionMessage(e)))
    })
    
    write.table(result\$df_sums,
                "${meta.id}_matrix_classification.tsv",
                sep = "\\t",
                row.names = FALSE,
                quote = FALSE)
    
    if (!"error" %in% colnames(result\$df_sums)) {
        saveRDS(result, "${meta.id}_matrix_details.rds")
    }
    
    version_info <- c(
        '"${task.process}":',
        paste0('    r-base: "', R.version.string, '"')
    )
    writeLines(version_info, "versions.yml")
    """
}

workflow MDALL {
    take:
    samples  // channel: [ val(meta), path(count), path(vcf), path(fc), path(cicero) ]
    reference_dir

    main:
    ch_versions = Channel.empty()
    
    MDALL_CLASSIFY(
        samples,
        reference_dir
    )
    ch_versions = ch_versions.mix(MDALL_CLASSIFY.out.versions)

    emit:
    classification = MDALL_CLASSIFY.out.classification
    details = MDALL_CLASSIFY.out.details
    versions = ch_versions
}