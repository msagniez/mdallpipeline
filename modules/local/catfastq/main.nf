process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'
    label "process_low_time"       // allow 2h
    label "process_medium_high_memory" //increased memory for drac memory alloc
    label "process_low_cpu"         // Label for mpgi drac 1 cpu usage

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("${meta.id}_concatenated.fastq.gz"), emit: fastq

    script:
    def is_compressed = fastq_files[0].name.endsWith('.gz')
    
    if (is_compressed) {
        """
        cat ${fastq_files} > ${meta.id}_concatenated.fastq.gz
        """
    } else {
        """
        cat ${fastq_files} | gzip > ${meta.id}_concatenated.fastq.gz
        """
    }
}