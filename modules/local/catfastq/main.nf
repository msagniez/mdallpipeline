process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'

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