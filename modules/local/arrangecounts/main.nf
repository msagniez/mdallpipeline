process ARRANGE_COUNTS {
    tag "$meta.id"
    label 'process_low'
    label 'process_low_cpu'

    input:
    tuple val(meta), path(count_table)

    output:
    tuple val(meta), path("${meta.id}.featureCounts_arranged.tsv"), emit: arrcounts

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract first column (Geneid) and last column (read count)
    # Skip comment lines and header, then select columns
    awk 'BEGIN {OFS="\\t"} 
         /^#/ {next} 
         NR==1 {print \$1, \$NF; next} 
         {print \$1, \$NF}' ${count_table} > ${prefix}.featureCounts_arranged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/^.*Awk //; s/,.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts_arranged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/^.*Awk //; s/,.*\$//')
    END_VERSIONS
    """
}