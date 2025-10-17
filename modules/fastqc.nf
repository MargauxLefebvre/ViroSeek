/*
 * Process: FastQC before trimming
 */
process fastqc {
    label 'fastqc'
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/${outpath}", mode: 'copy'

    input:
        tuple val(meta), path(reads)
        val outpath
        val suffix

    output:
        path ("*logs")


    script:

        def add_suffix = suffix ? "_${suffix}_" : '_'

        """
        mkdir ${meta.id}${add_suffix}logs
        fastqc -t ${task.cpus} -o ${meta.id}${add_suffix}logs -q ${reads}
        """
}