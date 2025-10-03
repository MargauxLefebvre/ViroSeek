/*
 * Process: Build DIAMOND database
 */
process build_diamond_db {
    label 'diamond'
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/fastQC/fastqc_pretrim", mode: 'copy'

    input:
        tuple val(meta), path(reads)


    output:
        path ("*logs")


    script:

        def sample_id = meta.id

        """
        mkdir -p ${sample_id}_fastqc_pre
        fastqc -t ${task.cpus} -q ${reads} -o ${sample_id}_fastqc_pre
        """
}