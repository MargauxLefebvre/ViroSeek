/*
 * Process: Taxonomic assignation with diamond
 */
process taxo_assign_diamond {
    tag "$meta.id"
    label 'diamond'

    publishDir "${params.outdir}/$meta.id", mode: 'copy'

    input:
        tuple val(meta), path(assemby)
        path(diamond_db)

    output:
        tuple val(meta), path("*.tsv")

    script:

    """
    diamond blastx -p ${task.cpus} -d ${diamond_db} -q ${assemby} \
        -o ${meta.id}_diamond_blastx.tsv --max-target-seqs 1 -e ${params.diam_evalue} --id ${params.diam_id} \
        --query-cover ${params.diam_querycov} --range-culling -F 15 ${params.diam_sensi}
    """
}
