/*
 * Process: Quantification by mapping cleaned reads to the assembly
 */
process quantification_minimap2 {
    label 'minimap2'
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads), path(assembly)

    output:
        tuple val(meta), path("*.sam"), emit: sam_file

    script:
        def read_args_in = meta.paired ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
        """
        minimap2 -t ${task.cpus} \
                -a ${assembly} \
                ${read_args_in} \
                > ${meta.id}.sam
        """
}
