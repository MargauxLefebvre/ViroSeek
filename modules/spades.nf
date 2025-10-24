/*
 * Process: Asssembly with Spades
 */
process assembly_spades {
    label 'spades'
    tag "${meta.id}"
    
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path(reads), path("assembly_spades")

    script:
        def read_args_in = meta.paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
        """
        mkdir -p assembly_spades
        
        spades.py \
        --rnaviral \
        ${read_args_in} \
        --threads ${task.cpus} \
        --memory ${task.memory.toGiga()} \
        -o assembly_spades
        """
}