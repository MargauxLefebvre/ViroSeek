/*
 * Process: Asssembly with Spades
 */
process filter_assembly_bioawk{
    label 'bioawk'
    tag "$meta.id"
    
    publishDir "${params.outdir}/$meta.id", mode: 'copy'

    input:
        tuple val(meta), path(reads), path(dir_spades)

    output:
        tuple val(meta), path(reads), path("*.filtered.fasta"), emit: read_assembly
        tuple val(meta), path("*.filtered.fasta"), emit: assembly_only

    script:
        """
        # Filter the assembly
        bioawk -c fastx '{ if(length(\$seq) > ${params.length_seq}) { print ">"\$name; print \$seq }}' \
        assembly_spades/contigs.fasta > ${meta.id}_contigs.filtered.fasta
        """
}