/*
 * Process: Quantification by mapping cleaned reads to the assembly
 */
process samtools_sam2sortedbam {
    label 'samtools'
    tag "$meta.id"

    input:
        tuple val(meta), path(sam_file)

    output:
        tuple val(meta), path("*_sorted.bam"), emit: bam_file

    script:
    """
    samtools view -@ ${task.cpus} -S -b ${sam_file} \
      | samtools sort -@ ${task.cpus} -o ${meta.id}_sorted.bam
    """
}

process samtools_markdup {
    label 'samtools'
    tag "$meta.id"

    input:
        tuple val(meta), path(bam_file)

    output:
        tuple val(meta), path("*_dedup.bam"), emit: bam_dedup_file

    script:
    """
    samtools markdup -@ ${task.cpus} -r ${bam_file} ${meta.id}_dedup.bam
    """
}

process samtools_idxstats {
    label 'samtools'
    tag "$meta.id"
    
    publishDir "${params.outdir}/$meta.id", mode: 'copy'

    input:
        tuple val(meta), path(bam_file)

    output:
        tuple val(meta), path("*_contigs_reads.tsv"), emit: idxstats

    script:
    """
    samtools index ${bam_file}
    samtools idxstats ${bam_file} > ${meta.id}_contigs_reads.tsv
    """
}