/*
 * Process: Filtration with BBduk
 */
process filtering_bbduk {
    label 'bbmap'
    tag "$meta.id"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*_rmrdna.fastq.gz")

    script:

        def read_args_in = meta.paired ? "in1=${reads[0]} in2=${reads[1]}" : "in=${reads[0]}"
        def read_args_out = meta.paired ? "out1=${meta.id}_R1_rmrdna.fastq.gz out2=${meta.id}_R2_rmrdna.fastq.gz" : "out=${meta.id}_rmrdna.fastq.gz"
        """
        bbduk.sh \
            -Xmx${task.memory.toGiga()}g \
            ${read_args_in} \
            ${read_args_out} \
            threads=${task.cpus} \
            ref=${params.silva_ref}
        """
}