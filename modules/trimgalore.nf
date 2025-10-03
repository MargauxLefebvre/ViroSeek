/*
 * Process: Trimming by Trimgalore
 */
process trimming {
    label 'trimgalore'
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta.id), path("trimgalore/${sample_id}/*.fq.gz"), emit: trimmed

    script:
        
        def paired_opt = meta.paired ? "--paired" : ""
        def read_args = meta.paired ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
        """
        mkdir -p trimgalore/${meta.id}
        trim_galore ${paired_opt} --cores ${task.cpus} --gzip \
            --output_dir trimgalore/${sample_id} \
            ${read_args}
        """
}