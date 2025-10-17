/*
Here are described all processes related to fastp
fastp is a tool designed to provide fast all-in-one preprocessing for reads files.
See https://github.com/OpenGene/fastp
*/
 
 
process trim_fastp {
    label 'fastp'
    tag "${meta.id}"
   
    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*_trim.reads.gz"), emit: trim
        path("${meta.id}_fastp_report.html"), emit: trim_report
   
    script:

        // set input/output according to short_paired parameter
        def input = "-i ${reads[0]}" 
        def readsBase0 = ViroSeekUtils.getCleanName(reads[0])
        def output = "-o ${readsBase0}_trim.reads.gz" 
        if ( meta.paired ){
            def readsBase1 = ViroSeekUtils.getCleanName(reads[1])
            input = "-i ${reads[0]} -I ${reads[1]}"
            output = "-o ${readsBase0}_trim.reads.gz -O ${readsBase1}_trim.reads.gz"
        }

        """
        fastp $input \\
              $output \\
              --thread ${task.cpus} \\
              --html ${meta.id}_fastp_report.html
        """
}