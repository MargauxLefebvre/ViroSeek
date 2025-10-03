#!/usr/bin/env nextflow
 
nextflow.enable.dsl=2

//*************************************************
// STEP 0 - parameters
//*************************************************    

/* Default value in case removed from nextflow.config */
input          = ''        // Path to a CSV file with 3 columns
outdir         = 'results' // Path to output results

// Database
silva_ref      = ''         // Path to the SILVA reference file for BBduk (fasta format)
diamond_db     = ''         // Path to the Diamond database file (ncbi-nr.taxonomy.dmnd)
prot_accession = ''         // Path to the prot.accession2taxid.txt file from NCBI (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz)
taxonkit_dir      = ''         // Path to the TaxonKit files (download and uncompress from https://bioinf.shenwei.me/taxonkit/)

// Trim
trimming       = ''          // Tools to use to trim [trimgalore, fastp]. Empty is accepted and means no trimming.
trimming_opt   = ''          // option to set to trimming tool choosen

// Assembly
length_seq     = '0'         // Minimum contig length to keep after assembly (default 0, no filtering)

// Diamond options
diam_sensi     = ''          // Sensitivity option for Diamond (default empty). Options: --sensitive, --more-sensitive, --very-sensitive
diam_evalue    = '0.001'     // E-value threshold for Diamond (default 0.001)
diam_id        = '0'         // Minimum percentage identity for Diamond (default 0)
diam_querycov  = '0'  

// Other
trimming_tools = [ 'trimgalore', 'fastp' ]
valid_sensi_flags = ["--sensitive", "--more-sensitive", "--very-sensitive"]
params.help = null
params.debug = false

//*************************************************
// STEP 1 - check parameters
//*************************************************

// ---- check SILVA file
def silva_ref_file        = ViroSeekUtils.resolveFile(params.silva_ref, "silva_ref")
def diamond_db_file       = ViroSeekUtils.resolveFile(params.diamond_db, "diamond_db")
def prot_accession_file   = ViroSeekUtils.resolveFile(params.prot_accession, "prot_accession")
def taxonkit_dir_file     = ViroSeekUtils.resolveFile(params.taxonkit_dir, "taxonkit_dir")
def nr_db_file            = ViroSeekUtils.resolveFile(params.nr_db, "nr_db")     

//*************************************************
// STEP 2 - Include needed modules
//*************************************************
include {build_diamond_db} from "$baseDir/modules/diamond.nf"
include {fastqc as fastqc_raw; fastqc as fastqc_posttrim} from "$baseDir/modules/fastqc.nf"
include {trimming} from "$baseDir/modules/trimgalore.nf"

/*
 * Process: Filtration with BBduk
 */
process bbduk_filtering {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(trimmed_read1), path(trimmed_read2)

    output:
    tuple val(sample_id), 
          path("${sample_id}_R1_rmrdna.fastq.gz"), 
          path("${sample_id}_R2_rmrdna.fastq.gz")

    script:
    """
    bbduk.sh \
        -Xmx34g \
        in1=${trimmed_read1} \
        in2=${trimmed_read2} \
        out1=${sample_id}_R1_rmrdna.fastq.gz \
        out2=${sample_id}_R2_rmrdna.fastq.gz \
        threads=8 \
        ref=${params.silva_ref}
    """
}

/*
 * Process: Merge the fastq for Spades
 */
process merging {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(postBBduk_read1), path(postBBduk_read2)

    output:
        tuple val(sample_id), 
          path("${sample_id}_inter_rmrdna.fastq.gz")

    script:
    """
    seqtk mergepe ${postBBduk_read1} ${postBBduk_read2} | bgzip > ${sample_id}_inter_rmrdna.fastq.gz
    """
}

/*
 * Process: FastQC after BBduk
 */
process fastqc_postBBduk {

    tag "$sample_id"
    
    publishDir "${params.outdir}/$sample_id/fastQC/fastqc_postBBduk", mode: 'copy'

    input:
    tuple val(sample_id), path(postBBduk_read1), path(postBBduk_read2)

    output:
    path("${sample_id}_fastqc_postbbduk")

    script:
    """
    mkdir -p ${sample_id}_fastqc_postbbduk
    zcat ${postBBduk_read1} ${postBBduk_read2} | fastqc -t 4 stdin:${sample_id}.trim --outdir=${sample_id}_fastqc_postbbduk
    """
}

/*
 * Process: Asssembly with Spades
 */
process assembly {

    tag "$sample_id"
    
    publishDir "${params.outdir}/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(postBBduk_inter)

    output:
        tuple val(sample_id), 
          path("assembly_spades")

    script:
    """
    mkdir -p assembly_spades
    
    ${params.spadesbin}/spades.py --rnaviral \
    --12 ${postBBduk_inter} \
    --threads 12 \
    --memory 72 \
    -o assembly_spades
    
    # Filter the assembly
    bioawk -c fastx '{ if(length(\$seq) > ${params.length_seq}) { print ">"\$name; print \$seq }}' \
    assembly_spades/contigs.fasta > assembly_spades/contigs.filtered.fasta
    """
}

/*
 * Process: Quantification by mapping cleaned reads to the assembly
 */
process quantification {

    tag "$sample_id"
    
    publishDir "${params.outdir}/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(postBBduk_read1), path(postBBduk_read2), path(dir_spades)

    output:
        tuple val(sample_id), 
          path("${sample_id}_contigs_reads.tsv")

    script:
    """
    minimap2 -t 16 -a ${dir_spades}/contigs.filtered.fasta \
            ${postBBduk_read1} ${postBBduk_read2} \
            | samtools view -@ 16 -S -b - > ${sample_id}_map.bam

    samtools sort -@ 16 ${sample_id}_map.bam -o ${sample_id}_sorted.bam
    samtools markdup -@ 16 -r ${sample_id}_sorted.bam ${sample_id}_dedup.bam
    samtools index ${sample_id}_dedup.bam
    samtools idxstats ${sample_id}_dedup.bam > ${sample_id}_contigs_reads.tsv
    """
}

/*
 * Process: Taxonomic assignation with diamond
 */
process taxo_assign {

    tag "$sample_id"

    publishDir "${params.outdir}/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(dir_spades)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")

    script:
    // Safely build sensitivity flag only if valid
    def valid_sensi_flags = [
        '--faster', '--fast', '--mid-sensitive',
        '--sensitive', '--more-sensitive', '--very-sensitive', '--ultra-sensitive'
    ]

    def sensi_flag = valid_sensi_flags.contains(params.diam_sensi) ? params.diam_sensi : ''
    """
    diamond blastx -p 32 -d ${params.diamond_db} -q ${dir_spades}/contigs.filtered.fasta \
        -o ${sample_id}.tsv --max-target-seqs 1 -e ${params.diam_evalue} --id ${params.diam_id} \
        --query-cover ${params.diam_querycov} --range-culling -F 15 ${sensi_flag}
    """
}


/*
 * Process: Generate the taxonomic table
 */
process taxo_table {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(output_diamond)

    output:
        tuple val(sample_id),
        path("${sample_id}.accession_taxid.txt"),
          path("${sample_id}_taxonomy_table.txt")

    script:
    """
  # Retrieve the accession IDs column from the Diamond file
      cut -f2 ${output_diamond} > ${sample_id}_accession.txt
  
  # Associate taxIDs with accession IDs
      zgrep -F -f ${sample_id}_accession.txt ${params.prot_accession} > ${sample_id}.accession_taxid.txt
  
  # Extract taxIDs
      cut -f2 ${sample_id}.accession_taxid.txt > ${sample_id}.taxid.txt
  
  # Use TaxonKit to obtain taxonomy from taxIDs
      taxonkit lineage ${sample_id}.taxid.txt --data-dir ${params.taxonkit_dir} > ${sample_id}_taxonomy_table.txt
    """
}

/*
 * Process: Add the quantification to the taxo table
 */
process taxo_quanti {

    tag "$sample_id"
    
    publishDir "${params.outdir}/$sample_id", mode: 'copy'

    input:
    tuple val(sample_id), path(quanti_stats), path(output_diamond), path(taxid), path(taxo_table)

    output:
        tuple val(sample_id), 
          path("${sample_id}_taxonomy_viral.clean.txt"), 
          path("${sample_id}_all_viral_taxonomy.txt")

    script:
    """
    # Merge contigs_reads.tsv and X.tsv files by contig number
        join -1 1 -2 1 <(sort "${quanti_stats}") <(sort "${output_diamond}") > ${sample_id}_contig_reads_accession.txt
    
    # Add taxIDs
        join -1 5 -2 1 <(sort -k5,5 "${sample_id}_contig_reads_accession.txt") <(sort -k1,1 "${taxid}") > ${sample_id}_contig_reads_accession_taxid.txt
    
    # Merge with taxonomy table
        join -1 16 -2 1 <(sort -k16,16 "${sample_id}_contig_reads_accession_taxid.txt") <(sort -k1,1 "${taxo_table}") > ${sample_id}_final_assembly.txt
    
    # Filter to keep only lines containing “virus”.
        grep -i "virus" ${sample_id}_final_assembly.txt > ${sample_id}_virus_taxonomy.txt
    uniq ${sample_id}_virus_taxonomy.txt > ${sample_id}_filter_viral_taxonomy.txt
    
    # Keep only reads count and the taxonomy in a tab delimited file
      awk '{printf "%s\t", \$5; for (i=17; i<=NF; i++) printf "%s%s", \$i, (i<NF?" ":"\\n")}' ${sample_id}_filter_viral_taxonomy.txt > ${sample_id}_taxonomy_viral.txt
     sed 's/;/\t/g' ${sample_id}_taxonomy_viral.txt > ${sample_id}_taxonomy_viral.clean.txt
     cp -r ${sample_id}_filter_viral_taxonomy.txt ${sample_id}_all_viral_taxonomy.txt
    """
}

//*************************************************
// STEP 1 - HELP
//*************************************************

println header()
if (params.help) { exit 0, helpMSG() }

/*************************************************
/ STEP 1 - PARAMS CHECK
/************************************************/

/*
 * Main workflow
 */
workflow {

    // -------------------------------- INPUT --------------------------------
    // ---- CSV INPUT FILE
    samples_ch = Channel.fromPath(params.input)
            .splitCsv(header: true, sep: ',')
            .map { row ->

                // Check sample column
                if ( row.sample == null ){ 
                    error "The input ${input} file does not contain a 'sample' column!\n" 
                } 
                def sample_id    = row.sample

                // Check input1/fastq1/read1 column
                def fastq1 = row.input1?.trim() ?: row.fastq1?.trim() ?: row.read1?.trim()
                if (!fastq1)
                    error "The input ${params.input} file does not contain a 'input1' or 'fastq1' or 'read1' column!\n"
                else
                    fastq1 = file(fastq1)   
                if (!ViroSeekUtils.is_url(fastq1) && !fastq1.exists())
                    error "The input ${fastq1} file does not exist!\n"

                if (ViroSeekUtils.is_url(fastq1))
                    log.info "This fastq input is an URL: ${fastq1}"
                
                // Check input2/fastq2/read2 column
                def paired = false;
                def fastq2 = row.input2?.trim() ?: row.fastq2?.trim() ?: row.read2?.trim()
                fastq2 = file(fastq2) 
                if (fastq2) {
                    paired = true
                    if (!ViroSeekUtils.is_url(fastq2) && !fastq2.exists())
                        error "The input ${fastq2} file does not exist!\n"
                    if (ViroSeekUtils.is_url(fastq2))
                        log.info "This fastq input is an URL: ${fastq2}"
                }

                // Create a tuple with metadata and reads
                def meta = [ id: sample_id, paired: paired ]
                def reads = paired ? [fastq1, fastq2] : fastq1

                // Return only if the fastq file(s) extension are valid
                if ( ViroSeekUtils.is_fastq(fastq1) && (!fastq2 || ViroSeekUtils.is_fastq(fastq2)) )
                    return tuple(meta, reads)
                else
                    error "File(s) for sample ${sample_id} do not look like FASTQ"
            }
    // ---- DATABASES
    silva_ref = Channel.fromPath(params.silva_ref, checkIfExists: true)
                        .ifEmpty { exit 1, "Cannot find silva_ref file matching ${params.silva_ref}!\n" }
    if (diamond_db_file ){
        diamond_db = Channel.fromPath(params.diamond_db, checkIfExists: true)
                            .ifEmpty { exit 1, "Cannot find diamond_db matching ${params.diamond_db}!\n" }
    } else {
        log.info "Diamond_db not provided, ViroSeek will prepare it for you!"
        if (nr_db_file) {
            nr_db = Channel.fromPath(params.nr_db, checkIfExists: true)
                            .ifEmpty { exit 1, "Cannot find nr_db matching ${params.nr_db}!\n" }
        } else {
            exit 1, "Error: nr_db mandatory to prepare the DIAMOND database! You must provide either --diamond_db or --nr_db parameter!\n"
        }
    }
    prot_accession = Channel.fromPath(params.prot_accession, checkIfExists: true)
                        .ifEmpty { exit 1, "Cannot find prot_accession matching ${params.prot_accession}!\n" }
    taxonkit_dir = Channel.fromPath(params.taxonkit_dir, checkIfExists: true)
                        .ifEmpty { exit 1, "Cannot find taxonkit_dir matching ${params.taxonkit_dir}!\n" }



    // Run FastQC pre-trimming
    samples_ch.view()
    fastqc_raw(samples_ch, "fastQC/fastqc_pretrim", "raw")

    // Run trimming and capture output channel
   // trimmed_ch = samples_ch | trimming

/*
    // FastQC post trimming using trimmed reads
    trimmed_ch | fastqc_posttrim

    // BBduk filtering using trimmed reads
    BBduk_ch = trimmed_ch | bbduk_filtering
    
    // FastQC post BBduk
    BBduk_ch | fastqc_postBBduk
    
    //Make the assembly
    assembly_ch = BBduk_ch | merging | assembly
    
    //Map on the assembly
    joinfastqassembly_ch = BBduk_ch.join(assembly_ch, by: 0)
    quanti_ch = joinfastqassembly_ch | quantification
    
    //Taxonomic assignation
    diamond_ch = assembly_ch | taxo_assign 
    taxo_ch = diamond_ch | taxo_table
    
    jointaxoquanti_ch = quanti_ch
    .join(diamond_ch, by: 0)
    .map { sample_id, quant_stats, diamond_file -> tuple(sample_id, quant_stats, diamond_file) }
    .join(taxo_ch, by: 0)
    .map { sample_id, quant_stats, diamond_file, taxid, taxo_table -> tuple(sample_id, quant_stats, diamond_file, taxid, taxo_table) }

    jointaxoquanti_ch | taxo_quanti
    */
}

//*************************************************
def header(){
    // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";

    return """
    -${c_dim}--------------------------------------------------${c_reset}-
    ${c_blue}.-./`) ${c_white}.-------.    ${c_red} ______${c_reset}
    ${c_blue}\\ .-.')${c_white}|  _ _   \\  ${c_red} |    _ `''.${c_reset}     French National   
    ${c_blue}/ `-' \\${c_white}| ( ' )  |  ${c_red} | _ | ) _  \\${c_reset}    
    ${c_blue} `-'`\"`${c_white}|(_ o _) /  ${c_red} |( ''_'  ) |${c_reset}    Research Institute for    
    ${c_blue} .---. ${c_white}| (_,_).' __ ${c_red}| . (_) `. |${c_reset}
    ${c_blue} |   | ${c_white}|  |\\ \\  |  |${c_red}|(_    ._) '${c_reset}    Sustainable Development
    ${c_blue} |   | ${c_white}|  | \\ `'   /${c_red}|  (_.\\.' /${c_reset}
    ${c_blue} |   | ${c_white}|  |  \\    / ${c_red}|       .'${c_reset}
    ${c_blue} '---' ${c_white}''-'   `'-'  ${c_red}'-----'`${c_reset}
    ${c_purple}ViroSeek - Viral detection pipeline for second-generation sequencing - v${workflow.manifest.version}${c_reset}
    ${c_white}${workflow.manifest.author}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

// Help Message
def helpMSG() {
    log.info """
    ViroSeek - Viral detection pipeline for second-generation sequencing - v${workflow.manifest.version}

        ViroSeek is a pipeline designed for the analysis of target-enriched
        libraries, with specific optimization for managing high PCR duplicate
        rates and performing per-sample assembly. The primary inputs are
        paired-end FASTQ files generated by next-generation sequencing (NGS),
        while external databases for ribosomal RNA filtering and taxonomic
        assignment must be provided by the user.

        Usage example:
        nextflow run ViroSeek.nf --input /path/to/file.csv --trimming trimgalore -profile singularity

        --help                      prints the help section

    Mandatory Parameters
        General
            --input                 Path to a CSV file that expects 3 columns: `sample_id,read1.fastq.gz,read2.fastq.gz`. Each row represents a paired-end sample.
    
        Database
            --silva_ref             Path to the SILVA reference file for BBduk (fasta format)
            --diamond_db            Path to the DIAMOND database file (ncbi-nr.taxonomy.dmnd)
            --prot_accession        Path to the prot.accession2taxid.txt file from NCBI (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz)
            --taxonkit_dir          Path to the TaxonKit files (download and uncompress from https://bioinf.shenwei.me/taxonkit/)

    Optional parameters
        General
            --outdir                Path to output results directory. Default: results

        Trimming
            --trimming              Tools to use to trim among this list ${trimming_tools} [trimgalore, fastp]. Empty is accepted and means no trimming.
            --trimming_opt          option to be used by the trimming tool choosen

        Assembly
            --length_seq             Minimum contig length to keep after assembly (default 0, no filtering)

        Taxonomic assignation
            --diam_sensi            Sensitivity mode for DIAMOND among this list ${valid_sensi_flags} (default: empty, which means no sensitivity flag)
            --diam_evalue           E-value threshold for DIAMOND (default: 1e-5)
            --diam_id               Minimum percent identity for DIAMOND (default: 0)
            --diam_querycov         Minimum query coverage for DIAMOND (default: 0)
  
        Other
            --monochrome_logs        Set to true to disable color in logs (default: false)

    """
}


/**************         onComplete         ***************/

workflow.onComplete {

    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.success) {
        log.info "\n${c_green}    ViroSeek pipeline complete!${c_reset}"
    } else {
        log.error "${c_red}Oops .. something went wrong${c_reset}"
    }

    log.info "    The results are available in the ‘${params.outdir}’ directory."
    log.info """
    ViroSeek Pipeline execution summary
    --------------------------------------
    Completed at : ${workflow.complete}
    UUID         : ${workflow.sessionId}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Exit Status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """
}
