#!/usr/bin/env nextflow
 
nextflow.enable.dsl=2

/*************************************************
/ STEP 0 - parameters
/*************************************************/

/* Default value in case removed from nextflow.config */
input          = ''        // Path to a CSV file with 3 columns
outdir         = 'results' // Path to output results

// Database
    conta_ref      = 'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz,https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz'         // Path to the contamination reference file for BBduk (fasta format). By default, non-viral rRNA references from the SILVA rRNA database, release 138
    taxonkit_dir   = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'         // Path to the uncompressed taxdump folder (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
    diamond_db     = 'https://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz'         // Path to the Diamond database file (ncbi-nr.taxonomy.dmnd). Made using nr_db, prot_accession and taxdump
    prot_accession = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz' // URL to prot.accession2taxid.txt from NCBI, used to map accession numbers to Taxonomy IDs.

// Trim
trim           = ''          // Tools to use to trim [trimgalore, fastp]. Empty is accepted and means no trimming.
trim_opt       = ''          // option to set to trimming tool chosen

// Assembly
length_seq     = '0'         // Minimum contig length to keep after assembly (default 0, no filtering)

// Diamond options
diam_sensi     = ''          // Sensitivity option for Diamond (default empty). Options:  --faster, --fast, --mid-sensitive, --sensitive, --more-sensitive, --very-sensitive, --ultra-sensitive
diam_evalue    = '0.001'     // E-value threshold for Diamond (default 0.001)
diam_id        = '0'         // Minimum percentage identity for Diamond (default 0)
diam_querycov  = '0'  

// Other
params.skip_dedup = false
trimming_tools = [ 'trimgalore', 'fastp' ]
valid_sensi_flags = [ '--faster', '--fast', '--mid-sensitive', '--sensitive', '--more-sensitive', '--very-sensitive', '--ultra-sensitive' ]
params.help = null
params.debug = false

/*************************************************
/ STEP 0 HELP
/*************************************************/

println header()
if (params.help) { exit 0, helpMSG() }

/*************************************************
/ STEP 1 - check parameters
/*************************************************/

// check trimming tool
if (params.trim && !trimming_tools.contains(params.trim)) {
    exit 1, "Error: trimming tool ${params.trim} not recognized! Choose among this list: ${trimming_tools}\n"
}

// chck sensi flag for diamond
if (params.diam_sensi && !valid_sensi_flags.contains(params.diam_sensi)) {
    exit 1, "Error: DIAMOND sensi flag (--diam_sensi) ${params.diam_sensi} not recognized! Choose among this list: ${valid_sensi_flags} or let it empty!\n"
}


// Parameter message
log.info """

General Parameters
    input                      : ${params.input}
    outdir                     : ${params.outdir}

Contamination Filtering Parameters
    conta_ref                  : ${params.conta_ref ?: "default SILVA rRNA database (release 138)"}

Assembly Parameters
    length_seq                 : ${params.length_seq}
    skip_dedup                 : ${params.skip_dedup ? "yes" : "no"}

Taxonomic Assignment Parameters
    diamond_db                 : ${params.diamond_db ?: "default NCBI RefSeq viral protein DB"}
    taxonkit_dir               : ${params.taxonkit_dir ?: "default (auto-downloaded)"}
    diam_sensi                 : ${params.diam_sensi ?: "none (default)"}
    diam_evalue                : ${params.diam_evalue}
    diam_id                    : ${params.diam_id}
    diam_querycov              : ${params.diam_querycov}

Trimming Parameters
    trim                       : ${params.trim ?: "none (skipped)"}
    trim_opt                   : ${params.trim_opt ?: "none"}

"""

/*************************************************
/ STEP 2 - Include needed modules
/*************************************************/
include {prepare_taxonkit; taxo_quanti} from "$baseDir/modules/bash.nf"
include {filtering_bbduk} from "$baseDir/modules/bbmap.nf"
include {filter_assembly_bioawk} from "$baseDir/modules/bioawk.nf"
include {taxo_assign_diamond} from "$baseDir/modules/diamond.nf"
include {trim_fastp} from "$baseDir/modules/fastp.nf"
include {fastqc as fastqc_raw; fastqc as fastqc_posttrim; fastqc as fastqc_postBBduk; fastqc as fastqc_align; fastqc as fastqc_dedup;} from "$baseDir/modules/fastqc.nf"
include {quantification_minimap2} from "$baseDir/modules/minimap2.nf"
include {samtools_sam2sortedbam; samtools_markdup; samtools_idxstats} from "$baseDir/modules/samtools.nf"
include {assembly_spades} from "$baseDir/modules/spades.nf"
include {taxo_table_taxonkit} from "$baseDir/modules/taxonkit.nf"
include {trim_trimgalore} from "$baseDir/modules/trimgalore.nf"
include {prepare_silva_DB; prepare_silva_DB_list; dwnload_diamond_DB; dwnload_taxonkit_DB; prepare_viral_accessions; prepare_accession2taxid; subset_diamond_DB; prepare_diamond_DB} from "$baseDir/modules/database.nf"

/*************************************************
/ STEP 4 - MAIN WORKFLOW
/*************************************************/
workflow {

    // Initialize channels
    Channel.empty().set{logs}

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

// -------------------------------- WORKFLOW STEPS --------------------------------
    // ---- DATABASES
    
    // ---- check conta sequence
    if (params.conta_ref) {

        // Parse params
        def conta_list = params.conta_ref instanceof List ?
            params.conta_ref.collect { it.toString().trim() } :
            params.conta_ref.toString().tokenize(',').collect { it.trim() }

        // Validate
        def invalid_files = conta_list.findAll { ref ->
            !ViroSeekUtils.is_url(ref) &&
            !ViroSeekUtils.resolveFile(ref, "conta_ref")
        }

        if (invalid_files) {
            exit 1, "Error: contamination sequence file(s) not found:\n${invalid_files.join('\n')}"
        }

        // Local files are used directly and remote URLs are downloaded.
        def conta_ch = Channel.fromList(
            conta_list.collect { ref -> file(ref) }
        )

        // Group all references into one list
        conta_ch
            .collect()
            .set { conta_files_ch }

        // Concatenate all contamination references
        silva_ref = prepare_silva_DB(conta_files_ch)
    }


// ---- check DIAMOND resources only if BOTH are provided
    if (params.diamond_db && params.taxonkit_dir) {
      def found = false
      
      if (ViroSeekUtils.is_url(params.diamond_db) && ViroSeekUtils.is_url(params.taxonkit_dir) && ViroSeekUtils.is_url(params.prot_accession)){
        found = true
        log.info "No Diamond database and/or Taxonkit directory provided, using NCBI Ref-Seq non-redundant protein database restricted to viral sequences"
        dwnload_diamond_DB()
        dwnload_taxonkit_DB()
        prepare_viral_accessions(dwnload_taxonkit_DB.out.taxonkit_dir)
        prepare_accession2taxid(prepare_viral_accessions.out.viral_taxids)
        subset_diamond_DB(prepare_accession2taxid.out.viral_accessions, dwnload_diamond_DB.out.nr_database)
        prepare_diamond_DB(prepare_accession2taxid.out.accession2taxidFULL, subset_diamond_DB.out.viral_ref, dwnload_taxonkit_DB.out.taxonkit_dir)
        diamond_db   = prepare_diamond_DB.out.diamond_db
        taxonkit_dir = dwnload_taxonkit_DB.out.taxonkit_dir}
        
      else if (ViroSeekUtils.resolveFile(params.diamond_db, "diamond_db")){
        found = true
        log.info "Use the datatbase provided by the user."
        diamond_db   = Channel.fromPath(params.diamond_db, checkIfExists: true)
        taxonkit_dir = Channel.fromPath(params.taxonkit_dir, checkIfExists: true)}
      
      if (!found) {exit 1, "Error: no database provided or incomplete.\n Please check input files or URL provided.\n"}
    }

// check trimming tool
if (params.trim && !trimming_tools.contains(params.trim)) {
    exit 1, "Error: trimming tool ${params.trim} not recognized! Choose among this list: ${trimming_tools}\n"
}
    
    
    // Run FastQC pre-trimming
    fastqc_raw(samples_ch, "fastQC/fastqc_pretrim", "raw")
    logs.concat(fastqc_raw.out).set{logs} // save log

    // Run trimming and capture output channel
    if (params.trim && params.trim == 'trimgalore') {
        trimmed_ch = samples_ch | trim_trimgalore
    } else if (params.trim && params.trim == 'fastp') {
        trimmed_ch = samples_ch | trim_fastp
        logs.concat(trimmed_ch.trim_report).set{logs} // save log
    } else {
        log.info "No trimming selected, continuing without trimming step."
        trimmed_ch = samples_ch
    }

    // FastQC post trimming using trimmed reads
    fastqc_posttrim(trimmed_ch.trim, "fastQC/fastqc_posttrim", "posttrim")
    logs.concat(fastqc_posttrim.out).set{logs} // save log

    // BBduk filtering using trimmed reads
    BBduk_ch = trimmed_ch.trim.combine(silva_ref) | filtering_bbduk

    // FastQC post BBduk
    fastqc_postBBduk(BBduk_ch, "fastQC/fastqc_postBBduk", "postBBduk")
    logs.concat(fastqc_postBBduk.out).set{logs} // save log

    //Make the assembly
    assembly_ch = BBduk_ch | assembly_spades | filter_assembly_bioawk

    // Quantification
    quantification_minimap2(assembly_ch.read_assembly)
    samtools_sam2sortedbam(quantification_minimap2.out.sam_file)
    fastqc_align(samtools_sam2sortedbam.out.bam_file, "fastQC/fastqc_align", "align")
    logs.concat(fastqc_align.out).set{logs} // save log
    // Deduplicatiom option
    if( !params.skip_dedup ) {
    samtools_markdup(samtools_sam2sortedbam.out.bam_file)
    fastqc_dedup(samtools_markdup.out.bam_dedup_file, "fastQC/fastqc_dedup", "dedup")
    logs.concat(fastqc_dedup.out).set{logs} // save log
    samtools_idxstats(samtools_markdup.out.bam_dedup_file)
    } else {
      samtools_idxstats(samtools_sam2sortedbam.out.bam_file)
    }

    //Taxonomic assignation
    diamond_ch = taxo_assign_diamond(assembly_ch.assembly_only, diamond_db.collect())
    prepare_taxonkit(diamond_ch)
    taxo_table_taxonkit(prepare_taxonkit.out.taxid, taxonkit_dir.collect())

    taxo_table_taxonkit.out.taxo_table.map { meta, accession_taxid_txt, taxo_table_txt -> tuple(meta.id, meta, accession_taxid_txt, taxo_table_txt) }
        .join(diamond_ch.map { meta, diamond_file_tsv -> tuple(meta.id, meta, diamond_file_tsv) })
        .join(samtools_idxstats.out.idxstats.map { meta, idxstats_file_tsv -> tuple(meta.id, meta, idxstats_file_tsv) })
        .map { id, meta, accession_taxid_txt, taxo_table_txt, meta2, diamond_file_tsv, meta3, idxstats_file_tsv -> tuple(meta, idxstats_file_tsv, diamond_file_tsv,  accession_taxid_txt, taxo_table_txt) }
        .set { jointaxoquanti_ch }
    
    jointaxoquanti_ch | taxo_quanti
    
/* FYI
idxstats_file_tsv = quant_stats
diamond_file_tsv = diamond_file    
taxid = *.accession_taxid.txt"
taxo_table_txt=     *_taxonomy_table.txt"
*/
}

/*************************************************
/               EXTRA - FUCTIONS
/*************************************************/

// Function to print the header
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
            --input                 Path to a CSV file that expects 3 columns (with a header): `sample,read1,read2.`. Each row represents a paired-end sample.

    Optional parameters
        General
            --outdir                Path to output results directory. Default: results

        Trimming
            --trim                  Tools to use to trim among this list ${trimming_tools} [trimgalore, fastp]. Empty is accepted and means no trimming.
            --trim_opt              Option to be used by the trimming tool chosen

        Assembly
            --length_seq            Minimum contig length to keep after assembly (default 0, no filtering)
            
        Deduplication
            --skip_dedup            Skip deduplication step (samtools markdup), default: false
    
        Database
            --conta_ref             Path to the contamination reference file for BBduk (fasta format)
            --diamond_db            Path to the DIAMOND database file (ncbi-nr.taxonomy.dmnd)
            --taxonkit_dir          Path to the TaxonKit files (download and uncompress from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)

        Taxonomic assignation
            --diam_sensi            Sensitivity mode for DIAMOND among this list ${valid_sensi_flags} (default: empty, which means no sensitivity flag)
            --diam_evalue           E-value threshold for DIAMOND (default: 1e-5)
            --diam_id               Minimum percent identity for DIAMOND (default: 0)
            --diam_querycov         Minimum query coverage for DIAMOND (default: 0)
  
        Other
            --monochrome_logs        Set to true to disable color in logs (default: false)

    """
}

// When the workflow is complete, print a message
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
