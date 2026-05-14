/*
 * Process: Prepare the taxonkit input file
 */
process prepare_taxonkit {
    label 'bash'
    tag "$meta.id"

    input:
        tuple val(meta), path(output_diamond)

    output:
        tuple val(meta), path("*_taxid.txt"), emit: taxid

    script:
        """
        # Retrieve the accession IDs column from the Diamond file
        # and Associate taxIDs with accession IDs
        # Remove duplicates
        cut -f13 ${output_diamond} | sort -u > ${meta.id}_taxid.txt
        """
}

/*
 * Process: Add the quantification to the taxo table
 * !!! join remove lines without correspondance in both files !!! join -a 1 to keep all lines from the first file
 */
process taxo_quanti {
    label 'bash'
    tag "$meta.id"
    
    publishDir "${params.outdir}/$meta.id", mode: 'copy'

    input:
    tuple val(meta), path(quanti_stats), path(output_diamond), path(taxid), path(taxo_table)

    output:
        tuple val(meta), 
          path("*_taxonomy_viral.clean.txt"), 
          path("*_all_viral_taxonomy.txt")

    script:
    """
        # Merge contigs_reads.tsv and X.tsv files by contig number
        join -1 1 -2 1 <(sort "${quanti_stats}") <(sort "${output_diamond}") > ${meta.id}_contig_reads_accession.txt
    
        # Merge with taxonomy table
        join -1 16 -2 1 <(sort -k16,16 "${meta.id}_contig_reads_accession.txt") <(sort -k1,1 "${taxo_table}") > ${meta.id}_final_assembly.txt
    meta.id=SIM
        # Filter to keep only lines containing “virus”.
        grep -i "virus" ${meta.id}_final_assembly.txt > ${meta.id}_virus_taxonomy.txt
        uniq ${meta.id}_virus_taxonomy.txt > ${meta.id}_filter_viral_taxonomy.txt
    
        # Keep only reads count and the taxonomy in a tab delimited file
        awk '{printf "%s\t", \$4; for (i=17; i<=NF; i++) printf "%s%s", \$i, (i<NF?" ":"\\n")}' ${meta.id}_filter_viral_taxonomy.txt > ${meta.id}_taxonomy_viral.txt
        sed 's/;/\t/g' ${meta.id}_taxonomy_viral.txt > ${meta.id}_taxonomy_viral.temp.txt
 
        #Put clear headers
        cp -r ${meta.id}_filter_viral_taxonomy.txt ${meta.id}_all_viral_taxonomy.txt
    """
}

