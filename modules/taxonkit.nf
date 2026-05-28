/*
 * Process: Generate the taxonomic table
 */
process taxo_table_taxonkit {
    label 'taxonkit'
    tag "$meta.id"

    input:
        tuple val(meta), path(accession_taxid)
        path(taxonkit_dir)

    output:
        tuple val(meta), path(accession_taxid), path("*_taxonomy_table.txt"), emit: taxo_table

    script:
    """
        # Use TaxonKit to obtain taxonomy from taxIDs
        taxonkit reformat2 ${accession_taxid} -I 1 \
        --data-dir ${taxonkit_dir} \
        -f "{domain|acellular root|superkingdom};{phylum};{class};{order};{family};{genus};{species};{subspecies|strain|no rank}" \
        -r "NA" > ${meta.id}_taxonomy_table.txt
    """
}