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
        # Extract taxIDs
        cut -f2 ${accession_taxid} > ${meta.id}.taxid.txt

        # Use TaxonKit to obtain taxonomy from taxIDs
        taxonkit lineage ${meta.id}.taxid.txt --data-dir ${taxonkit_dir} > ${meta.id}_taxonomy_table.txt
    """
}