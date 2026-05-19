/*
 * Process: create the SILVA reference database used for filtering
 */
process prepare_silva_DB {
    label 'wget'
    
    publishDir "${params.outdir}/Database", mode: 'copy'

    output:
        path("SILVA.fa.gz"), emit: fasta

    script:
        """
        # Download the references database from https://www.arb-silva.de
        wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz
        wget -c https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
        
        # Merge the database together (>650Mo)
        zcat SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz | gzip > SILVA.fa.gz
        
        # Remove the downloaded fasta to save space 
        rm SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz 
        rm SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz 
        """
}

/*
 * Process: download the TaxonKit database
 */
process dwnload_taxonkit_DB {
    label 'wget'
    
    publishDir "${params.outdir}/Database", mode: 'copy'
        
    output:
      path("taxonkit"), emit: taxonkit_dir
      
    script:
        """
      # Download and untar the taxonomy dump (essential for linking protein sequences to taxonomic information with Taxonkit).
      mkdir -p taxonkit
      wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
      tar -xzf taxdump.tar.gz -C taxonkit
        """
}

/*
 * Process: Get viral taxids
 */
process prepare_viral_accessions {
    label 'taxonkit'
    
    input:
      path taxonkit_dir
        
    output:
      path("viral_taxids.txt"), emit: viral_taxids
      
    script:
        """
      # Keep all viral taxids (10239 is the id of viruses, at the root)
      taxonkit list --ids 10239 --data-dir ${taxonkit_dir} | sed 's/^[[:space:]]*//' > viral_taxids.txt
        """
}


/*
 * Process: download and subset accession2taxid file
 */
process prepare_accession2taxid {
    label 'wget'
    
    input:
      path viral_taxids    
    
    output:
      path("viral_accessions.txt"), emit: viral_accessions
      path("prot.accession2taxid.FULL.gz"), emit: accession2taxidFULL
      
    script:
        """
      # Download the accession-to-taxonomy mapping (this file maps protein accession numbers to NCBI Taxonomy IDs) >150Go uncompressed. Uncompress and keep compressed one too for DIAMOND
      wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
      
      # Extract the viral accessions
      zcat prot.accession2taxid.FULL.gz \
      | awk -F '\t' 'NR==FNR {tax[\$1]; next} FNR>1 && (\$2 in tax) {print \$1}' \
      ${viral_taxids} - > viral_accessions.txt
        """
}

/*
 * Process: download the DIAMOND databases
 */
process dwnload_diamond_DB {
    label 'wget'
        
    output:
      path("nr.gz"), emit: nr_database
      
    script:
        """
      # Fetch the non-redundant (NR) protein database in FASTA format from NCBI. >186Go
      wget -c https://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
        """
}

/*
 * Process: Subset only the viral sequence from the dataset
 */
process subset_diamond_DB {
    label 'seqkit'
    
    input:
        path viral_accessions
        path nr_database

    output:
      path("viral_nr.faa.gz"), emit: viral_ref

    script:
        """
      # Get only the viral sequences
      seqkit grep -f ${viral_accessions} ${nr_database} -j ${task.cpus} | gzip > viral_nr.faa.gz
        """
}

/*
 * Process: Subset only the viral sequence from the dataset
 */
process prepare_diamond_DB {
    label 'diamond'
    
    publishDir "${params.outdir}/Database", mode: 'copy'

    input:
      path accession2taxidFULL
      path viral_ref
      path taxonkit_dir

    output:
      path("ncbi-nr.viral.taxonomy.dmnd"), emit: diamond_db

    script:
        """
      # ----- Create a DIAMOND database with taxonomy mapping ~ 3h with 16 CPU -----
      # use diamond >= v2.1.12 to avoid new NCBI taxonomic ranks "cellular root", "acellular root", "domain" and "realm".
      
      gunzip -c ${viral_ref} | sed '/^>/s/ .*//' | diamond makedb --threads ${task.cpus} \
        --taxonmap ${accession2taxidFULL} \
        --taxonnames ${taxonkit_dir}/names.dmp \
        --taxonnodes ${taxonkit_dir}/nodes.dmp \
        --db ncbi-nr.viral.taxonomy.dmnd #name of the output DIAMOND database
        """
}