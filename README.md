ViroSeek
================
Margaux Lefebvre and Audric Berger
2025-06-11

# Install the software for the pipeline

Install all the tools, with the [yml
file](https://github.com/MargauxLefebvre/ViroSeek/blob/main/ViroSeek_env.yml)…
(easiest and adviced option)

``` bash
conda env create -f ViroSeek_env.yml
```

… or each tool in conda

``` bash
conda create -n ViroSeek python=3.8
conda activate ViroSeek
conda install minimap2=2.29 -c bioconda
conda install samtools=1.21 -c bioconda
conda install taxonkit=0.9.0 -c bioconda
conda install diamond=2.1.10 -c bioconda
conda install trim-galore=0.6.10 -c bioconda
conda install bbmap=39.18 -c bioconda
conda install seqtk=1.4 -c bioconda
```

*Note: python v3.8 is the minimum version for Spades.*

Most of the tools can be installed with conda but we need the binaries
for [SPAdes
v4.0.0](https://github.com/ablab/spades/releases/tag/v4.0.0).

# Versions and manual for each software

- [FastQC
  v0.12.1](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trim-galore
  v0.6.10](https://github.com/FelixKrueger/TrimGalore?tab=readme-ov-file)
- [BBmap v39.18 (for BBduk and
  BBnorm)](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/)
- [SPAdes v4.0.0](https://github.com/ablab/spades)
- [Minimap2 v2.29](https://github.com/lh3/minimap2)
- [Diamond v2.1.10](https://github.com/bbuchfink/diamond)
- [Taxonkit v0.9.0](https://bioinf.shenwei.me/taxonkit/)
- [Samtools v1.21](http://www.htslib.org/doc/1.21/samtools.html)
- [Seqtk v1.4](https://github.com/lh3/seqtk)

# Generate the references databases

## SILVA

Reads corresponding to ribosomal DNA sequences are filtered out using
reference data from the [SILVA database](https://www.arb-silva.de/).
This pipeline was developed and tested using the latest available
release at the time (version 138.2).

Below, we explain how to create the SILVA reference database used for
this filtering.

``` bash
# Download the references database
wget -c https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/ARB_files/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb.gz
wget -c https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/ARB_files/SILVA_138.2_LSURef_NR99_03_07_24_opt.arb.gz
# Merge the database together
zcat SILVA_138.2*.gz > SILVA.fasta
```

## Diamond and Taxonkit

Taxonomic assignment requires a reference database built with DIAMOND
along with corresponding taxonomy mapping files.

Below, we provide instructions for generating the DIAMOND database and
preparing the necessary files used by TaxonKit to enable taxonomy
assignment.

``` bash
conda activate ViroSeek # to have the same Diamond version than the one used in the pipeline

# Fetch the non-redundant (NR) protein database in FASTA format from NCBI.
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
 
# Download and untar the taxonomy dump (essential for linking protein sequences to taxonomic information with Taxonkit).
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xfz taxdump.tar.gz
 
# Download the accession-to-taxonomy mapping (this file maps protein accession numbers to NCBI Taxonomy IDs).
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
 
# Create a DIAMOND database with taxonomy mapping
gunzip -c nr.gz | sed '/^>/s/ .*//' | diamond makedb --threads 16 \
  --taxonmap prot.accession2taxid.FULL.gz \
  --taxonnames names.dmp \
  --taxonnodes nodes.dmp \
  --db ncbi-nr.taxonomy.dmnd #name of the output DIAMOND database
  
# Make the directory for Taxonkit and transfert the accession-to-taxonomy mapping files
mkdir $TAXONKITDIR
mv *.dmp $TAXONKITDIR

# Decompress the prot.accession2taxid.FULL.gz file
gzip -d prot.accession2taxid.FULL.gz
```

**Note:**

When building the DIAMOND database in June 2025, we encountered a
compatibility issue between the NCBI taxonomy files and DIAMOND (see
[issue \#352](https://github.com/bbuchfink/diamond/issues/352)).
Specifically, the NCBI taxonomy (nodes.dmp) included taxonomic ranks
such as domain and realm, which were not recognized by DIAMOND. To
resolve this, we replaced these terms before creating the database:
`sed -i 's/domain/superkingdom/g' nodes.dmp` and
`sed -i 's/realm/kingdom/g' nodes.dmp`. These replacements ensure
compatibility with DIAMOND expected taxonomy rank names during database
creation.

# How to run it

*TO DO*
