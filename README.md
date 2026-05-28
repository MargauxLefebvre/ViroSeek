ViroSeek
================
<img align="right" src="docs/img/IRD.png" width="225" height="66" /> <img align="right" src="docs/img/MIVEGEC.png" width="100" height="66" />

ViroSeek is a pipeline designed for the analysis of viromes derived from target-enriched libraries, with a particular emphasis on eukaryotic viruses, including arboviruses of public health concern.

# Table of Contents

   * [Installation](#installation)
      * [Nextflow](#nextflow)
      * [Container platform](#container-platform)
        * [Docker](#docker)
        * [Singularity](#singularity)
      * [Softwares](#softwares)
      * [Databases](#databases)
        * [Diamond and Taxonkit](#diamond-and-taxonkit)
  * [Usage](#usage)
    * [Parameters](#parameters)
      * [Mandatory parameters](#mandatory-parameters)
      * [Optional parameters](#optional-parameters)
    * [Example](#example)
    * [Profiles](#profiles)

# Installation

The prerequisites to run the pipeline are:  

  * [Nextflow](https://www.nextflow.io/)  >= 22.04.0
  * [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/)  

## Nextflow 

  * Via conda 

    <details>
      <summary>See here</summary>
      
      ```bash
      conda create -n nextflow
      conda activate nextflow
      conda install bioconda::nextflow
      ```  
    </details>

  * Manually
    <details>
      <summary>See here</summary>
      Nextflow runs on most POSIX systems (Linux, macOS, etc) and can typically be installed by running these commands:

      ```bash
      # Make sure 11 or later is installed on your computer by using the command:
      java -version
      
      # Install Nextflow by entering this command in your terminal(it creates a file nextflow in the current dir):
      curl -s https://get.nextflow.io | bash 
      
      # Add Nextflow binary to your user's PATH:
      mv nextflow ~/bin/
      # OR system-wide installation:
      # sudo mv nextflow /usr/local/bin
      ```
    </details>

## Container platform

To run the workflow you will need a container platform: docker or singularity.

### Docker

Please follow the instructions at the [Docker website](https://docs.docker.com/desktop/)

### Singularity

Please follow the instructions at the [Singularity website](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)

## Softwares

Software is provided as containers, which help ensure reproducibility and are automatically retrieved by the pipeline.  
To access information about software versions, please refer to `config/softwares.config`.  
You can easily update the version of any software by modifying the corresponding software image in that file.

## Databases

If you don't want the default ViroSeek database, you can prepare your reference database.
    
### Diamond and Taxonkit

By default, ViroSeek uses a database derived from the NCBI reference sequence non-redundant protein database, restricted to viral sequences. However, users may specify an alternative database. 

Below, we provide instructions for generating the alternative DIAMOND database and preparing the necessary files used by TaxonKit to enable taxonomy assignment.

  <details>
    <summary>See here</summary>
      
      ```bash
    # Download and untar the taxonomy dump (essential for linking protein sequences to taxonomic information with Taxonkit).
    mkdir -p taxonkit
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz -C taxonkit
    
    # Download the accession-to-taxonomy mapping (this file maps protein accession numbers to NCBI Taxonomy IDs) >150Go uncompressed. Uncompress and keep compressed one too for DIAMOND
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
    
    # ----- Create a DIAMOND database with taxonomy mapping ~ 3h with 32 CPU - ~365 Go -----
    # --taxonmap accepts compressed and uncompressed input
    # /!\ Be careful to use same DIAMOND version than the one used by ViroSeek (see config/softwares.config file).
    
    # First fetch the non-redundant (NR) protein database in FASTA format from NCBI. >186Go
    wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
    # Now make the diamon DB
    gunzip -c nr.gz | sed '/^>/s/ .*//' | diamond makedb --threads 32 \
      --taxonmap prot.accession2taxid.FULL.gz \
      --taxonnames taxonkit/names.dmp \
      --taxonnodes taxonkit/nodes.dmp \
      --db ncbi-nr.taxonomy.dmnd #name of the output DIAMOND database
      ```  
  </details>
    
# Usage

## Parameters

### Mandatory parameters

| Parameter | Description |
|----|----|
| `--input` | Path to a **CSV file** with 3 columns (with a header): `sample,read1,read2.`. Each row represents a paired-end sample. |

Here an example of an input:

```
sample,read1,read2
SIM,/path/to/fastq/SIM_R1.fastq.gz,/path/to/fastq/SIM_R2.fastq.gz
AltMix,/path/to/fastq/AltMix_R1.fastq.gz,/path/to/fastq/AltMix_R2.fastq.gz
MixA,/path/to/fastq/MixA_R1.fastq.gz,/path/to/fastq/MixA_R2.fastq.gz
MixB,/path/to/fastq/MixB_R1.fastq.gz,/path/to/fastq/MixB_R2.fastq.gz
MixC,/path/to/fastq/MixC_R1.fastq.gz,/path/to/fastq/MixC_R2.fastq.gz
```

### Optional parameters

| Parameter | Default value | Description |
|----|----|----|
| `-work-dir` | `./work` | Directory for Nextflow’s temporary working files. |
| `--outdir` | `./results` | Directory where output results are saved. |
| `--trim` | (empty) | Tools to use to trim: `trimgalore`, `fastp`. Empty is accepted and means no trimming. |
| `--trim_opt` | (empty) | Option to set to trimming tool chosen. |
| `--skip_dedup` | false | Skip the deduplication step with `samtools markdup`. |
| `--length_seq` | `0` | Minimum length (in bp) to keep contigs after SPAdes assembly. Use `0` to **keep all** contigs without filtering. |
| `--conta_ref` | Non viral ribosomal DNA (16S/18S and 23S/28S) sequences from SILVA rRNA database (release 138). | Path to a reference file containing sequences to exclude as potential contaminants (e.g. a host genome), provided in FASTA format (compressed or uncompressed). Alternatively, a comma-separated list of URLs may be supplied to automatically download and use the reference sequences. |
| `--diamond_db` | NCBI RefSeq non-redundant protein database restricted to viral sequences. | Path to the **DIAMOND-formatted database** for protein sequence taxonomic assignment. ([DIAMOND setup example](https://github.com/MargauxLefebvre/ViroSeek#diamond-and-taxonkit)) |
| `--taxonkit_dir` | NCBI taxonomy files. | Path to a directory containing the **NCBI taxonomy dump files** required by TaxonKit (`names.dmp`, `nodes.dmp`, etc.). ([DIAMOND setup example](https://github.com/MargauxLefebvre/ViroSeek#diamond-and-taxonkit)) |
| `--diam_evalue` | `0.001` | Maximum **e-value** for DIAMOND alignment hits. Lower values increase stringency. |
| `--diam_id` | `0` | Minimum **percent identity** required for DIAMOND hits (0–100). |
| `--diam_querycov` | `0` | Minimum **query coverage** percentage required for DIAMOND hits. |
| `--diam_sensi` | (empty) | The sensitivity modes of Diamond. The accepted values are: `--faster`, `--fast`, `--mid-sensitive`, `--sensitive`, `--more-sensitive`, `--very-sensitive` and `--ultra-sensitive`. For more details, see [DIAMOND wiki](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes). |

**Note:** it is also possible to specify which configuration file to use
with the `-c' option '/path/to/the/nextflow.config'`.

## Example

To get help run:

```bash
nextflow run MargauxLefebvre/ViroSeek -r v0.0.2 --help
```

A typical command might look like the following:

``` bash
nextflow run MargauxLefebvre/ViroSeek -r v0.0.2 --input '/path/to/samples.csv' -profile singularity,slurm \
  -work-dir '/your/work_dir/' \
  --outdir '/your/results/dir/' \
  --trim 'trimgalore' \
  --length_seq 140 --diam_evalue 1e-200 --diam_id 96.0 --diam_querycov 51.0 \
  --diam_sensi='--fast'
```

*Replace file paths and profiles as needed*

**Note:** 

The command `nextflow run MargauxLefebvre/ViroSeek -r v0.0.2` allows you to run the pipeline directly from the v0.0.2 release without any manual installation (apart from the required dependencies, i.e. Nextflow and a supported container platform; see the [Installation](#installation) section). Nextflow automatically downloads and manages the pipeline code locally.

Alternatively, you can clone the repository and run the pipeline from a local copy. In this case, the command becomes `nextflow run ViroSeek.nf`. Follow the steps below:

```
# clone the repository
git clone https://github.com/MargauxLefebvre/ViroSeek.git

# move into the repository
cd ViroSeek

# run the pipeline (example: display help)
nextflow run ViroSeek.nf --help
```

## Profiles 

### Profiles for the containers

- For Docker, use `-profile docker`
- For singularity, use `-profile singularity`

### Profiles for HPC environments

- For SLURM, use `-profile slurm`
- For SGE, use `-profile sge`
- For local use with more limited resources,  use `-profile local`

## What if a job stops and I want to restart it?

One of Nextflow key strengths is its ability to resume from where it
left off after a failure or interruption. By default your pipeline will restart
without re-running completed tasks. To deactive this behaviour you can set 
the `resume` parameter in the `nextflow.config` file to `false`.

If your job fails due to resource limitations (e.g., RAM or CPU), you
can adjust the resource settings in the `nextflow.config` file. Modify
the CPU or memory allocation for the affected process to better suit
your system before rerunning.

## Contributing

Contributions from the community are welcome ! See the [Contributing guidelines](https://github.com/MargauxLefebvre/ViroSeek/blob/main/CONTRIBUTING.md).

## How to cite

If you use ViroSeek, please cite: **ViroSeek: a viral detection pipeline for second-generation sequencing** by Audric Berger, Margaux J. M. Lefebvre, Jacques Dainat, Davy Jiolle, Isabelle Conclois, Loic Talignani, Emilio Mastriani, Sylvie Cornelie, Nicolas Berthet, Christophe Paupy. [Preprint available](https://doi.org/10.64898/2026.03.04.706323)
