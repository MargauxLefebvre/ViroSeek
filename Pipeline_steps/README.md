Steps of the pipeline
================
Margaux Lefebvre and Audric Berger
2025-05-12

# Set up directories

``` bash
# Set the variables
SAMPLEID=
FASTQDIR=
SPADESDIR=
DBDIAMOND=
PROTACCESION=
SILVAREF=
TAXONKITDIR=

TEMPDIR=./temp_$SAMPLEID
RESULTDIR=./taxo_$SAMPLEID

# Create the TEMPDIR and RESULTDIR
mkdir -p $TEMPDIR
mkdir -p $RESULTDIR
```

# Fastq filtering

## FastQC

Check the quality of the fastq.

``` bash
echo "--> FastQC for $SampleID before filtering"
mkdir $TEMPDIR/fastQC
zcat $FASTQDIR/${SampleID}_R1.fastq.gz $FASTQDIR/${SampleID}_R2.fastq.gz | fastqc stdin:${SampleID}.pretrim --outdir=$TEMPDIR/fastQC
```

## Trimgalore

Trim the adaptaters and low quality bases.

``` bash
echo "--> Trimgalore for $SampleID"
mkdir $TEMPDIR/trimgalore_fastq

trim_galore --paired --cores 8 --gzip \
--output_dir $TEMPDIR/trimgalore_fastq \
$FASTQDIR/${SampleID}_R1.fastq.gz \
$FASTQDIR/${SampleID}_R2.fastq.gz

# FastQC to check trimming effect
zcat $FASTQDIR/${SampleID}_R1.fastq.gz $FASTQDIR/${SampleID}_R2.fastq.gz | fastqc stdin:${SampleID}.trim --outdir=$TEMPDIR/fastQC
```

## BBduk and BBnorm

Remove ribosomal DNA

``` bash
echo "--> BBduk for $SampleID"
mkdir $TEMPDIR/BBduk

bbduk.sh \
    -Xmx34g \
    in1=$TEMPDIR/trimgalore_fastq/${SampleID}_R1.fastq.gz \
    in2=$TEMPDIR/trimgalore_fastq/${SampleID}_R2.fastq.gz  \
    out1=$TEMPDIR/BBduk/${SampleID}_R1.fastq.gz out2=$TEMPDIR/BBduk/${SampleID}_R1.fastq.gz \
    threads=6 \
    ref=$SILVAREF

# FastQC to check BBduk effect
zcat $FASTQDIR/${SampleID}_R1.fastq.gz $FASTQDIR/${SampleID}_R2.fastq.gz | fastqc stdin:${SampleID}.BBduk --outdir=$TEMPDIR/fastQC
```

Normalization of the reads count.

``` bash
echo "--> BBnorm for $SampleID"
mkdir $TEMPDIR/BBnorm

# Read 1
bbnorm.sh \
    in=$TEMPDIR/BBduk/${SampleID}_R1.fastq.gz \
    out=$TEMPDIR/BBnorm/${SampleID}_R1.fastq.gz \
    target=100 min=5 \
    threads=6 \
    -Xmx34g 

# Read 2
bbnorm.sh \
    in=$TEMPDIR/BBduk/${SampleID}_R2.fastq.gz \
    out=$TEMPDIR/BBnorm/${SampleID}_R2.fastq.gz \
    target=100 min=5 \
    threads=6 \
    -Xmx34g 

# FastQC to check BBnorm effect
zcat $FASTQDIR/${SampleID}_R1.fastq.gz $FASTQDIR/${SampleID}_R2.fastq.gz | fastqc stdin:${SampleID}.BBnorm --outdir=$TEMPDIR/fastQC
```

*Digital normalization is applied solely for the assembly step, while
the complete set of sequences is retained for quantification.*

# Create the assembly

``` bash
echo "--> Spades for $SampleID"

mkdir $TEMPDIR/assembly_spades
spades.py --rnaviral \
    -1 $TEMPDIR/BBnorm/${SampleID}_R1.fastq.gz \
    -2 $TEMPDIR/BBnorm/${SampleID}_R2.fastq.gz \
    --threads 12 \
    --memory 72 \
    -o $TEMPDIR/assembly_spades
```

# Quantification by mapping cleaned reads to the assembly

``` bash
echo "--> Minimap for $SampleID"
mkdir $TEMPDIR/Bam

minimap2 -t 16 -a $TEMPDIR/assembly_spades/spades.contigs.fa.gz \
        $TEMPDIR/BBduk/${SampleID}_R1.fastq.gz \
        $TEMPDIR/BBduk/${SampleID}_R1.fastq.gz \
        | samtools view -@ 16 -S -b - > $TEMPDIR/Bam/${SampleID}_map.bam

echo "--> Quantify reads count for $SampleID"

samtools sort -@ 16 $TEMPDIR/Bam/${SampleID}_map.bam -o $TEMPDIR/Bam/${SampleID}_sorted.bam

samtools markdup -@ 16 -r $TEMPDIR/Bam/${SampleID}_sorted.bam $TEMPDIR/Bam/${SampleID}_dedup.bam

samtools index $TEMPDIR/Bam/${SampleID}_dedup.bam

samtools idxstats $TEMPDIR/Bam/${SampleID}_dedup.bam > $TEMPDIR/Bam/${SampleID}_contigs_reads.tsv
```

# Taxonomic assignation

``` bash
echo "--> Diamond for $SampleID"
mkdir $TEMPDIR/Taxo

diamond blastx -p 16 -d $DBDIAMOND -q $TEMPDIR/assembly_spades/spades.contigs.fa.gz \
        -o $TEMPDIR/Taxo/${SampleID}.tsv --max-target-seqs 1

# Retrieve the accession IDs column from the Diamond file
    cut -f2 $TEMPDIR/Taxo/${SampleID}.tsv > $TEMPDIR/Taxo/${SampleID}_accession.txt

# Associate taxIDs with accession IDs
    grep -F -f $TEMPDIR/Taxo/${SampleID}_accession.txt $PROTACCESION > $TEMPDIR/Taxo/${SampleID}.accession_taxid.txt

# Extract taxIDs
    cut -f2 $TEMPDIR/Taxo/${SampleID}.accession_taxid.txt > $TEMPDIR/Taxo/${SampleID}.taxid.txt

# Use TaxonKit to obtain taxonomy from taxIDs
echo "--> TaxonKit for $SampleID"
    taxonkit lineage $TEMPDIR/Taxo/${SampleID}.taxid.txt --data-dir $TAXONKITDIR > $TEMPDIR/Taxo/${SampleID}_taxonomy_table.txt
```

# Clean the taxonomy table

``` bash
echo "--> Clean the taxonomy table for $SampleID"

# Merge contigs_reads.tsv and X.tsv files by contig number
    join -1 1 -2 1 <(sort "$TEMPDIR/Bam/${SampleID}_contigs_reads.tsv") <(sort "$TEMPDIR/Taxo/${SampleID}.tsv") > "$TEMPDIR/Taxo/${SampleID}_contig_reads_accession.txt"

# Add taxIDs
    join -1 5 -2 1 <(sort -k5,5 "$TEMPDIR/Taxo/${SampleID}_contig_reads_accession.txt") <(sort -k1,1 "$TEMPDIR/Taxo/${SampleID}.accession_taxid.txt") > "$TEMPDIR/Taxo/${SampleID}_contig_reads_accession_taxid.txt"

# Merge with taxonomy table
    join -1 16 -2 1 <(sort -k16,16 "$TEMPDIR/Taxo/${SampleID}_contig_reads_accession_taxid.txt") <(sort -k1,1 "$TEMPDIR/Taxo/${SampleID}_taxonomy_table.txt") > "$TEMPDIR/Taxo/${SampleID}_final_assembly.txt"

# Filter to keep only lines containing “virus”.
    grep -i "virus" "$TEMPDIR/Taxo/${SampleID}_final_assembly.txt" > "$TEMPDIR/Taxo/${SampleID}_filter_viral_taxonomy.txt"
```
