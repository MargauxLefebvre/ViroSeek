Steps of the pipeline
================
Margaux Lefebvre and Audric Berger
2025-07-05

# Set up directories and environment

``` bash
# Activate the environment
conda activate ViroSeek

# Set the variables
SampleID=sample_ID
FASTQDIR=/path/to/fastq/directory
TEMPDIR=./temp_$SampleID
RESULTDIR=./results/taxo_$SampleID

SILVAREF=/path/to/SILVA/fasta
DBDIAMOND=/path/to/ncbi-nr.taxonomy.dmnd
PROTACCESION=/path/to/prot.accession2taxid.txt  # this file must be decompressed
TAXONKITDIR=/path/to/Taxonkit/files
SPADESDIR=/path/to/SPades/bin

LENGTH_SEQ='0'
DIAM_EVALUE='0.001'
DIAM_ID='0'
DIAM_QUERYCOV='0'
DIAM_SENSI=''

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
zcat $TEMPDIR/trimgalore_fastq/${SampleID}_R1_val_1.fq.gz $TEMPDIR/trimgalore_fastq/${SampleID}_R2_val_2.fq.gz | fastqc stdin:${SampleID}.trim --outdir=$TEMPDIR/fastQC
```

## BBduk

Remove ribosomal DNA

``` bash
echo "--> BBduk for $SampleID"
mkdir $TEMPDIR/BBduk

bbduk.sh \
    -Xmx34g \
    in1=$TEMPDIR/trimgalore_fastq/${SampleID}_R1_val_1.fq.gz \
    in2=$TEMPDIR/trimgalore_fastq/${SampleID}_R2_val_2.fq.gz  \
    out1=$TEMPDIR/BBduk/${SampleID}_R1_rmrdna.fastq.gz out2=$TEMPDIR/BBduk/${SampleID}_R2_rmrdna.fastq.gz \
    threads=6 \
    ref=$SILVAREF

# FastQC to check BBduk effect
zcat $TEMPDIR/BBduk/${SampleID}_R1_rmrdna.fastq.gz $TEMPDIR/BBduk/${SampleID}_R2_rmrdna.fastq.gz | fastqc stdin:${SampleID}.BBduk --outdir=$TEMPDIR/fastQC
```

# Create the assembly

``` bash
echo "--> Spades for $SampleID"

# Merge paired-end FASTQ files into a single interleaved FASTQ file.
seqtk mergepe $TEMPDIR/BBduk/${SampleID}_R1_rmrdna.fastq.gz $TEMPDIR/BBduk/${SampleID}_R2_rmrdna.fastq.gz | bgzip > $TEMPDIR/BBduk/${SampleID}_inter_rmrdna.fastq.gz

# Make the assembly
mkdir $TEMPDIR/assembly_spades
$SPADESDIR/spades.py --rnaviral \
    --12 $TEMPDIR/BBduk/${SampleID}_inter_rmrdna.fastq.gz \
    --threads 12 \
    --memory 72 \
    -o $TEMPDIR/assembly_spades

# Filter the assembly
bioawk -c fastx "{ if(length($seq) > $LENGTH_SEQ) { print '>'$name; print $seq }}" $TEMPDIR/assembly_spades/contigs.fasta > $TEMPDIR/assembly_spades/contigs.filtered.fasta
```

# Quantification by mapping cleaned reads to the assembly

``` bash
echo "--> Minimap for $SampleID"
mkdir $TEMPDIR/Bam

minimap2 -t 16 -a $TEMPDIR/assembly_spades/contigs.filtered.fasta \
        $TEMPDIR/BBduk/${SampleID}_R1_rmrdna.fastq.gz \
        $TEMPDIR/BBduk/${SampleID}_R2_rmrdna.fastq.gz \
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

diamond blastx -p 16 -F 15 -d $DBDIAMOND -q $TEMPDIR/assembly_spades/contigs.filtered.fasta \
        -o $TEMPDIR/Taxo/${SampleID}.tsv --max-target-seqs 1 -e $DIAM_EVALUE --id $DIAM_ID --query-cover $DIAM_QUERYCOV \
        --very-sensitive --range-culling $DIAM_SENSI

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
   grep -i "virus" "$TEMPDIR/Taxo/${SampleID}_final_assembly.txt" > "$TEMPDIR/Taxo/${SampleID}_virus_taxonomy.txt"
   uniq "$TEMPDIR/Taxo/${SampleID}_virus_taxonomy.txt" > "$TEMPDIR/Taxo/${SampleID}_filter_viral_taxonomy.txt"

# Keep only reads count and the taxonomy in a tab delimited file
  awk '{printf "%s\t", $5; for (i=17; i<=NF; i++) printf "%s%s", $i, (i<NF?" ":"\n")}' "$TEMPDIR/Taxo/${SampleID}_filter_viral_taxonomy.txt" > "$TEMPDIR/Taxo/${SampleID}_taxonomy_viral.txt"
 sed -i 's/;/\t/g' "$TEMPDIR/Taxo/${SampleID}_taxonomy_viral.txt"
```

# Finalize the pipeline and clean the intermediate files

``` bash
# Get the final files
cp -r $TEMPDIR/Bam/${SampleID}_contigs_reads.tsv $RESULTDIR
cp -r $TEMPDIR/fastQC/ $RESULTDIR
cp -r $TEMPDIR/Taxo/${SampleID}_taxonomy_viral.txt $RESULTDIR/${SampleID}_taxonomy_viral.clean.txt
cp -r $TEMPDIR/assembly_spades $RESULTDIR
cp -r $TEMPDIR/Taxo/${SampleID}_filter_viral_taxonomy.txt $RESULTDIR/${SampleID}_all_viral_taxonomy.txt
cp -r $TEMPDIR/Taxo/${SampleID}.tsv $RESULTDIR
# Delete intermediate files
rm -rf $TEMPDIR
```
