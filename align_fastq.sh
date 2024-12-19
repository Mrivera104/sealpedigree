#!/bin/bash

# Define input, output, and reference directories
INPUT_DIR="/home/migriver/eseal_raw_fastq/fastp_output"
OUTPUT_DIR="/home/migriver/eseal_batch1/batch1_alignment"
REFERENCE="/home/migriver/eseal_batch1/20230202.mMirAng1.NCBI.hap1.fasta"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Index the reference genome (if not already indexed)
if [ ! -f "$REFERENCE.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$REFERENCE" || { echo "Failed to index $REFERENCE"; exit 1; }
fi

# Loop through all cleaned R1 files in the input directory
for R1 in "$INPUT_DIR"/*_R1_clean.fastq.gz; do
    # Derive the R2 filename from the R1 filename
    R2="${R1/_R1_clean.fastq.gz/_R2_clean.fastq.gz}"

    # Check if R2 exists
    if [ ! -f "$R2" ]; then
        echo "Paired file $R2 not found for $R1. Skipping..."
        continue
    fi

    # Extract the base name for output files
    BASENAME=$(basename "$R1" _R1_clean.fastq.gz)

    # Define output filenames
    SAM_FILE="$OUTPUT_DIR/${BASENAME}.sam"
    SORTED_BAM_FILE="$OUTPUT_DIR/${BASENAME}_sorted.bam"

    # Align reads to the reference genome using BWA-MEM
    echo "Aligning $BASENAME..."
    bwa mem -t 80 "$REFERENCE" "$R1" "$R2" > "$SAM_FILE" || { echo "Alignment failed for $BASENAME"; exit 1; }

    # Convert SAM to BAM, sort, and remove the SAM file
    echo "Converting and sorting $BASENAME..."
    samtools view -@ 80 -b "$SAM_FILE" | samtools sort -@ 80 -o "$SORTED_BAM_FILE" || { echo "SAM to BAM conversion failed for $BASENAME"; exit 1; }
    rm "$SAM_FILE"

    # Index the sorted BAM file
    samtools index "$SORTED_BAM_FILE" || { echo "Indexing failed for $BASENAME"; exit 1; }

    echo "$BASENAME processing complete."
done

echo "All files aligned and processed. Output saved to $OUTPUT_DIR."

