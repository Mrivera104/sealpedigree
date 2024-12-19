#!/bin/bash

# Define input and output directories
INPUT_DIR="/home/migriver/eseal_raw_fastq"
OUTPUT_DIR="fastp_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all R1 files in the input directory
for R1 in "$INPUT_DIR"/*_R1_001.fastq.gz; do
    # Derive the R2 filename from the R1 filename
    R2="${R1/_R1.fastq.gz/_R2_001.fastq.gz}"

    # Extract the base name for output files
    BASENAME=$(basename "$R1" _R1_001.fastq.gz)

    # Define output filenames
    OUTPUT_R1="$OUTPUT_DIR/${BASENAME}_R1_clean.fastq.gz"
    OUTPUT_R2="$OUTPUT_DIR/${BASENAME}_R2_clean.fastq.gz"
    HTML_REPORT="$OUTPUT_DIR/${BASENAME}_fastp.html"
    JSON_REPORT="$OUTPUT_DIR/${BASENAME}_fastp.json"

    # Run fastp
    fastp \
        -i "$R1" -I "$R2" \
        -o "$OUTPUT_R1" -O "$OUTPUT_R2" \
        --html "$HTML_REPORT" \
        --json "$JSON_REPORT" \
        --thread 4

    echo "Processed: $BASENAME"
done

echo "All files processed. Output saved to $OUTPUT_DIR."
