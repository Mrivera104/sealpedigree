#!/bin/bash

# Directory where the FASTQ files are located
FASTQ_DIR="/home/migriver/eseal_raw_fastq"
# Directory to save the FastQC results
OUTPUT_DIR="/home/migriver/eseal_raw_fastq/fastqc_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC on each FASTQ file
for file in "$FASTQ_DIR"/*_S*_R*_001.fastq.gz
do
    echo "Running FastQC on $file..."
    fastqc "$file" -o "$OUTPUT_DIR"
done

echo "FastQC analysis complete! Results are saved in $OUTPUT_DIR."
