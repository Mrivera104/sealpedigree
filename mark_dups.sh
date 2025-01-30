#!/bin/bash

# Define input and output directories
INPUT_DIR=/public/groups/meyerlab/eseal/batch1_bams/batch1_alignment # Replace with the path to your BAM files
OUTPUT_DIR="${INPUT_DIR}/duplicate_marked"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to process a single BAM file
process_bam() {
  BAM_FILE="$1"
  OUTPUT_DIR="$2"
  BASENAME=$(basename "$BAM_FILE" .bam)
  OUTPUT_BAM="${OUTPUT_DIR}/${BASENAME}_dupMarked.bam"
  METRICS_FILE="${OUTPUT_DIR}/${BASENAME}_dupMetrics.txt"
  
  # Run GATK MarkDuplicates
  echo "Processing $BAM_FILE..."
  gatk MarkDuplicates \
    -I "$BAM_FILE" \
    -O "$OUTPUT_BAM" \
    -M "$METRICS_FILE"
  
  if [ $? -eq 0 ]; then
    echo "Finished processing $BAM_FILE. Output saved to $OUTPUT_BAM."
  else
    echo "Error processing $BAM_FILE."
  fi
}

export -f process_bam

# Find all BAM files and process them in parallel
find "$INPUT_DIR" -name "*.bam" | parallel -j 64 process_bam {} "$OUTPUT_DIR"

echo "All BAM files have been processed. Duplicate-marked BAM files are in $OUTPUT_DIR."
