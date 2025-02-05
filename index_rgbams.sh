#!/bin/bash

# Define input and output directories
INPUT_DIR=/public/groups/meyerlab/eseal/batch1_bams/batch1_alignment/readgroup_bams # Replace with the path to your BAM files

# Function to process a single BAM file
process_bam() {
  BAM_FILE="$1"
  BASENAME=$(basename "$BAM_FILE" .bam)
  
  # Run samtools Index
  echo "Processing $BAM_FILE..."
  samtools index \
    "$BAM_FILE" 
  
  if [ $? -eq 0 ]; then
    echo "Finished processing $BAM_FILE."
  else
    echo "Error processing $BAM_FILE."
  fi
}

export -f process_bam

# Find all BAM files and process them in parallel
find "$INPUT_DIR" -name "*.bam" | parallel -j 64 process_bam {} "$OUTPUT_DIR"

echo "All BAM files have been processed."
