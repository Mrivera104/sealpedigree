#!/bin/bash

# Define input and output directories
INPUT_DIR="/public/groups/meyerlab/eseal/batch1_bams/batch1_alignment/readgroup_bams"  # Replace with the path to your BAM files
OUTPUT_DIR="/public/groups/meyerlab/eseal/batch1_bams/vcf_files"

# Function to process a single BAM file
process_bam() {
  BAM_FILE="$1"
  OUTPUT_DIR="$2"
  
  # Extract the file name and remove the unwanted suffix
  FILE=$(basename "$BAM_FILE")
  BASENAME=${FILE%_sorted_dupMarked.bam_RG.bam}
  
  OUTPUT_VCF="${OUTPUT_DIR}/${BASENAME}.g.vcf.gz"
  
  # Run GATK HaplotypeCaller
  echo "Processing $BAM_FILE..."
  gatk HaplotypeCaller \
    -I "$BAM_FILE" \
    -R /public/groups/meyerlab/eseal/batch1_bams/batch1_alignment/readgroup_bams/20230202.mMirAng1.NCBI.hap1.fasta \
    -ERC GVCF \
    -O "$OUTPUT_VCF" 
  
  if [ $? -eq 0 ]; then
    echo "Finished processing $BAM_FILE. Output saved to $OUTPUT_VCF."
  else
    echo "Error processing $BAM_FILE."
  fi
}

export -f process_bam

# Find all BAM files and process them in parallel using GNU parallel
find "$INPUT_DIR" -name "*.bam" | parallel -j 64 process_bam {} "$OUTPUT_DIR"

echo "All BAM files have been processed. GVCF files are in $OUTPUT_DIR."
