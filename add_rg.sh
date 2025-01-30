#!/bin/bash

# Set paths to GATK and input directory (adjust if necessary)
INPUT_DIR=/public/groups/meyerlab/eseal/batch1_bams/batch1_alignment/duplicate_marked  # Change to your BAM file directory
OUTPUT_DIR=/public/groups/meyerlab/eseal/batch1_bams/batch1_alignment/readgroup_bams  # Change if you want a different output folder

# Loop through all duplicate-marked BAM files
for BAM in "$INPUT_DIR"/*_dupMarked.bam; do
    # Extract filename without extension
    BASENAME=$(basename "$BAM" "_duplMarked.bam")
    
    # Define output file name
    OUTPUT_BAM="$OUTPUT_DIR/${BASENAME}_RG.bam"
    
    echo "Processing $BAM -> $OUTPUT_BAM"

    # Run GATK AddOrReplaceReadGroups
    gatk AddOrReplaceReadGroups \
        -I "$BAM" \
        -O "$OUTPUT_BAM" \
        -RGID "$BASENAME" \
        -RGLB "lib1" \
        -RGPL "ILLUMINA" \
        -RGPU "unit1" \
        -RGSM "$BASENAME"

    echo "Finished processing $BAM"
done

echo "All BAM files processed successfully."
