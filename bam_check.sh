#!/bin/bash

# Define the directory containing BAM files
BAM_DIR="/home/migriver/eseal_batch1/batch1_alignment"
OUTPUT_DIR="/home/migriver/eseal_batch1/bam_check"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all BAM files in the directory
for BAM in "$BAM_DIR"/*.bam; do
    # Extract the base name for output files
    BASENAME=$(basename "$BAM" .bam)
    RESULT_FILE="$OUTPUT_DIR/${BASENAME}_check_results.txt"

    echo "Checking $BAM..." > "$RESULT_FILE"

    # Check if BAM file is indexed
    if [ ! -f "${BAM}.bai" ]; then
        echo "Indexing $BAM..." >> "$RESULT_FILE"
        samtools index "$BAM" || { echo "Failed to index $BAM" >> "$RESULT_FILE"; continue; }
    fi

    # Check the BAM header
    echo "Header information:" >> "$RESULT_FILE"
    samtools view -H "$BAM" >> "$RESULT_FILE" || echo "Failed to read header for $BAM" >> "$RESULT_FILE"

    # Run flagstat to get alignment summary
    echo "Alignment statistics (flagstat):" >> "$RESULT_FILE"
    samtools flagstat "$BAM" >> "$RESULT_FILE" || echo "Failed to run flagstat for $BAM" >> "$RESULT_FILE"

    # Compute average coverage across the genome
    echo "Calculating average coverage..." >> "$RESULT_FILE"
    AVG_COVERAGE=$(samtools depth "$BAM" | awk '{sum+=$3; cnt++} END {if (cnt > 0) print sum/cnt; else print "No coverage data"}')
    echo "Average coverage across the genome: ${AVG_COVERAGE}x" >> "$RESULT_FILE"

    # Detailed alignment statistics
    echo "Detailed alignment statistics:" >> "$RESULT_FILE"
    samtools stats "$BAM" >> "$RESULT_FILE" || echo "Failed to compute detailed stats for $BAM" >> "$RESULT_FILE"

    # Optional: Generate plots (requires gnuplot)
    echo "Generating plots..." >> "$RESULT_FILE"
    samtools stats "$BAM" > "$OUTPUT_DIR/${BASENAME}.stats"
    plot-bamstats -p "$OUTPUT_DIR/${BASENAME}" "$OUTPUT_DIR/${BASENAME}.stats" >> "$RESULT_FILE" 2>/dev/null || echo "Failed to generate plots for $BAM" >> "$RESULT_FILE"

    echo "$BAM check completed. Results saved to $RESULT_FILE"
done

echo "All BAM files checked. Results saved to $OUTPUT_DIR."
