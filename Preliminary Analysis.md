# Starting notes
I began this thread on preliminary analysis on October 30th, 2024. I am going to be looking at the first 20 northern elephant seals I've resequenced (10 mother-offspring pairs, whole genome sequencing). This preliminary look will help me determine lots of things, including how many more individuals I should aim to sequence and how much coverage I'll need per individual. 

I will be putting code here that will help me streamline the QC of my samples and prepare them for analyses. This will be a slow process, and hopefully I will be learning a lot. I'm going to try a lot of this on my own, but will reach out for help when needed. I hope I can create a robust bioinformatic pipeline that will serve me in the future. LET'S GOOOOO

Putting this as a to-do list for myself, but: 
1.) I need to detail the methods for how I created the genome libraries for this initial batch of resequencing and how they were sequenced. 

Here is a neat little figure from the CD Genomics website (https://bioinfo.cd-genomics.com/whole-genome-resequencing-analysis.html): 
![Whole-Genome-Resequencing-Analysis-picture-3](https://github.com/user-attachments/assets/2474b351-d0c0-4924-8453-c702489a23af)

This should thus be the data analysis workflow I am aiming to replicate. :) 

# Step 1: Quality Checking Fastq Files using FastQC

I used fastQC to check my raw fastq files. It outputs stats like adapter content and per base sequence content. 

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

Here is an example of the basic stats reported by fastqc: 

![image](https://github.com/user-attachments/assets/39e34ac0-de8e-4d0e-ae71-29c3513c3413)

There were one or two sequences that are a little iffy, but for the most part, they got the greenlight aside from adapter content which is expected. Now it's time to trim and filter the raw fastq files to prepare for alignment. 

# Step 2: Quality Control

I used fastp to go through and clean up/trim my fastq files. I used the default parameters and haave HTML/JSON reports for all. 

    # Run fastp
    fastp \
        -i "$R1" -I "$R2" \
        -o "$OUTPUT_R1" -O "$OUTPUT_R2" \
        --html "$HTML_REPORT" \
        --json "$JSON_REPORT" \
        --thread 4
        
Here is an example of the results: 

![Screenshot 2024-12-19 140928](https://github.com/user-attachments/assets/fc91d582-d866-4bdf-b2b7-37c9441430fb)


