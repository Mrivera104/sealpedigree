# Batch 1 

I will be performing a genome alignment using the batch 1 resequenced eseal individuals and creating bam (binary alignment map) files as outputs. However, there is one sample that has very few reads (if any at all), and I will be resequencing that one when I start on batch 2 of resequencing. Just wanted to make a note of it :P 

For this alignment, I am using the BWA MEM aligner to align each forward and reverse fastq file to the reference genome. Here, I am using hap1 (primary) of the CCGP reference genome. 

# STEP 1: Alignment and BAM files

The first step in BWA MEM alignment is to index the reference fasta. This alignment process outputs a sam (sequence alignment/map format) file first, then converts that into a bam file. SAM files are typically converted to BAM (binary) files because they save storage space and are easier to manipulate. 

The BAM file then needs to be sorted. This step uses information from the reference genome to add chromosome/scaffold information. After sorting, the BAM file is indexed for easy lookup and organization. 

The bash script is named and run as the following: 

    bash align_fastq.sh

I went ahead and ran the bam files through MultiQC to do a quality check on them (https://seqera.io/multiqc/). My results are here: file:///C:/Users/Millie/Downloads/multiqc_report.html

I also had ChatGPT create a bash script for me to do this manually. 

    bash bam_check.sh





