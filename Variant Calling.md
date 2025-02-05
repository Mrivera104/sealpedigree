# Variant Calling

The next step is to mark duplicates in each BAM file using GATK. This is what the SMSC class has to say about the duplicate marking step: 

"Potential PCR duplicates need to be marked. Marking duplicates makes sense even if you used a PCR-free library preparation procedure because reads identified as duplicates are not removed and can be included in the subsequent analyses if needed." 

    bash mark_dups.sh

After marking duplicates, we have to add read groups to the dup marked bam files. From the GATK website: 

"There is no formal definition of what a 'read group' is, however in practice this term refers to a set of reads that are generated from a single run of a sequencing instrument. "

"Many tools (Picard and GATK for example) require or assume the presence of at least one RG tag, defining a "read-group" to which each read can be assigned (as specified in the RG tag in the SAM record)."

    bash add_rg.sh

Next, we need to index the duplicate marked + rg-added bam files. First, we have to make a sequence dictionary for our fasta file that we can use to then index our bam files. This is an explanation from the GATK website: 

"Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no SAMRecords, and the header contains only sequence records. "

"Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records."

    gatk CreateSequenceDictionary -R=20230202.mMirAng1.NCBI.hap1.fasta -O=20230202.mMirAng1.NCBI.hap1.dict

We have to have the fasta file and the dict file in the same directory and call from that directory when indexing our bam file. 

    bash index_rgbams.sh

Now, we can finally move on to calling variants using gatk HaplotypeCaller. This is from the GATK website: 

"The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper." 

https://gatk.broadinstitute.org/hc/en-us/articles/9570334998171-HaplotypeCaller

    bash run_hapcaller.sh

