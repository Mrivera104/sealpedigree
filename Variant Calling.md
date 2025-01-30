# Variant Calling

The next step is to mark duplicates in each BAM file using GATK. This is what the SMSC class has to say about the duplicate marking step: 

"Potential PCR duplicates need to be marked. Marking duplicates makes sense even if you used a PCR-free library preparation procedure because reads identified as duplicates are not removed and can be included in the subsequent analyses if needed." 
