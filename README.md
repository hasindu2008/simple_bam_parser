# simple_bam_parser

**A simple BAM file parser**


BAM is the binary version of a SAM file. BAM is space efficient and allows fast random access, but not easy to parse like a ASCII based SAM file.
Luckily libraries such as htslib exists. Yet, the documentation is minimal that for a beginner it takes time to get familiar with.
This is a simple program that uses htslib to perform random access to a BAM file. First it,
- opens a sorted BAM file and its index,  then
- Extracts and prints the reads aligned to the genomic region  provided as a command line argument
