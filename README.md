# simple_bam_parser

**A simple BAM file parser**

BAM is the binary version of a SAM file. BAM is space efficient and allows fast random access, but not easy to parse like a ASCII based SAM file.
Luckily libraries such as htslib exists. Yet, the documentation is minimal that for a beginner it takes time to get familiar with.
Here are few simple programs that explains how to use those functions in htslib

At the moment there are 3 example programs : 

1. randomaccess : A simple program that uses htslib to perform random access to a BAM file.
- Opens a sorted BAM file and its index,  then
- Extracts and prints the reads aligned to the genomic region  provided as a command line argument

2. sequentialaccess : A simple program that performs sequential access of a bam file
- Opens a BAM file
- Print the reads to the standard output

3. write : A simple program that make a copy of a BAM file
- Opens a BAM file 
- Read read abd read and write to another BAM file


The source files are  

randomacess.c           - source code for 1
sequentialaccess.c      - source code for 2
write.c                 - source code for 3

common.h and common.c   : common functions used in example programs


##Running

First clone the repository and build htslib
```
got clone https://github.com/hasindu2008/simple_bam_parser.git
cd simple_bam_parser/htslib && make clean && make
```

Then build the example programs
```
cd .. && make
```

Then run a program that you wish,

./randomacess input.bam chr:start_pos-stop_pos 
./sequentialaccess input.bam
./write input.bam output.bam

Some examples : 
```
./randomacess test/sample.bam 1:10000-10005
./sequentialaccess test/sample.bam
./write test/sample.bam out.bam
```



