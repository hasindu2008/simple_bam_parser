./../write sample.bam out.bam 
samtools view -h sample.bam > expected.txt
samtools view -h out.bam > out.txt
diff expected.txt out.txt

