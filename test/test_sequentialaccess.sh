samtools view sample.bam | awk '{print $1,$2,$3,$4,$10}' > expected.txt
./../sequentialaccess sample.bam | grep -v "Chromosome ID" | awk '{print $1,$2,$3,$4,$10}' > out.txt
diff expected.txt out.txt

