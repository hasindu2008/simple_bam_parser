samtools view sample.bam 1:10000-10005  | awk '{print $1,$2,$3,$4,$10}' > expected.txt
./../randomacess sample.bam 1:10000-10005 | grep -v "Chromosome ID" | awk '{print $1,$2,$3,$4,$10}' > out.txt
diff expected.txt out.txt

