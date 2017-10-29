valgrind --tool=memcheck --leak-check=full ./../randomacess sample.bam 1:10000-10010 | awk '{print $1,$2,$4,$10}' > out.txt

