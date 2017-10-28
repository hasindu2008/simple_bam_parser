samtools view ~/store1/Platypus_human_WGS/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam 1:10000-10010 | awk '{print $1,$2,$4,$10}' > sams.txt
./randomacess ~/store1/Platypus_human_WGS/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam 1:10000-10010 | awk '{print $1,$2,$4,$10}' > ours.txt
diff sams.txt ours.txt

