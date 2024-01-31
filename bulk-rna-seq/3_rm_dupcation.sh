#Example, using samtools to remove duplication reads.

for ((i=1; i<4; i++))
do
    echo $i
    samtools fixmate --threads 8 -m -O bam ../2-mapped/dsR${i}-clean.sam ../2-mapped/dsR${i}-clean-fixmate.bam
    samtools sort --threads 8 -O bam -o ../2-mapped/dsR${i}-clean-fixmate-sorted.bam ../2-mapped/dsR${i}-clean-fixmate.bam
    samtools markdup -O bam --threads 8 -r ../2-mapped/dsR${i}-clean-fixmate-sorted.bam ../2-mapped/dsR${i}-clean-fixmate-sorted-rmdup.bam 
    samtools sort --threads 8 -n -O sam  -o ../2-mapped/dsR${i}-clean-fixmate-sorted-rmdup-sorted.sam ../2-mapped/dsR${i}-clean-fixmate-sorted-rmdup.bam 

done
