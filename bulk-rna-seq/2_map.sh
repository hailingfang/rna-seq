#Example, using bowtie2 to mapped pair end read to transcriptom.

for ((i=1; i<4; i++))
do
    echo $i
    bowtie2 -x ../../../refdata/rna-from-genomic-bowtie2/GCF_000002035.6_GRCz11_rna_from_genomic \
        -1 ../1-clean-data/dsR${i}-clean_L1_1.fq.gz \
        -2 ../1-clean-data/dsR${i}-clean_L1_2.fq.gz \
        -S ../2-mapped/dsR${i}-clean.sam -k 3 --threads 8

done

