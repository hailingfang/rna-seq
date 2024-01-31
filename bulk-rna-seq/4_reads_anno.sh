# Example, using annoread to annotate mapped reads (sam file is required).

for ((i=1; i<4; i++))
do
    echo $i
    python3 /home/fanghl/work/git_repo/annoread/annoread.py \
        --ref tran \
        --ref-type ncbi-genome-rna \
        --gtf /home/comm/data/reference-genome/Danio-rerio/GRCz11/GCF_000002035.6_GRCz11_genomic.gtf \
        --fasta /home/comm/data/reference-genome/Danio-rerio/GRCz11/GCF_000002035.6_GRCz11_rna_from_genomic.fna \
        --seq-type pe \
        --out ../3-reads-annotation/dsR${i}-rmdup-sorted.rda \
        ../2-mapped/dsR${i}-clean-fixmate-sorted-rmdup-sorted.sam

done

