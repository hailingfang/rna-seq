#An example using cutadapt to remove adapters.

for((i=1; i<4; i++))
do
    echo $i
    cutadapt -a AGATCGGAAGAGCACACG -A AGATCGGAAGAGCGTCGT \
        --minimum-length 50 \
        -q 15,10 \
        --trim-n \
        -o ../1-clean-data/dsR${i}-clean_L1_1.fq.gz \
        -p ../1-clean-data/dsR${i}-clean_L1_2.fq.gz \
        ../raw-data/dsR${i}_L1_1.fq.gz \
        ../raw-data/dsR${i}_L1_2.fq.gz
done

