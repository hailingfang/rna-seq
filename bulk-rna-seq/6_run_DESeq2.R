#! /usr/bin/env Rscript
library(DESeq2)

args <- commandArgs()
print(args)

mtx <- read.csv(args[6], sep="\t", row.names="geneid")
design <- read.csv(args[7], sep="\t", row.names="sample")
design$condition <- factor(design$condition, levels=c("ctl", "exp"))

dds <- DESeqDataSetFromMatrix(countData=mtx, colData=design, design =~condition)
dds <- DESeq(dds)
res <- results(dds)
norm_count <- counts(dds, normalized=TRUE)
write.csv(res, file=args[8])
write.csv(norm_count, file=args[9])
