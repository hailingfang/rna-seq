#!/usr/bin/env python3
"""
This script is dedicated to normalize the read count matrix by RPKM/FPKM and TPM.

Note that, the RPKM/FPKM and TPM is meaing only within the sample. For the comparion
accross the samples, please use Differential Expression Analysis software.
"""
import argparse
from biobrary.bioparse import FASTA
import numpy as np
import re


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-file", \
                        help="RNA fasta file to get the gene length", \
                        default=None)
    parser.add_argument("--gene-len", help="a tab seprated gene length file", \
                        default=None)
    parser.add_argument("--out", help="filename of output")
    parser.add_argument("countmtx", help="count matrix file")

    args = parser.parse_args()
    if (not args.fasta_file) and (not args.gene_len):
        print("At lest one of --fastq-file and --gene-len be given.")
        exit()
    return args.fasta_file, args.gene_len, args.out, args.countmtx


def parse_seqinfo(info_line):
    re_tmp = re.compile(r'\[(.+?=.+?)\]')
    finding = re_tmp.findall(info_line)
    info = {}
    for ele in finding:
        key, value = ele.split("=")
        info[key] = value
    return info


def get_gene_len(fasta_file, gene_len_file):
    gene_len_dic = {}
    if gene_len_file:
        fin = open(gene_len_file, "r")
        for line in fin:
            line = line.rstrip().split("\t")
            geneid, length = line[0], int(line[1])
            if geneid not in gene_len_dic:
                gene_len_dic[geneid] = length
            else:
                print(f"Warning, The the {geneid} already exists in the dictionary")
        return gene_len_dic
    
    else:
        fasta = FASTA(fasta_file)
        for seq in fasta:
            seqinfo = seq.seqid_append
            seqlen = seq.seqlen
            seqinfo = parse_seqinfo(seqinfo)
            gene = seqinfo["gene"]
            if gene not in gene_len_dic:
                gene_len_dic[gene] = [seqlen]
            else:
                gene_len_dic[gene].append(seqlen)

        tmp = {}
        for gene in gene_len_dic:
            tmp[gene] = int(sum(gene_len_dic[gene]) / len(gene_len_dic[gene]))
        gene_len_dic = tmp
        return gene_len_dic


def read_countmtx(countmtx):
    fin = open(countmtx, "r")
    head = fin.readline().rstrip().split("\t")
    genes = []
    data = []
    for line in fin:
        line = line.rstrip().split("\t")
        genes.append(line[0])
        data.append([int(e) for e in line[1:]])
    data = np.asarray(data, dtype="float32")
    
    return head, genes, data


def FPKM_norm(gene_len_dic, genes, data_mtx):
    data_sum = np.sum(data_mtx, axis=0)
    data_mtx = data_mtx * 1e9
    gene_len = []
    for gg in genes:
        gene_len.append(gene_len_dic[gg])
    gene_len = np.asarray(gene_len)
    data_mtx = data_mtx.T / gene_len
    data_mtx = data_mtx.T
    data_mtx = data_mtx / data_sum
    return data_mtx


def TPM_norm(gene_len_dic, genes, data_mtx):
    data_mtx = data_mtx * 1e6
    gene_len = []
    for gg in genes:
        gene_len.append(gene_len_dic[gg])
    
    data_mtx = data_mtx.T / np.asarray(gene_len, dtype=np.float32)
    data_mtx = data_mtx.T
    data_sum = np.sum(data_mtx, axis=0)
    data_sum = data_sum / 1e6
    data_mtx = data_mtx / data_sum
    return data_mtx


def output_data(head, genes, data_norm, out):
    fout = open(out, "w")
    print("\t".join(head), file=fout)
    for i, gg in enumerate(genes):
        print("\t".join([gg] + [format(e, ".6f") for e in data_norm[i]]), file=fout)
    fout.close()
    
    return


def main():
    fasta_file, gene_len_file, out, countmtx = getargs()
    gene_len_dic = get_gene_len(fasta_file, gene_len_file)

    head, genes, data = read_countmtx(countmtx)
    data_norm_fpkm = FPKM_norm(gene_len_dic, genes, data)
    output_data(head, genes, data_norm_fpkm, out + "_fpkm.tsv")
    data_norm_tpm = TPM_norm(gene_len_dic, genes, data)
    output_data(head, genes, data_norm_tpm, out + "_tpm.tsv")

    return 0


if __name__ == "__main__":
    main()
