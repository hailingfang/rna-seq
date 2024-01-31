#!

import argparse
import matplotlib.pyplot as plt
import numpy as np


def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-count", help="raw count")
    parser.add_argument("--norm-count", help="normalized count")
    parser.add_argument("--stat-res", help="deseq2 result file")
    parser.add_argument("--out", default="out")

    args = parser.parse_args()

    return args.raw_count, args.norm_count, args.stat_res, args.out


def read_raw_count(raw_count):
    genes = []
    data = []
    fin = open(raw_count, "r")
    head = fin.readline().rstrip().split("\t")[1:]
    for line in fin:
        line = line.rstrip().split("\t")
        genes.append(line[0])
        data.append([float(e) for e in line[1:]])
    return head, genes, data


def read_norm_count(norm_count):
    genes = []
    data = []
    fin = open(norm_count, "r")
    head = fin.readline().rstrip().split(",")[1:]
    head = [e.strip('"') + "_norm" for e in head]
    for line in fin:
        line = line.rstrip().split(",")
        genes.append(line[0].strip('"'))
        data.append([float(e.strip('"')) for e in line[1:]])
    return head, genes, data


def read_stat_res(stat_res):
    head = []
    genes = []
    data = []
    fin = open(stat_res, "r")
    head = fin.readline().rstrip().split(",")[1:]
    head = [e.strip('"') for e in head]
    for line in fin:
        line = line.rstrip().split(",")
        genes.append(line[0].strip('"'))
        line = [e.strip('"') for e in line[1:]]
        new_line = []
        for e in line:
            if e == "NA":
                new_line.append(np.NAN)
            else:
                new_line.append(float(e))
        data.append(new_line)
    return head, genes, data


def combine_data(head_r, genes_r, data_r, head_n, genes_n, data_n, head_s, genes_s, data_s):
    head = head_r + head_n + head_s
    for ele in zip(genes_r, genes_n, genes_s):
        if len(set(ele)) != 1:
            print(ele, "not consistent")
            exit()
    genes = genes_r
    data = np.hstack([np.asarray(data_r, dtype=np.float32), np.asarray(data_n, dtype=np.float32), np.asarray(data_s, dtype=np.float32)])

    return head, genes, data


def rm_nan(genes, data):
    genes_out = []
    data_out = []

    for i, row in enumerate(data):
        if not any(np.isnan(row)):
            genes_out.append(genes[i])
            data_out.append(row)
    data_out = np.asarray(data_out, dtype=np.float32)

    return genes_out, data_out

def plot_coor(data, c1, c2, out, prefix):
    fig, ax = plt.subplots()
    ax.scatter(data.T[c1], data.T[c2], s=3)
    fig.savefig(prefix + "_" + out + ".pdf")


def plot_scatter(data, out):
    x_mean_exp = []
    y_l2fc = []
    c = []
    alp = []
    for row in data:
        mean_exp = row[14]
        l2fc = row[15]
        padj = row[19]
        if (l2fc <= -1 or l2fc >= 1) and padj < .001:
            c.append("red")
            alp.append(1)
        else:
            c.append("gray")
            alp.append(.5)
        x_mean_exp.append(mean_exp)
        y_l2fc.append(l2fc) 

    fig, ax = plt.subplots(layout="constrained")
    ax.set_xscale("log", base=10)
    ax.scatter(x_mean_exp, y_l2fc, color=c, alpha=alp, s=3, linewidths=0)
    ax.set_xlabel("mean of normalized counts, $log_{10}$")
    ax.set_ylabel("log fold change, $log_2$")
    fig.savefig(out + "_mean_l2fc.pdf")


def plot_bar(head, genes, data, item_genes, out):
    item_idx = []
    for gg in item_genes:
        item_idx.append(genes.index(gg))
    data = data[item_idx]
    exp_mean = []
    ctl_mean = []
    exp_std = []
    ctl_std = []
    padj = []
    l2fc = []
    fout = open(out + "_item_exp.tsv", "w")
    print("\t".join(["geneid"] + head), file=fout)
    i = 0
    for row in data:
        print("\t".join([item_genes[i]] + [format(e, "1.3e") for e in row]), file=fout)
        i += 1
        exp = row[7: 11]
        ctl = row[11: 14]
        exp_mean.append(np.mean(exp))
        exp_std.append(np.std(exp))
        ctl_mean.append(np.mean(ctl))
        ctl_std.append(np.std(ctl))
        padj.append(row[19])
        l2fc.append(row[15])
    fout.close()
    fig, ax = plt.subplots()
    x = np.arange(len(item_genes))
    ax.bar(x - .2, ctl_mean, width=.4, label="ctl")
    ax.errorbar(x - .2, ctl_mean, yerr=ctl_std, linestyle="")
    ax.bar(x + .2, exp_mean, width=.4, label="exp")
    ax.errorbar(x + .2, exp_mean, yerr=exp_std, linestyle="")
    y = np.max([ctl_mean, exp_mean], axis=0)
    
    y += 100
    y[1] += 300
    l2fc = [format(e, ".3f") for e in l2fc]
    for xx, yy, ff in zip(x, y, l2fc):
        tt = ax.text(xx, yy, ff, fontsize=7)
    
    y += 100
    y[1] += 800
    padj = [format(e, "1.3e") for e in padj]
    for xx, yy, pp in zip(x - .2, y, padj):
        tt = ax.text(xx, yy, pp, fontsize=7)

    ax.legend()
    ax.set_xticks(x, item_genes)
    ylim = ax.get_ylim()

    ax.set_yscale("log", base=10)

    fig.savefig(out + "_bar.pdf")


def main():
    raw_count, norm_count, stat_res, out = getargs()
    head_r, genes_r, data_r = read_raw_count(raw_count)
    head_n, genes_n, data_n = read_norm_count(norm_count)
    head_s, genes_s, data_s = read_stat_res(stat_res)

    head, genes, data = combine_data(head_r, genes_r, data_r, head_n, genes_n, data_n, head_s, genes_s, data_s)
    print(data)
    genes, data = rm_nan(genes, data)
    print(head)
    print(data.shape)
    plot_coor(data, 0, 1, out, "m1_m2")
    plot_coor(data, 3, 4, out, "c1_c2")
    plot_scatter(data, out)
    item_genes = ["nanos3", "ddx4", "dnd1", "tdrd7a", "ca15b"]
    plot_bar(head, genes, data, item_genes, out)


if __name__ == "__main__":
    main()