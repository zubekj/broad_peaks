from os import path
import numpy as np
import pandas as pd
import glob

import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages

fnames = ["FOXP1", "GATAD1", "NFYB", "CHD2", "BAP18", "POL3", "P68", "ZNF274",
          "H3.1", "H3.3", "H2B_K120", "HEY1", "HNF4A", "HNF4G", "THAP1",
          "FOS-G", "HDAC8", "JUNB", "METHYLCAP", "ICE1", "FDR5", "ELL2",
          "MLL4", "MECP2", "ZC3H8", "POL_II_SERUM", "SER2P_POLII", "SIN3A"]


def translate_index(idx):
    nidx = []
    for s in idx:

        s = s.upper()

        if s[:3] == "GSM":
            s = s.split("_")[-1]
        elif s[:2] == "H9" or s[:4] == "A549" or s[:6] == "HCT116":
            s = s.split("_")[1]

        for name in fnames:
            if name in s:
                s = name

        s = s.replace(".", "")
        s = s.replace("-", "")
        s = s.replace("_", "")

        nidx.append("{0}".format(s))

    return nidx


def select_nnan(x, y):
    mask = (~np.isnan(x)) & (~np.isnan(y))
    return x[mask], y[mask]


def cosine(x, y):
    return sum(x*y)/(np.sqrt(sum(x**2))*np.sqrt(sum(y**2)))


def euclidean(x, y):
    return np.sqrt(sum((x-y)**2))


def produce_clustermap(rankings, fname, metric):
    sns.set()
    sns.set_context("paper")

    pp = PdfPages(fname)
    sns.clustermap(rankings, col_cluster=False, metric=metric)
    pp.savefig()
    pp.close()


def produce_sorted_heatmap(rankings, fname):
    sns.set()
    sns.set_context("paper")

    pp = PdfPages(fname)
    plt.figure(figsize=(10, 12))
    rankings = rankings.iloc[rankings.mean(axis=1).argsort()[::-1]]

    sns.heatmap(rankings, mask=rankings.isnull())

    pp.savefig()
    pp.close()


def produce_distance_heatmap(rankings, fname, metric):
    dist_matrix = []
    for i1 in rankings.index:
        dist_matrix.append([])
        for i2 in rankings.index:
            dist_matrix[-1].append(metric(rankings.loc[i1], rankings.loc[i2]))

    dist_matrix = pd.DataFrame(dist_matrix)
    dist_matrix.columns = rankings.index
    dist_matrix.index = rankings.index

    sns.set()
    plt.figure()
    pp = PdfPages(fname)
    cmap = sns.cubehelix_palette(as_cmap=True, reverse=True)
    sns.heatmap(dist_matrix, cmap=cmap)
    pp.savefig()
    pp.close()


def process_cell_line(path_pattern, mod_name, na_strategy="none",
                      normalize=False):
    rankings = pd.DataFrame()
    for f in glob.glob(path_pattern):

        if "MCF7" in f or "HepG2" in f:
            continue

        d = pd.read_csv(f, index_col=1, usecols=[0, 1])

        split_path = path.basename(f).split("_")
        if split_path[2][:3] == "chr":
            d.columns = ["{0}_{1}".format(split_path[0], split_path[2])]
        else:
            d.columns = [split_path[0]]

        d.index = translate_index(d.index)
        rankings = rankings.join(d, how="outer")

    rankings = rankings[rankings.notnull().sum(axis=1) > 1]

    with open("{0}_feature_overlap.txt".format(mod_name), "w") as f:
        for i1, c1 in enumerate(rankings.columns):
            for c2 in rankings.columns[i1:]:
                f.write("{0}, {1}:\n".format(c1, c2))
                f.write(str(list(rankings.index[rankings[c1].notnull()
                                                & rankings[c2].notnull()])))
                f.write("\n\n")

    pd.Series(rankings.index).to_csv('{0}_features.csv'.format(mod_name),
                                     index=False)

    if normalize:
        rankings /= rankings.max(axis=0)

    base_metric = euclidean
    metric = base_metric

    if na_strategy == "drop":
        rankings = rankings.dropna()
    elif na_strategy == "fill":
        rankings = rankings.fillna(0)
    elif na_strategy == "overlap":
        metric = lambda x, y: base_metric(*select_nnan(x, y))
    else:
        pass

    rankings = rankings.sort_index(axis=0)
    rankings = rankings.sort_index(axis=1)

    produce_sorted_heatmap(rankings, '{0}_sorted_heatmap.pdf'.format(mod_name))
    #produce_distance_heatmap(rankings, '{0}_heatmap.pdf'.format(mod_name), metric)


process_cell_line("../data/processed_peaks/feature_rankings/*H3K4me3_chr*",
                  "H3K4me3_chr1_chr20")
#process_cell_line("../data/processed_peaks/feature_rankings/*H3K27ac*",
#                 "H3K27ac")
#process_cell_line("../data/processed_peaks/feature_rankings/*H3K4me3*",
#                  "H3K4me3")
#process_cell_line("../data/processed_peaks/feature_rankings/*H3K27ac*",
#                  "H3K27ac_dropna", "drop")
#process_cell_line("../data/processed_peaks/feature_rankings/*H3K4me3*",
#                  "H3K4me3_dropna", "drop")
#process_cell_line("../data/processed_peaks/feature_rankings/*H3K27ac*",
#                  "H3K27ac_overlap", "overlap")
#process_cell_line("../data/processed_peaks/feature_rankings/*H3K4me3*",
#                  "H3K4me3_overlap", "overlap")
