import numpy as np
from scipy.stats import entropy
from itertools import product

def convert_hmmer(dom_file, out_file):
    """Converts HMMER --domtblout output to regular TSV"""

    header = ["protein_id", "pfam", "i_Evalue", "domain_start", "domain_end"]

    with open(dom_file, "rt") as f:
        with open(out_file, "wt") as fout:
            fout.write("\t".join(header) + "\n")
            for line in f:
                if not line.startswith("#"):
                    line = line.split()
                    row = [line[0]] + [line[4]] + [line[12]] + line[17:19]
                    fout.write("\t".join(row) + "\n")

def coerce_numeric(s):
    """Tries to coerce string to numeric and returns number if possible"""
    try:
        return float(s)
    except ValueError:
        return s

def jsd(mat, base=np.e):
    """
    Computes Janson-Shannon Divergence given a matrix with probability vectors.
    """
    dist_vec = []
    for p, q in product(mat, mat):
        p, q = np.asarray(p), np.asarray(q)

        # Normalize p, q to probabilities
        p, q = p/p.sum(), q/q.sum()
        m = 1./2*(p + q)

        dist = entropy(p, m, base=base)/2. + entropy(q, m, base=base)/2.
        dist_vec.append(dist)

    return np.array(dist_vec).reshape(len(mat), len(mat))
