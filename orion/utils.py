import numpy as np
from scipy.stats import entropy
from itertools import product

def coerce_numeric(s):
    """Tries to coerce string to numeric and returns number if possible"""
    try:
        return float(s)
    except ValueError:
        return s

def jsd_pairwise(p, q, base=2):
    """
    Computes Janson-Shannon Divergence given two probability vectors p and q.
    """
    p, q = np.asarray(p), np.asarray(q)

    # Normalize p, q to probabilities
    p, q = p/p.sum(), q/q.sum()
    m = 1./2*(p + q)
    return entropy(p, m, base=base)/2. + entropy(q, m, base=base)/2.

def jsd(mat, base=2):
    """
    Computes Janson-Shannon Divergence given a matrix with probability vectors.
    """
    dist_vec = []
    for p, q in product(mat, mat):
        dist = jsd_pairwise(p, q)
        dist_vec.append(dist)
    return np.array(dist_vec).reshape(len(mat), len(mat))

def tanimoto_pairwise(p, q):
    """
    Computes the cosine tanimoto coefficient given two probability vectors p and q.
    """
    pq = p * q
    p_square = p ** 2
    q_square = q ** 2
    return pq / (p_square + q_square - pq)


def tanimoto(mat):
    """
    Computes the cosine tanimoto coefficient a matrix with probability vectors.
    """
    dist_vec = []
    for p, q in product(mat, mat):
        dist = tanimoto_pairwise(p, q)
        dist_vec.append(dist)
    return np.array(dist_vec).reshape(len(mat), len(mat))
