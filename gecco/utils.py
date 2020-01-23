from itertools import product

import numpy
from scipy.stats import entropy
from scipy.spatial.distance import jensenshannon

def coerce_numeric(s):
    """Tries to coerce string to numeric and returns number if possible"""
    try:
        return float(s)
    except ValueError:
        return s

def jsd_pairwise(p, q, base=2):
    """
    Computes Janson-Shannon Distance given two probability vectors p and q.
    """
    return jensenshannon(p, q, base=2)

def jsd(mat, base=2):
    """
    Computes Janson-Shannon Divergence given a matrix with probability vectors.
    """
    dist_vec = []
    for p, q in product(mat, mat):
        dist = jsd_pairwise(p, q)
        dist_vec.append(dist)
    return numpy.array(dist_vec).reshape(len(mat), len(mat))

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
    return numpy.array(dist_vec).reshape(len(mat), len(mat))
