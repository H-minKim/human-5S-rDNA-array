# -*- coding: utf-8 -*-
# The common ruler for binning beta.
#   Because beta = k / n_cpg, it takes only n + 1 distinct values. Cutting 0 to 1 at a fixed width such as 0.02
#     · spreads those values unevenly across the bins and distorts the shape, and
#     · drawing through the bin centres ends short of 1.0 at the right-hand end (with 8 CpG sites the last centre is 0.95).
#   The boundaries are therefore chosen so that the bin centres fall at 0, 1/n, ..., 1. Above a count of CAP it falls back to fixed bins.
import numpy as np


def value_bins(n_cpg, cap=60):
    """Returns : (edges, centers), with bin centres at the values β can take"""
    n = int(np.clip(round(float(n_cpg)), 4, cap))
    edges = np.linspace(-0.5 / n, 1 + 0.5 / n, n + 2)
    centers = np.arange(n + 1) / float(n)
    return edges, centers


def density(b, n_cpg, cap=60):
    """density (fraction / bin width) and bin centres; 0 and 1 are the centres of the two end bins"""
    edges, centers = value_bins(n_cpg, cap)
    h, e = np.histogram(np.asarray(b, dtype=float), bins=edges)
    w = np.diff(e)
    tot = h.sum()
    return (h / tot) / w if tot else h * 0.0, centers


def fraction(b, n_cpg, cap=60):
    """fraction (read x units in the bin / all) and bin centres

    Difference from `density()` : `density()` is fraction / bin width, a density of area 1.
    The axis is named `Fraction of read-units`,
    so the plotted quantity is the fraction.
    The shape is the same : within a panel the bin width `w = 1/n` is constant, y is scaled to the maximum of each panel
    and no y ticks are written.
    `density()` is kept for the other figures that use it.
    """
    edges, centers = value_bins(n_cpg, cap)
    h, _e = np.histogram(np.asarray(b, dtype=float), bins=edges)
    tot = h.sum()
    return (h / tot if tot else h * 0.0), centers


def value_index(b, n_cpg, cap=50):
    """heat-map row index : β to the bin of the nearest attainable value (with large n the number of bins is fixed at cap)"""
    n = int(round(float(n_cpg)))
    b = np.asarray(b, dtype=float)
    if n <= cap:
        return np.clip(np.round(b * n).astype(int), 0, n), n + 1
    return np.clip((b * cap).astype(int), 0, cap - 1), cap
