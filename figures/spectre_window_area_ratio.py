"""Area of the window Z_0(G) (the Mystic-pair class) against Z_0(O), from the depth-7 point cloud at four tolerances.
Prediction (toymodel_v8 Stage 7): G/O -> 1/(3 + sqrt15) = 0.14550.   Usage: python3 figures/spectre_window_area_ratio.py"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from lakes_of_wada import spectre_window as sw
from scipy.spatial import cKDTree
P = {S: sw.points(S, 0, 7)[0] for S in "OG"}
lo = np.min([P[S].min(0) for S in "OG"], 0) - 0.05; hi = np.max([P[S].max(0) for S in "OG"], 0) + 0.05
trees = {S: cKDTree(P[S]) for S in "OG"}
for n, tols in ((1400, (0.012, 0.006)), (2800, (0.003, 0.0015))):
    xs = np.linspace(lo[0], hi[0], n); ys = np.linspace(lo[1], hi[1], n); px = (xs[1] - xs[0]) * (ys[1] - ys[0])
    A = {(S, t): 0.0 for S in "OG" for t in tols}
    for i in range(0, n, 200):
        G = np.stack(np.meshgrid(xs, ys[i:i + 200]), -1).reshape(-1, 2)
        for S in "OG":
            d = trees[S].query(G)[0]
            for t in tols: A[(S, t)] += (d < t).sum() * px
    for t in tols: print("n %d tol %.4f  O %.4f  G %.4f  G/O %.5f" % (n, t, A[("O", t)], A[("G", t)], A[("G", t)] / A[("O", t)]))
print("prediction 1/(3+sqrt15) = %.5f" % (1 / (3 + math.sqrt(15))))
