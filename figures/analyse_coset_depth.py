"""Purity of the flags 'coset 0 missing' / 'another coset missing' on pieces of Z_0(O) of depth 1, 2, 3, against a shuffle null that keeps the rate within each child.
Usage: python3 figures/analyse_coset_depth.py [json]"""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from scipy.spatial import cKDTree
from lakes_of_wada import spectre_window as sw
fn = sys.argv[1] if len(sys.argv) > 1 else "figures/spectre_coset_label_q-6.4_0.3.json"
rows = json.load(open(fn))["rows"]; d = sw._data(); lam = float(d["lam"])
def lab_points(S, p, depth):
    if depth == 0: return np.zeros((1, 2)), [()]
    pts, labs = [], []
    for k, (i, s2) in enumerate(d[f"{S}{p}"]):
        sub, sl = lab_points("OG"[s2], 1 - p, depth - 1)
        pts.append((d["GA"][i, p] + sub / lam) @ d["O"][i, p]); labs += [(k,) + t for t in sl]
    return np.concatenate(pts), labs
pts, labs = lab_points("O", 0, 6); tree = cKDTree(pts)
Z = np.array([[r[0], r[1]] for r in rows]); L = [labs[j] for j in tree.query(Z)[1]]
flag = np.array([2 * ([0, 0, 0, 0] in r[3]) + (len(r[3]) > 0 and [0, 0, 0, 0] not in r[3]) for r in rows])   # 0 none, 1 other, 2 coset 0
def purity(depth, fl):
    cells = collections.defaultdict(list)
    for l, f in zip(L, fl): cells[l[:depth]].append(f)
    return sum(len(v) for v in cells.values() if len(set(v)) == 1) / len(fl), len(cells)
rng = np.random.default_rng(0); kids = np.array([l[0] for l in L])
print("depth  pieces(occupied)  purity   shuffle-null mean +- sd (rates kept within each child)")
for dep in (1, 2, 3):
    p, nc = purity(dep, flag); nulls = []
    for _ in range(300):
        fl = flag.copy()
        for k in set(kids): idx = np.where(kids == k)[0]; fl[idx] = rng.permutation(flag[idx])
        nulls.append(purity(dep, fl)[0])
    print("%3d  %8d  %9.3f   %.3f +- %.3f" % (dep, nc, p, np.mean(nulls), np.std(nulls)))
