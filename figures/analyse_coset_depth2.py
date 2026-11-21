"""Is 'the tile of coset 0 covering q is missing' constant on pieces of Z_0(O)?  Labels phases by the two first maps (child, grandchild) of the exact
graph-directed system, then tabulates the missing fraction per piece.  Usage: python3 figures/analyse_coset_depth2.py [json]"""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
from scipy.spatial import cKDTree
from lakes_of_wada import spectre_window as sw
fn = sys.argv[1] if len(sys.argv) > 1 else "figures/spectre_coset_label_q-6.4_0.3.json"
rows = json.load(open(fn))["rows"]
d = sw._data(); lam = float(d["lam"])
pts, lab = [], []
for k, (i, s2) in enumerate(d["O0"]):
    sub, sk = sw._level("OG"[s2], 1, 5, {})                       # depth 6 in all
    img = (d["GA"][i, 0] + sub / lam) @ d["O"][i, 0]
    pts.append(img); lab += [(k, int(j)) for j in sk]
pts = np.concatenate(pts); tree = cKDTree(pts)
Z = np.array([[r[0], r[1]] for r in rows]); L = [lab[j] for j in tree.query(Z)[1]]
miss0 = np.array([[0, 0, 0, 0] in r[3] for r in rows]); other = np.array([len(r[3]) > 0 and [0, 0, 0, 0] not in r[3] for r in rows])
cell = collections.defaultdict(list)
for l, m, o in zip(L, miss0, other): cell[l].append((m, o))
print("%-8s %5s %8s %8s" % ("piece", "n", "miss0", "other"))
pure = 0; tot = 0; mixed = []
for l in sorted(cell):
    a = np.array(cell[l]); n = len(a); f0, fo = a[:, 0].mean(), a[:, 1].mean()
    tot += n; pure += n if (f0 in (0.0, 1.0) and fo in (0.0, 1.0)) else 0
    if f0 > 0 or fo > 0: print("%-8s %5d %8.3f %8.3f%s" % (l, n, f0, fo, "" if (f0 in (0.0, 1.0) and fo in (0.0, 1.0)) else "   mixed"))
    if not (f0 in (0.0, 1.0) and fo in (0.0, 1.0)): mixed.append((l, n, round(f0, 3), round(fo, 3)))
print("phases in pieces on which both flags are constant: %d of %d (%.3f); pieces %d, mixed pieces %d" % (pure, tot, pure / tot, len(cell), len(mixed)))
