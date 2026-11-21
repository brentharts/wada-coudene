"""Langevin version of the Higgs toy: probability of each final well from a start at rest, on the grid of figures/higgs_lam1_gamma*.npy
(every second point, 100 x 100), against the deterministic labels.  For each gamma and temperature T: mean probability of ending in a
well other than the deterministic one ('flip'), fraction of starts whose largest probability is below 0.9 ('mixed'), fraction unsettled.
Usage: python3 figures/analyse_langevin.py GAMMA T [T ...]   (env LANG_M = realisations per start, default 12)"""
import sys, os, csv, time
sys.path.insert(0, ".")
import numpy as np
from lakes_of_wada import yangmills as Y

gamma = float(sys.argv[1]); Ts = [float(t) for t in sys.argv[2:]]
M = int(os.environ.get("LANG_M", "12"))
det = np.load("figures/higgs_lam1_gamma%g.npy" % gamma).astype(int)
n = det.shape[0]
a = (np.arange(n) + 0.5) / n * 3.6 - 1.8
sub = slice(None, None, 2 if n == 200 else 1)
X, Yg = np.meshgrid(a[sub], a[sub], indexing="ij")
D = det[sub, sub]
rows = []
for T in Ts:
    t0 = time.time()
    p = Y.higgs_langevin(X, Yg, T=T, gamma=gamma, M=M, seed=7)
    ok = D > 0
    pd = np.take_along_axis(p, D[..., None], axis=-1)[..., 0]
    flip = float(np.mean(1 - pd[ok] - p[..., 0][ok]))
    mixed = float(np.mean(p[..., 1:].max(axis=-1)[ok] < 0.9))
    uns = float(np.mean(p[..., 0]))
    print(f"gamma={gamma} T={T:g}: flip {flip:.4f}  mixed {mixed:.4f}  unsettled {uns:.4f}  ({time.time()-t0:.0f}s)", flush=True)
    rows.append((gamma, T, M, X.shape[0], flip, mixed, uns))
path = "figures/langevin_analysis.csv"
new = not os.path.exists(path)
with open(path, "a", newline="") as f:
    w = csv.writer(f)
    if new:
        w.writerow(["gamma", "T", "M", "n", "flip", "mixed", "unsettled"])
    w.writerows(rows)
