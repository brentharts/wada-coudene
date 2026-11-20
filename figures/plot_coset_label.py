"""Plot of figures/spectre_coset_label_q*.json: for phases z in the window Z_0(O) of the tile at the origin, whether the family of tiles in the same
relative lattice coset (kappa = 0) has a tile covering the far point q, and how many of the 36 cosets do.
Usage: python3 figures/plot_coset_label.py [json]      -> figures/spectre_coset_label.png"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lakes_of_wada import spectre_window as sw
fn = sys.argv[1] if len(sys.argv) > 1 else "figures/spectre_coset_label_q-6.4_0.3.json"
d = json.load(open(fn)); rows = d["rows"]
x = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows]); kid = np.array([r[2] for r in rows])
miss0 = np.array([[0, 0, 0, 0] in r[3] for r in rows]); other = np.array([(len(r[3]) > 0) and ([0, 0, 0, 0] not in r[3]) for r in rows])
npres = np.array([r[4] for r in rows])
pts, kd = sw.points("O", 0, 5)
fig, ax = plt.subplots(1, 3, figsize=(17, 5.6))
a = ax[0]
a.scatter(pts[:, 0], pts[:, 1], c=kd, s=0.3, cmap="tab20", alpha=0.25, rasterized=True)
a.scatter(x[~miss0 & ~other], y[~miss0 & ~other], s=7, c="0.35", label="all cosets covered (%d)" % int((~miss0 & ~other).sum()))
a.scatter(x[other], y[other], s=16, c="tab:blue", label="another coset missing (%d)" % int(other.sum()))
a.scatter(x[miss0], y[miss0], s=16, c="tab:red", label="coset 0 missing: the hole (%d)" % int(miss0.sum()))
a.set_aspect("equal"); a.set_title("phases z in Z_0(O), point q = (%.1f, %.1f)" % tuple(d["q"])); a.legend(loc="lower left", fontsize=8)
b = ax[1]
frac = [(miss0[kid == k].mean() if (kid == k).any() else np.nan) for k in range(13)]
cnts = [int((kid == k).sum()) for k in range(13)]
b.bar(range(13), frac, color="tab:red"); b.set_xlabel("child piece of Z_0(O) (first-map index)"); b.set_ylabel("fraction of phases with coset 0 missing")
for k in range(13): b.text(k, (frac[k] if frac[k] == frac[k] else 0) + 0.01, "n=%d" % cnts[k], ha="center", fontsize=7)
b.set_title("is the hole constant on the 13 children?")
c = ax[2]
vals, ct = np.unique(npres, return_counts=True); c.bar(vals, ct, color="0.4"); c.set_xlabel("number of the 36 relative cosets with a present tile covering q"); c.set_ylabel("phases")
c.set_title("covering tiles per phase")
plt.tight_layout(); plt.savefig("figures/spectre_coset_label.png", dpi=130); print("wrote figures/spectre_coset_label.png")
