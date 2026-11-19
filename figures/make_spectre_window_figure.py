"""The four exact Spectre window pieces Z_p(S) (S = O, G; p = 0, 1) coloured by the child piece (first map used) -> figures/spectre_windows.png.
Rasterised at depth 7 from lakes_of_wada/spectre_window.py (maps of brentharts/spectre boundary_certificate/spectre_window.py)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lakes_of_wada import spectre_window as sw

DEPTH = 7
fig, axs = plt.subplots(1, 4, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.3, 1.0, 1.1, 0.5]})
for ax, (S, p) in zip(axs, [("O", 0), ("O", 1), ("G", 0), ("G", 1)]):
    pts, kid = sw.points(S, p, DEPTH)
    ax.scatter(pts[:, 0], pts[:, 1], c=kid, cmap="tab20", s=0.15, marker=".", linewidths=0, rasterized=True)
    ax.set_aspect("equal"); ax.set_title(f"$Z_{p}({S})$, {sw.count(S, p, DEPTH)} points, {kid.max() + 1} children", fontsize=9)
fig.suptitle("Spectre window pieces (graph-directed attractor, depth %d); colour = child piece; seam dimension 1.1110" % DEPTH, fontsize=10)
fig.tight_layout()
fig.savefig(Path(__file__).with_name("spectre_windows.png"), dpi=110)
