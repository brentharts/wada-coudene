#!/usr/bin/env python3
"""Figure of the damped gauge-Higgs basins from the arrays written by analyse_higgs.py (lam = 1)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
HERE = Path(__file__).resolve().parent
cm = ListedColormap(["#ddd", "#d1495b", "#2e86ab", "#edae49", "#5b8c5a"])
gs = [1.0, 0.3, 0.1, 0.05, 0.02]
fig, axs = plt.subplots(1, len(gs), figsize=(3.6 * len(gs), 3.9))
for ax, g in zip(axs, gs):
    a = np.load(HERE / ("higgs_lam1_gamma%g.npy" % g))
    ax.imshow(a.T, origin="lower", extent=(-1.8, 1.8, -1.8, 1.8), cmap=cm, vmin=0, vmax=4, interpolation="nearest")
    ax.set_title("friction %g (%d x %d)" % (g, a.shape[0], a.shape[0]), fontsize=9); ax.set_xlabel("x0"); ax.set_ylabel("y0")
fig.suptitle("final well from rest: red +x, blue +y, yellow -x, green -y   (V = x^2 y^2/2 + (x^2 + y^2 - 1)^2/4)", fontsize=9)
fig.tight_layout(); fig.savefig(HERE / "higgs_basins.png", dpi=105)
