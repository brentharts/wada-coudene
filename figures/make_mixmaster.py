#!/usr/bin/env python3
"""Mixmaster exit basins (lakes_of_wada.mixmaster): figures/mixmaster.png and figures/mixmaster_profiles.csv.

    python3 figures/make_mixmaster.py        (about 3 minutes)
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import mixmaster as M, wadatest as W      # noqa: E402

COLS = ["#999999", "#ffffff", "#d1495b", "#2e86ab", "#edae49"]      # invalid, undecided, axis a, b, c
cmap = ListedColormap(COLS)

fig, axs = plt.subplots(1, 3, figsize=(17, 5.6))
rows = []
HALF = 1.2
for n in (100, 200, 400):
    pic = M.mixmaster_picture(n=n, half=HALF)
    prof = W.wada_profile(pic.label, (1, 2, 4, 8, 16))
    rad = W.wada_radius(pic.label, percentiles=(50, 90, 99))["radius"]
    rows.append([n, prof["K"], prof["boundary_pixels"], int((pic.label == 0).sum())]
                + [round(prof["profile"][r], 4) for r in (1, 2, 4, 8, 16)]
                + [round(rad[q] * 2 * HALF / n, 4) for q in (50, 90, 99)])
    if n == 400:
        axs[0].imshow(pic.label.T, origin="lower", extent=(-HALF, HALF, -HALF, HALF), cmap=cmap, vmin=-1, vmax=3,
                      interpolation="nearest")
axs[0].set_title("exit axis of the maximum-expansion state, 400 x 400\n(grey: V >= 0, not an initial state; white: no exit)")
axs[0].set_xlabel("beta+"); axs[0].set_ylabel("beta-")
axs[0].axhline(0, color="k", lw=0.5, ls=":")
axs[0].text(-1.15, 0.04, "dotted: b = c (Berger / Hopf line)", fontsize=8)

# zoom
zx, zy, zh = 0.30, 0.12, 0.12
n = 300
ax_ = (np.arange(n) + 0.5) / n * 2 * zh
BP, BM = np.meshgrid(zx - zh + ax_, zy - zh + ax_, indexing="ij")
lab = M.exit_basins(BP, BM)
axs[1].imshow(lab.T, origin="lower", extent=(zx - zh, zx + zh, zy - zh, zy + zh), cmap=cmap, vmin=-1, vmax=3,
              interpolation="nearest")
axs[1].set_title("zoom, 300 x 300 on a window 10 times smaller")
axs[1].set_xlabel("beta+"); axs[1].set_ylabel("beta-")
zp = W.wada_profile(lab, (1, 2, 4, 8, 16))
rows.append(["zoom300", zp["K"], zp["boundary_pixels"], int((lab == 0).sum())] + [round(zp["profile"][r], 4) for r in (1, 2, 4, 8, 16)] + ["", "", ""])

for r in rows[:3]:
    axs[2].plot((1, 2, 4, 8, 16), r[4:9], "o-", label="%d x %d" % (r[0], r[0]))
axs[2].set_xscale("log", base=2); axs[2].set_ylim(0, 1.02)
axs[2].set_xlabel("radius r (pixels)"); axs[2].set_ylabel("fraction of boundary pixels seeing all 3 exits within r")
axs[2].set_title("Wada profile"); axs[2].legend(); axs[2].grid(alpha=0.3)
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(HERE / "mixmaster.png", dpi=120); plt.close(fig)
with open(HERE / "mixmaster_profiles.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["n", "K", "boundary_pixels", "undecided", "r1", "r2", "r4", "r8", "r16", "radius_p50", "radius_p90", "radius_p99"])
    w.writerows(rows)
print(open(HERE / "mixmaster_profiles.csv").read())
