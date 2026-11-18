#!/usr/bin/env python3
"""Boundary dimension and stability islands of the Mixmaster exit basins (lakes_of_wada.mixmaster).

    python3 figures/analyse_mixmaster.py        (about 25 minutes)   writes figures/mixmaster_analysis.png, mixmaster_analysis.csv

A  uncertainty exponent: fraction f(eps) of initial states whose exit axis changes when the state is moved by eps in a random
   direction; f ~ eps^alpha, boundary dimension D = 2 - alpha.  Repeated for four step sizes: the numerical noise inflates f.
B  survival in the number of Kasner epochs (exponential = hyperbolic, no islands) and the undecided fraction against the step budget.
C  undecided pixels of a zoom after a 10 times longer integration: do they persist (islands) or exit (slow escape)?
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import ndimage as ndi

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import mixmaster as M            # noqa: E402

rng = np.random.default_rng(2026)
rows = []
log = lambda *a: print(*a, flush=True)

# ---------------------------------------------------------------- A: uncertainty exponent at several step sizes
def valid_points(n, half=1.0):
    pts = rng.uniform(-half, half, size=(3 * n, 2))
    _, _, V = M.turnaround_state(pts[:, 0], pts[:, 1])
    return pts[V < 0][:n]

N = 4000
P = valid_points(N)
EPS = np.array([0.02, 0.01, 0.005, 0.002, 0.001])
PHI = {e: rng.uniform(0, 2 * np.pi, N) for e in EPS}           # the same perturbations for every step size
fits = {}
t0 = time.time()
for h0 in (0.01, 0.005, 0.0025, 0.00125):
    nm = int(30000 * 0.01 / h0)
    lab0 = M.exit_basins(P[:, 0], P[:, 1], h0=h0, nmax=nm)
    f, ferr = [], []
    for e in EPS:
        Q = P + e * np.stack([np.cos(PHI[e]), np.sin(PHI[e])], axis=1)
        lab1 = M.exit_basins(Q[:, 0], Q[:, 1], h0=h0, nmax=nm)
        m = (lab0 > 0) & (lab1 > 0)
        k = int((lab0[m] != lab1[m]).sum())
        f.append(k / m.sum()); ferr.append(np.sqrt(max(k, 1)) / m.sum())
        rows.append(["uncertainty_h%g" % h0, e, f[-1], ferr[-1], int(m.sum())])
    floor = float("nan")
    if h0 > 0.002:
        labh = M.exit_basins(P[:, 0], P[:, 1], h0=h0 / 2, nmax=2 * nm)
        m = (lab0 > 0) & (labh > 0)
        floor = float((lab0[m] != labh[m]).mean())
        rows.append(["step_floor_h%g" % h0, h0 / 2, floor, "", int(m.sum())])
    f = np.array(f)
    alpha, c = np.polyfit(np.log(EPS), np.log(f), 1)
    fits[h0] = (f, np.array(ferr), alpha, c, floor)
    rows.append(["alpha_h%g" % h0, "", alpha, "", len(EPS)])
    log("h0 = %-8g f(eps) = %s   floor(h0 vs h0/2) = %.3f   alpha = %.3f  D = %.3f   (%.0fs)" %
        (h0, f.round(3), floor, alpha, 2 - alpha, time.time() - t0))

# ---------------------------------------------------------------- B: survival in Kasner epochs, and the undecided fraction against the step budget
NS = 30000
PS = valid_points(NS)
t0 = time.time()
labS, tauS, nep = M.exit_basins(PS[:, 0], PS[:, 1], nmax=100000, return_time=True, return_epochs=True)
dec = labS > 0
log("survival run %.0fs, no exit in %d of %d" % (time.time() - t0, int((~dec).sum()), NS))
ns = np.arange(0, 40)
Sd = np.array([np.mean(nep[dec] > n) for n in ns])
sel = (ns >= 3) & (Sd > 1e-3)
lam, c0 = np.polyfit(ns[sel], np.log(Sd[sel]), 1)
log("among decided, S(epochs n) ~ exp(%.3f n);  epochs completed by the undecided: median %d, max %d; tau reached by them: median %.0f" %
    (lam, np.median(nep[~dec]), nep[~dec].max(), np.median(tauS[~dec])))
for n, s_ in zip(ns, Sd):
    rows.append(["survival_epochs_decided", n, s_, "", int(dec.sum())])
rows.append(["epoch_decay_rate", "", lam, "", int(sel.sum())])
sub = PS[:8000]
for nm in (30000, 100000, 300000):
    l_ = M.exit_basins(sub[:, 0], sub[:, 1], nmax=nm)
    rows.append(["undecided_fraction_steps", nm, float((l_ == 0).mean()), "", len(sub)])
    log("undecided fraction after %d steps: %.4f" % (nm, (l_ == 0).mean()))

# ---------------------------------------------------------------- C: zoom, undecided pixels after a longer run
zx, zy, zh, nz = 0.30, 0.12, 0.12, 300
a_ = (np.arange(nz) + 0.5) / nz * 2 * zh
BP, BM = np.meshgrid(zx - zh + a_, zy - zh + a_, indexing="ij")
t0 = time.time()
labZ, tauZ = M.exit_basins(BP, BM, nmax=30000, return_time=True)
und = labZ == 0
log("zoom: undecided at 30000 steps: %d of %d valid (%.0fs)" % (und.sum(), (labZ >= 0).sum(), time.time() - t0))
idx = np.argwhere(und)
labL, tauL = M.exit_basins(BP[und], BM[und], nmax=300000, return_time=True)
still = labL == 0
log("rerun %d undecided pixels with 10 x more steps: still undecided %d (%.1f %%)" % (len(idx), still.sum(), 100 * still.mean()))
rows += [["zoom_undecided_30k", "", int(und.sum()), "", int((labZ >= 0).sum())], ["zoom_undecided_300k", "", int(still.sum()), "", len(idx)]]
labZ2 = labZ.copy()
for (i, j), L in zip(idx, labL):
    labZ2[i, j] = L
und2 = labZ2 == 0
comp, ncomp = ndi.label(und2, structure=np.ones((3, 3)))
sizes = np.bincount(comp.ravel())[1:] if ncomp else np.array([])
log("components of the still-undecided set: %d, largest sizes %s" % (ncomp, np.sort(sizes)[::-1][:8]))
rows.append(["zoom_components_300k", "", ncomp, "", int(sizes.max()) if ncomp else 0])
tz = tauZ.copy()

# ---------------------------------------------------------------- figure
fig, axs = plt.subplots(2, 2, figsize=(13, 10.5))
ax = axs[0, 0]
for h0, (f_, fe_, al_, c_, fl_) in fits.items():
    ax.errorbar(EPS, f_, yerr=fe_, fmt="o-", ms=4, label="step %g: alpha = %.2f (D = %.2f)" % (h0, al_, 2 - al_))
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("eps (window units, window width 2)"); ax.set_ylabel("uncertain fraction f(eps)")
ax.set_title("A: uncertainty exponent; it has not converged with the step"); ax.legend(); ax.grid(alpha=0.3)
ax = axs[0, 1]
ax.semilogy(ns, np.maximum(Sd, 1e-5), "o", ms=4, label="decided states")
ax.semilogy(ns, np.exp(c0 + lam * ns), "-", label="exp(%.3f n)" % lam)
ax.set_xlabel("number of Kasner epochs n"); ax.set_ylabel("fraction still without exit after n epochs")
ax.set_title("B: survival in epochs (30 000 states): exponential"); ax.legend(); ax.grid(alpha=0.3)
ax = axs[1, 0]
im = ax.imshow(np.log10(np.where(np.isfinite(tz), tz, np.nan)).T, origin="lower", extent=(zx - zh, zx + zh, zy - zh, zy + zh), cmap="viridis")
plt.colorbar(im, ax=ax, label="log10 exit time tau"); ax.set_title("C: exit time in the zoom of MIXMASTER.md"); ax.set_xlabel("beta+"); ax.set_ylabel("beta-")
ax = axs[1, 1]
show = np.where(labZ2 > 0, 1, np.where(labZ2 == 0, 2, 0))
ax.imshow(show.T, origin="lower", extent=(zx - zh, zx + zh, zy - zh, zy + zh), cmap=ListedColormap(["#999", "#e8e8e8", "#c0392b"]), vmin=0, vmax=2, interpolation="nearest")
ax.set_title("C: red = still no exit after 300 000 steps (%d pixels)" % und2.sum()); ax.set_xlabel("beta+"); ax.set_ylabel("beta-")
fig.tight_layout(); fig.savefig(HERE / "mixmaster_analysis.png", dpi=110); plt.close(fig)
with open(HERE / "mixmaster_analysis.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["kind", "x", "y", "err", "n"]); w.writerows(rows)
log("done")
