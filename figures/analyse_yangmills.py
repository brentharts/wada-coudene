#!/usr/bin/env python3
"""Exit basins of homogeneous SU(2) Yang-Mills mechanics (lakes_of_wada.yangmills): structure, Wada profile, uncertainty exponent.

    python3 figures/analyse_yangmills.py        (about 8 minutes)   writes figures/yangmills.png and figures/yangmills_analysis.csv

For each exit radius R: the label picture on the section x = 0 (y in [-2, 2], p_y in (-sqrt 2, sqrt 2), E = 1); neighbour agreement of the
labels against a shuffled-label control; the Wada profile of the picture and of the shuffled picture; the step floor (labels at dt 0.001 against
0.0005); the uncertainty fraction f(eps) with its power-law fit.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import yangmills as Y, wadatest as W     # noqa: E402

rng = np.random.default_rng(11)
rows = [["kind", "R", "x", "value", "n"]]
log = lambda *a: print(*a, flush=True)
Rs = (2.0, 3.0, 4.0)
tmax = {2.0: 60, 3.0: 100, 4.0: 150}
EPS = np.array([1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4])
N = 6000
Y0 = rng.uniform(-2, 2, N); P0 = rng.uniform(-np.sqrt(2) * .999, np.sqrt(2) * .999, N); PH = rng.uniform(0, 2 * np.pi, N)
pics = {}
for R in Rs:
    t0 = time.time()
    lab, _, _ = Y.section_picture(n=200, R=R, t_max=tmax[R], dt=0.002)
    pics[R] = lab
    v = lab[lab > 0]
    ok = (lab[:-1] > 0) & (lab[1:] > 0)
    agree = float(np.mean(lab[:-1][ok] == lab[1:][ok]))
    sh = rng.permutation(v); null = float(np.mean(sh[:-1] == sh[1:]))
    wp = W.wada_profile(lab, radii=(1, 2, 3, 4, 8))["profile"]
    lab_s = lab.copy(); m = lab_s > 0; lab_s[m] = rng.permutation(lab_s[m])
    wps = W.wada_profile(lab_s, radii=(1, 2, 3, 4, 8))["profile"]
    a1 = Y.exit_basins(Y0, P0, R=R, dt=0.001, t_max=tmax[R]); a2 = Y.exit_basins(Y0, P0, R=R, dt=0.0005, t_max=tmax[R])
    m = (a1 > 0) & (a2 > 0); floor = float(np.mean(a1[m] != a2[m]))
    f = []
    for e in EPS:
        l1 = Y.exit_basins(Y0 + e * np.cos(PH), P0 + e * np.sin(PH), R=R, dt=0.0005, t_max=tmax[R])
        m = (a2 > 0) & (l1 > 0); f.append(float(np.mean(a2[m] != l1[m])))
    f = np.array(f)
    sel = f > 1.5 * floor
    al, c = np.polyfit(np.log(EPS[sel]), np.log(f[sel]), 1) if sel.sum() >= 3 else (np.nan, np.nan)
    log("R=%g: undecided %.3f; neighbour agreement %.3f (shuffled %.3f); Wada profile r=1..4,8: %s; shuffled: %s; step floor (dt .001 vs .0005) %.3f"
        % (R, np.mean(lab == 0), agree, null, [round(wp[r], 3) for r in (1, 2, 3, 4, 8)], [round(wps[r], 3) for r in (1, 2, 3, 4, 8)], floor))
    log("     f(eps) = %s ; alpha = %.3f (eps with f > 1.5 floor: %d), D = %.3f   (%.0fs)" % (f.round(3), al, sel.sum(), 2 - al, time.time() - t0))
    rows += [["agree", R, "", agree, ""], ["agree_shuffled", R, "", null, ""], ["step_floor", R, "", floor, N], ["alpha", R, "", al, int(sel.sum())]]
    rows += [["profile_r%d" % r, R, "", wp[r], ""] for r in (1, 2, 3, 4, 8)] + [["profile_shuffled_r%d" % r, R, "", wps[r], ""] for r in (1, 2, 3, 4, 8)]
    rows += [["f_eps", R, e, fv, N] for e, fv in zip(EPS, f)]
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
cm = ListedColormap(["#ddd", "#d1495b", "#2e86ab", "#edae49", "#5b8c5a"])
for ax, R in zip(axs, Rs):
    ax.imshow(pics[R].T, origin="lower", extent=(-2, 2, -np.sqrt(2), np.sqrt(2)), cmap=cm, vmin=0, vmax=4, aspect="auto")
    ax.set_title("exit radius R = %g (red +x, blue +y, yellow -x, green -y)" % R, fontsize=9); ax.set_xlabel("y"); ax.set_ylabel("p_y")
fig.tight_layout(); fig.savefig(HERE / "yangmills.png", dpi=110)
with open(HERE / "yangmills_analysis.csv", "w", newline="") as fh:
    csv.writer(fh).writerows(rows)
