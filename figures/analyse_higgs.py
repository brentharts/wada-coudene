#!/usr/bin/env python3
"""Basins of the damped gauge-Higgs toy V = x^2 y^2/2 + lam (x^2 + y^2 - 1)^2 / 4 (lakes_of_wada.yangmills.higgs_basins).

    python3 figures/analyse_higgs.py LAM GAMMA [GAMMA ...]     appends rows to figures/higgs_analysis.csv

For each friction GAMMA: the basin picture (released from rest on [-1.8, 1.8]^2, 200 x 200), neighbour agreement against shuffled labels, the Wada
profile of the picture and of the shuffled picture, the radius needed to see all four basins at n = 200 and 400, the step floor (dt 0.01 against 0.005)
and the uncertainty fraction f(eps) with its log-log fit.
"""
import csv, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import yangmills as Y, wadatest as W     # noqa: E402

import os
FAST = bool(os.environ.get('HIGGS_FAST'))      # grids 100 and 200 (instead of 200 and 400), 2000 states for the uncertainty fraction
N0, N1 = (100, 200) if FAST else (200, 400)
lam = float(sys.argv[1]); gammas = [float(g) for g in sys.argv[2:]]
rng = np.random.default_rng(5)
EPS = np.array([3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4]); N = 2000 if FAST else 4000
X0 = rng.uniform(-1.8, 1.8, N); Y0 = rng.uniform(-1.8, 1.8, N); PH = rng.uniform(0, 2 * np.pi, N)
rows = []
for g in gammas:
    t0 = time.time(); tm = max(25.0 / g, 40.0)
    kw = dict(lam=lam, gamma=g, t_max=tm)
    lab, _ = Y.higgs_picture(n=N0, dt=0.01, **kw)
    v = lab[lab > 0]; ok = (lab[:-1] > 0) & (lab[1:] > 0)
    agree = float(np.mean(lab[:-1][ok] == lab[1:][ok])); sh = rng.permutation(v); null = float(np.mean(sh[:-1] == sh[1:]))
    prof = W.wada_profile(lab, radii=(1, 2, 3, 4, 8))["profile"]
    ls = lab.copy(); m = ls > 0; ls[m] = rng.permutation(ls[m]); profs = W.wada_profile(ls, radii=(1, 2, 3, 4, 8))["profile"]
    rad = W.wada_radius(lab, percentiles=(50, 90, 99))["radius"]
    lab4, _ = Y.higgs_picture(n=N1, dt=0.01, **kw); rad4 = W.wada_radius(lab4, percentiles=(50, 90, 99))["radius"]
    a1 = Y.higgs_basins(X0, Y0, dt=0.01, **kw); a2 = Y.higgs_basins(X0, Y0, dt=0.005, **kw)
    m = (a1 > 0) & (a2 > 0); floor = float(np.mean(a1[m] != a2[m]))
    f = []
    for e in EPS:
        l1 = Y.higgs_basins(X0 + e * np.cos(PH), Y0 + e * np.sin(PH), dt=0.005, **kw)
        m = (a2 > 0) & (l1 > 0); f.append(float(np.mean(a2[m] != l1[m])))
    f = np.array(f); sel = f > 1.5 * floor
    al = np.polyfit(np.log(EPS[sel]), np.log(f[sel]), 1)[0] if sel.sum() >= 3 else float("nan")
    print("lam %g gamma %g: undecided %.3f; agreement %.3f (shuffled %.3f); profile r=1,2,3,4,8 %s (shuffled %s); radius px n=%d %s, n=%d %s; step floor %.3f; f %s alpha %.3f D %.3f (%d pts above floor)  (%.0fs)" %
          (lam, g, np.mean(lab == 0), agree, null, [round(prof[r], 3) for r in (1, 2, 3, 4, 8)], [round(profs[r], 3) for r in (1, 2, 3, 4, 8)],
           N0, [round(rad[k], 1) for k in (50, 90, 99)], N1, [round(rad4[k], 1) for k in (50, 90, 99)], floor, f.round(3), al, 2 - al, sel.sum(), time.time() - t0), flush=True)
    rows.append([lam, g, float(np.mean(lab == 0)), agree, null, floor, al, int(sel.sum())] + [prof[r] for r in (1, 2, 3, 4, 8)] + [profs[r] for r in (1, 2, 3, 4, 8)] + [rad[k] for k in (50, 90, 99)] + [rad4[k] for k in (50, 90, 99)] + list(f))
    np.save(HERE / ("higgs_lam%g_gamma%g.npy" % (lam, g)), lab.astype(np.int8))
f_ = HERE / "higgs_analysis.csv"; new = not f_.exists()
with open(f_, "a", newline="") as fh:
    w = csv.writer(fh)
    if new:
        w.writerow(["lam", "gamma", "undecided", "agree", "agree_shuffled", "step_floor", "alpha", "n_fit"] + ["prof_r%d" % r for r in (1, 2, 3, 4, 8)] + ["shuf_r%d" % r for r in (1, 2, 3, 4, 8)] + ["radA_p%d" % k for k in (50, 90, 99)] + ["radB_p%d" % k for k in (50, 90, 99)] + ["f_%g" % e for e in EPS])
    w.writerows(rows)
