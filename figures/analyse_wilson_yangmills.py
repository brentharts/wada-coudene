"""How far do the first-passage exit labels of x^2 y^2 / 2 change when the potential is the Wilson plaquette of spacing eps?
Output: fraction of (y, p_y) states (n x n grid, E = 1) whose label differs from the continuum label at the same step, and the step floor
(continuum dt against dt/2) for scale.  Usage: python3 figures/analyse_wilson_yangmills.py [n]"""
import sys, csv
sys.path.insert(0, ".")
import numpy as np
from lakes_of_wada import yangmills as Y

n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
rows = []
for R, tmax in ((2.0, 60.0), (3.0, 100.0)):
    dt = 0.002
    base, ys, ps = Y.section_picture(n=n, R=R, dt=dt, t_max=tmax)
    half, _, _ = Y.section_picture(n=n, R=R, dt=dt / 2, t_max=tmax)
    ok = base > 0
    floor = float(np.mean(base[ok] != half[ok]))
    print(f"R={R} continuum step floor (dt {dt} vs {dt/2}) = {floor:.4f}", flush=True)
    rows.append((R, "floor", floor))
    for eps in (0.05, 0.1, 0.2, 0.4):
        lab, _, _ = Y.section_picture(n=n, R=R, dt=dt, t_max=tmax, eps=eps)
        both = ok & (lab > 0)
        d = float(np.mean(base[both] != lab[both]))
        und = float(np.mean(lab[ok] == 0))
        print(f"  eps={eps}: label differs from continuum on {d:.4f} of states, undecided {und:.4f}", flush=True)
        rows.append((R, eps, d))
with open("figures/wilson_yangmills_analysis.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["R", "eps_or_floor", "fraction_labels_differ"]); w.writerows(rows)
