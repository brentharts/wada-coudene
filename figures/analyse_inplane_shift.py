"""Control for figures/analyse_three_variable.py: fraction of two-variable exit labels (R = 2, n = 100, E = 1) that change when the start is
shifted by eps inside the section, in y or in p_y (p_x from the energy).  Usage: python3 figures/analyse_inplane_shift.py"""
import sys, csv
sys.path.insert(0, ".")
import numpy as np
from lakes_of_wada import yangmills as Y

n, R, dt, tm = 100, 2.0, 0.002, 60.0
base, ys, ps = Y.section_picture(n=n, R=R, dt=dt, t_max=tm)
Yg, P = np.meshgrid(ys, ps, indexing="ij"); ok = base > 0
rows = []
for e in (0.01, 0.03, 0.1, 0.3):
    a = Y.exit_basins(Yg + e, P, R=R, dt=dt, t_max=tm); b = Y.exit_basins(Yg, P + e, R=R, dt=dt, t_max=tm)
    va = ok & (a > 0); vb = ok & (b > 0)
    ra, rb = float(np.mean(a[va] != base[va])), float(np.mean(b[vb] != base[vb]))
    print(f"shift {e}: y {ra:.4f}, p_y {rb:.4f}", flush=True); rows.append((e, ra, rb))
with open("figures/inplane_shift.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["eps", "frac_changed_y_shift", "frac_changed_py_shift"]); w.writerows(rows)
