"""Three-variable model H = p^2/2 + (x^2 y^2 + y^2 z^2 + z^2 x^2)/2, with the section x = 0 of the two-variable study kept and z released.
Starts: the (y, p_y) grid of the two-variable section (E = 1) with (z, p_z) = (0, zeta) [mode pz] or (zeta, 0) [mode z]; p_x from the energy.
Output: fraction of states ending on a z half axis, fraction whose x/y label differs from the two-variable label, invalid fraction.
Usage: python3 figures/analyse_three_variable.py [n]"""
import sys, csv
sys.path.insert(0, ".")
import numpy as np
from lakes_of_wada import yangmills as Y

n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
rows = []
for R, tmax in ((2.0, 60.0), (3.0, 100.0)):
    dt = 0.002
    base, ys, ps = Y.section_picture(n=n, R=R, dt=dt, t_max=tmax)
    Yg, P = np.meshgrid(ys, ps, indexing="ij")
    ok = base > 0
    h = Y.exit_basins3(Yg, P, 0.0, 0.0, R=R, dt=dt / 2, t_max=tmax)
    floor = float(np.mean(h[ok] != base[ok]))
    print(f"R={R}: step floor of the 3-variable code at zeta = 0 (dt {dt} vs {dt/2}): {floor:.4f}", flush=True)
    rows.append((R, "floor", 0.0, 0.0, 0.0, floor))
    for mode in ("pz", "z"):
        for zeta in (0.01, 0.03, 0.1, 0.3):
            z, pz = (0.0, zeta) if mode == "pz" else (zeta, 0.0)
            lab = Y.exit_basins3(Yg, P, z, pz, R=R, dt=dt, t_max=tmax)
            valid = ok & (lab != -1)
            onz = float(np.mean(np.isin(lab[valid], (5, 6))))
            xy = valid & ~np.isin(lab, (5, 6)) & (lab > 0)
            diff = float(np.mean(lab[xy] != base[xy]))
            inv = float(np.mean(lab[ok] == -1))
            print(f"  {mode}={zeta}: on z axes {onz:.4f}, x/y label differs from 2-variable on {diff:.4f} of the rest, invalid {inv:.4f}", flush=True)
            rows.append((R, mode, zeta, onz, diff, inv))
with open("figures/three_variable_analysis.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["R", "mode", "zeta", "frac_on_z_axes", "frac_xy_label_differs", "frac_invalid_or_floor"]); w.writerows(rows)
