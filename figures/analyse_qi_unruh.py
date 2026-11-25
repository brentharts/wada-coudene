"""Unruh-type velocity noise on Henon-Heiles with and without the inertia law (QINERTIA.md section 5).
Usage: python3 figures/analyse_qi_unruh.py -> figures/qinertia_unruh.csv"""
import csv, sys
import numpy as np
sys.path.insert(0, ".")
from lakes_of_wada import quantised_inertia as Q
from lakes_of_wada.scattering import harmonic_section, integrate_exits

E, n, REAL = 0.25, 50, 12
y0, valid = harmonic_section(3, n, E)
starts = y0[valid]
rows = []
for law, a0 in (("newton", 0.0), ("mond", 0.05)):
    ref, _, _ = integrate_exits(Q.hh_system(3, 5.0, law, a0), starts, 0.05, 200.0)
    for mode in ("unruh", "const"):
        for D0 in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3):
            p = Q.hh_noisy_labels(starts, law, a0, D0, mode, n_real=REAL)
            dec = ref > 0
            flip = float(np.mean(1.0 - p[np.arange(len(ref)), ref][dec]))
            mixed = float(np.mean(p[dec, 1:].max(axis=1) < 0.9))
            rows.append(dict(law=law, a0=a0, mode=mode, D0=D0, flip=flip, mixed=mixed, no_exit=float(p[:, 0].mean()), n_starts=len(starts), n_real=REAL))
            print(rows[-1], flush=True)
with open("figures/qinertia_unruh.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
