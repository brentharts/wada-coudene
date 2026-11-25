"""Henon-Heiles exit basins with the quantised-inertia laws (lakes_of_wada.quantised_inertia, QINERTIA.md).
Usage: python3 figures/analyse_quantised_inertia.py   -> figures/qinertia_analysis.csv  (about 10 minutes)"""
import csv, sys, time
import numpy as np
sys.path.insert(0, ".")
from lakes_of_wada import quantised_inertia as Q, wadatest as W

N, E = 160, 0.25
rows = []


def record(**kw):
    rows.append(kw)
    print(kw, flush=True)


ref = Q.hh_labels(n=N, energy=E, law="newton")
cases = [("newton", 0.0), ("mond", 0.005), ("mond", 0.02), ("mond", 0.05), ("mond", 0.2), ("qi1", 0.02), ("qi1", 0.2)]
for law, a0 in cases:
    t = time.time()
    lab = ref if law == "newton" else Q.hh_labels(n=N, energy=E, law=law, a0=a0)
    both = (ref > 0) & (lab > 0)
    decided = lab > 0
    wp = W.wada_profile(lab, radii=(1, 2, 4, 8))["profile"]
    # same test on decided pixels only (undecided relabelled invalid), so trapped regions do not count as boundary
    lab_dec = np.where(lab == 0, -1, lab)
    wpd = W.wada_profile(lab_dec, radii=(1, 2, 4, 8))["profile"]
    fine = Q.hh_labels(n=N, energy=E, law=law, a0=a0, dt=0.025)
    fb = (lab > 0) & (fine > 0)
    record(kind="picture", law=law, a0=a0, n=N, E=E, undecided=float(np.mean(lab[lab >= 0] == 0)),
           differs_from_newton=float(np.mean(ref[both] != lab[both])),
           wada_r2=wp[2], wada_r4=wp[4], wada_r8=wp[8], wada_dec_r2=wpd[2], wada_dec_r4=wpd[4], wada_dec_r8=wpd[8],
           step_floor=float(np.mean(lab[fb] != fine[fb])), seconds=round(time.time() - t, 1))
    if a0 <= 0.05:
        u = Q.hh_uncertainty(energy=E, law=law, a0=a0, n_samples=5000)
        record(kind="uncertainty", law=law, a0=a0, alpha=u["alpha"], D=u["dimension"], **{f"f_{e:g}": v for e, v in u["f"].items()})
for a0 in (0.2,):
    for tm in (200.0, 800.0):
        lab = Q.hh_labels(n=60, energy=E, law="mond", a0=a0, t_max=tm)
        record(kind="undecided_vs_tmax", law="mond", a0=a0, t_max=tm, undecided=float(np.mean(lab[lab >= 0] == 0)))
keys = sorted({k for r in rows for k in r})
with open("figures/qinertia_analysis.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)
