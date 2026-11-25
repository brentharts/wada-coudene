"""Are the undecided orbits of the mond law at a0 = 1/radius bounded, or just slow?
Usage: python3 figures/analyse_qi_trapped.py  -> figures/qinertia_trapped.csv"""
import csv, sys
import numpy as np
sys.path.insert(0, ".")
from lakes_of_wada import quantised_inertia as Q
from lakes_of_wada.scattering import integrate_exits
Q.integrate_exits = integrate_exits
from lakes_of_wada.scattering import harmonic_section

E, A0, R = 0.25, 0.2, 5.0
y0, valid = harmonic_section(3, 80, E)
s0 = y0[valid]
sysm = Q.hh_system(3, R, "mond", A0)
rng = np.random.default_rng(0)
lab, _, _ = Q.integrate_exits(sysm, s0, 0.05, 200.0) if hasattr(Q, "integrate_exits") else (None, None, None)
und = s0[lab == 0]
pick = und[rng.choice(len(und), size=min(300, len(und)), replace=False)]
print("undecided starts at t=200:", len(und), "of", len(s0), "using", len(pick), flush=True)


def run(s, dt, t_max, chunk=200.0):
    """Integrate in chunks; return exit time (nan if none), max radius, and r at the end."""
    y = s.copy(); m = len(y); texit = np.full(m, np.nan); rmax = np.hypot(y[:, 0], y[:, 1]); alive = np.ones(m, bool)
    t = 0.0
    while t < t_max and alive.any():
        l, _, f = Q.integrate_exits(sysm, y[alive], dt, chunk)
        idx = np.where(alive)[0]
        y[idx] = f
        # integrate_exits freezes exits at the step they happen; label>0 means exited in this chunk
        texit[idx[l > 0]] = t + chunk  # upper bound on the exit time (chunk resolution)
        rmax[idx] = np.maximum(rmax[idx], np.hypot(f[:, 0], f[:, 1]))
        alive[idx[l > 0]] = False
        t += chunk
    return texit, rmax, np.hypot(y[:, 0], y[:, 1])


rows = []
for dt, tm in ((0.05, 2000.0), (0.02, 2000.0)):
    te, rm, rend = run(pick, dt, tm)
    ex = np.isfinite(te)
    q = np.nanpercentile(te[ex], [25, 50, 75, 100]) if ex.any() else [np.nan] * 4
    rows.append(dict(dt=dt, t_max=tm, n=len(pick), exited=float(ex.mean()), exit_t_quartiles=";".join("%.0f" % v for v in q),
                     rmax_trapped_median=float(np.median(rm[~ex])) if (~ex).any() else np.nan,
                     rmax_trapped_max=float(np.max(rm[~ex])) if (~ex).any() else np.nan,
                     n_still_inside=int((~ex).sum())))
    print(rows[-1], flush=True)
with open("figures/qinertia_trapped.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
