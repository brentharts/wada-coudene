#!/usr/bin/env python3
"""Uncertainty fraction with the persistent exit rule, keeping the undecided states in the count.

    python3 figures/analyse_persistent.py TMAX [N]      appends rows to figures/mixmaster_persistent.csv

For states P and P + eps*(cos phi, sin phi), three counts (valid pairs only):
  differ      labels differ, undecided (0) counted as a fourth label (0 and 0 agree)
  either0     at least one end undecided, or labels differ (literal 'either end undecided')
  decided     labels differ among pairs with both ends decided (the count of MIXMASTER.md 2c, biased)
"""
import csv, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import mixmaster as M            # noqa: E402

tmax = float(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
rng = np.random.default_rng(77)
p = rng.uniform(-1, 1, size=(3 * N, 2)); _, _, V = M.turnaround_state(p[:, 0], p[:, 1]); P = p[V < 0][:N]
EPS = np.array([0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]); PHI = {e: rng.uniform(0, 2 * np.pi, N) for e in EPS}
t = time.time()
l0 = M.exit_basins_persistent(P[:, 0], P[:, 1], h0=0.1, tmax=tmax)
print("tmax %g: undecided fraction U = %.3f (%.0fs)" % (tmax, (l0 == 0).mean(), time.time() - t), flush=True)
rows = [["tmax", "eps", "differ", "either0", "decided", "U", "N"]]
for e in EPS:
    Q = P + e * np.stack([np.cos(PHI[e]), np.sin(PHI[e])], 1)
    l1 = M.exit_basins_persistent(Q[:, 0], Q[:, 1], h0=0.1, tmax=tmax)
    ok = (l1 >= 0)
    d = float((l0[ok] != l1[ok]).mean()); e0 = float(((l0[ok] != l1[ok]) | (l0[ok] == 0) | (l1[ok] == 0)).mean())
    m = (l0 > 0) & (l1 > 0); dc = float((l0[m] != l1[m]).mean())
    rows.append([tmax, e, d, e0, dc, float((l0 == 0).mean()), N])
    print("eps %-7g differ %.4f   either0 %.4f   decided-only %.4f" % (e, d, e0, dc), flush=True)
f = Path(HERE / "mixmaster_persistent.csv"); new = not f.exists()
with open(f, "a", newline="") as fh:
    w = csv.writer(fh)
    w.writerows(rows if new else rows[1:])
print("(%.0fs)" % (time.time() - t))
