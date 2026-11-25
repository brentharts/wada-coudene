"""Does the horizon scale of epoch 0 converge as the orbit is followed further towards the singularity?
Usage: python3 figures/analyse_mixmaster_horizon_budget.py  -> figures/mixmaster_horizon_budget.csv"""
import csv, sys
import numpy as np
sys.path.insert(0, ".")
from lakes_of_wada import mixmaster_horizon as MH
from lakes_of_wada.mixmaster import turnaround_state

rng = np.random.default_rng(5)
p = rng.uniform(-1, 1, (600, 2)); _, _, V = turnaround_state(p[:, 0], p[:, 1]); P = p[V < 0][:150]
res = {}
for tb in (10.0, 20.0, 40.0, 80.0, 160.0):
    orbs = MH.horizon_orbits(P[:, 0], P[:, 1], tau_budget=tb, tau_cap=tb + 80, thin=10, nmax=120000)
    out = {}
    for o in orbs:
        if not o["kasner_end"]:
            continue
        ep = MH.epochs(o)
        if ep:
            out[tuple(np.round(o["start"][0], 9))] = (ep[0]["R"], len(ep))
    res[tb] = out
    print("budget", tb, "orbits", len(out), flush=True)
keys = set.intersection(*[set(r) for r in res.values()])
print("common orbits", len(keys))
rows = []
ref = res[160.0]
for tb in (10.0, 20.0, 40.0, 80.0):
    ratio = np.array([np.max(res[tb][k][0]) / np.max(ref[k][0]) for k in keys])
    ratio_sum = np.array([np.sum(res[tb][k][0]) / np.sum(ref[k][0]) for k in keys])
    rows.append(dict(budget=tb, n=len(keys), median_ratio_Rmax=float(np.median(ratio)), p10=float(np.percentile(ratio, 10)), p90=float(np.percentile(ratio, 90)),
                     frac_within_5pct=float(np.mean(np.abs(ratio - 1) < 0.05)), median_ratio_sumR=float(np.median(ratio_sum))))
    print(rows[-1], flush=True)
with open("figures/mixmaster_horizon_budget.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
