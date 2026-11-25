"""Theta-like horizon scales along Mixmaster orbits (lakes_of_wada.mixmaster_horizon, QINERTIA.md section 7).
Usage: python3 figures/analyse_mixmaster_horizon.py [N] [tau_budget]  -> figures/mixmaster_horizon.csv (per-epoch rows)"""
import csv, sys, time
import numpy as np
sys.path.insert(0, ".")
from lakes_of_wada import mixmaster_horizon as MH
from lakes_of_wada.mixmaster import turnaround_state

N = int(sys.argv[1]) if len(sys.argv) > 1 else 150
TB = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
rng = np.random.default_rng(5)
p = rng.uniform(-1, 1, (4 * N, 2)); _, _, V = turnaround_state(p[:, 0], p[:, 1]); P = p[V < 0][:N]
t0 = time.time()
orbs = MH.horizon_orbits(P[:, 0], P[:, 1], tau_budget=TB, tau_cap=TB + 60, thin=10, nmax=60000)
print("orbits", len(orbs), "of", len(P), "kasner at end", np.mean([o["kasner_end"] for o in orbs]), "%.0fs" % (time.time() - t0), flush=True)
rows = []
seqs = []
for oi, o in enumerate(orbs):
    if not o["kasner_end"]:
        continue
    ep = MH.epochs(o)
    if len(ep) < 3:
        continue
    seqs.append(ep)
    for k, e in enumerate(ep):
        R = e["R"]; pe = e["p"]; w = int(np.argmax(R))
        rows.append(dict(orbit=oi, epoch=k, n_epochs=len(ep), tau=e["tau"], wide_axis=w, argmax_p=int(np.argmax(pe)), p_wide=float(pe[w]),
                         R_wide=float(R[w]), kasner_wide=float(e["kasner"][w]), excess_wide=float(R[w] / e["kasner"][w]),
                         R_min=float(R.min()), a0t_wide=float(1.0 / R[w])))
with open("figures/mixmaster_horizon.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
A = lambda k: np.array([r[k] for r in rows], float)
print("orbits used", len(seqs), "epochs", len(rows))
print("wide axis == argmax p:", np.mean(A("wide_axis") == A("argmax_p")).round(3))
ex = A("excess_wide"); print("excess R_wide*(1-p) percentiles 10/50/90/99:", np.percentile(ex, [10, 50, 90, 99]).round(3))
Rw = A("R_wide"); print("R_wide percentiles 10/50/90/99:", np.percentile(Rw, [10, 50, 90, 99]).round(2))
print("a0*t on widest axis percentiles 1/50/99:", np.percentile(A("a0t_wide"), [1, 50, 99]).round(3))
ep = A("epoch")
for k in range(0, 8):
    m = ep == k
    if m.sum() > 5:
        print("epoch %d: n=%d median R_wide %.2f  median log10 R_wide/R_min %.2f  mean excess %.2f" % (k, m.sum(), np.median(Rw[m]), np.median(np.log10(Rw[m] / A("R_min")[m])), ex[m].mean()))
# does the wide axis in epoch k predict the wide axis of the last epoch?
fin = np.array([int(np.argmax(s[-1]["R"])) for s in seqs])
for k in (0, 1, 2, 3, 5):
    ok = [i for i, s in enumerate(seqs) if len(s) > k + 1]
    a = np.mean([int(np.argmax(seqs[i][k]["R"])) == fin[i] for i in ok]); print("wide axis of epoch %d == final wide axis: %.3f (n=%d)" % (k, a, len(ok)))
for k in (1, 2, 3, 5):
    ok = [i for i, s in enumerate(seqs) if len(s) > k + 1]
    a = np.mean([int(np.argmax(seqs[i][-1 - k]["R"])) == fin[i] for i in ok]); print("wide axis %d epochs before the last == final: %.3f" % (k, a))
chg = np.mean([int(np.argmax(s[k]["R"])) != int(np.argmax(s[k + 1]["R"])) for s in seqs for k in range(len(s) - 1)])
print("wide axis changes between consecutive epochs:", round(chg, 3))
