"""Finite-time Lyapunov exponents of the Wilson-plaquette mechanics (eps = 0.4) from the section x = 0 at energies from 1 % of the barrier
height 8 / eps^4 to 1.5 times it, against the Hessian bound 2 / eps.  Also the map to the Wilson coupling: weight exp(-V/T) is
exp(N beta Re Tr U_p) with N = 2, beta = 1 / (eps^4 T), beta_W = N^2 beta = 4 / (eps^4 T).  Usage: python3 figures/analyse_lyapunov_highE.py [n] [T]"""
import sys, csv
sys.path.insert(0, "."); sys.path.insert(0, "figures")
import numpy as np
from analyse_lyapunov_bound import lyap

n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
T = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
eps = 0.4; B = 8 / eps ** 4; rows = []
print(f"eps {eps}: barrier 8/eps^4 = {B:.1f}; beta_W = 4/7 at T = 7/eps^4 = {7/eps**4:.1f} = {7/8:.3f} of the barrier; bound 2/eps = {2/eps}")
for f in (0.01, 0.25, 0.5, 0.75, 1.5):
    E = f * B
    lam, _ = lyap(eps, n, T, dt=0.0005, E=E, seed=1, ymax=np.pi / eps, renorm=200)
    print(f"E/barrier = {f}: E = {E:.1f}, Lyapunov mean {lam.mean():.2f}, max {lam.max():.2f}", flush=True)
    rows.append((eps, f, E, lam.mean(), lam.max(), 2 / eps))
with open("figures/lyapunov_highE.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["eps", "E_over_barrier", "E", "lyap_mean", "lyap_max", "hessian_bound"]); w.writerows(rows)
