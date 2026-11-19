"""Trivial-field Hessian of the Wilson action on the Spectre lift lattice A2 + A2 (brentharts/spectre, spectre_wilson.py, model H), the
analogue of Romanov (2026) Proposition 1.3: at the trivial field -Hess S_0(v, v) = sum_p w_p (d xi)_p^2, so its largest value over the unit
sphere of link variables (the top eigenvalue of the 6 x 6 Bloch matrix M(k) of `spectre_wilson.hessian`, over k) is a lower bound for
any field-uniform Hessian bound of this action.  Usage: SPECTRE_DIR=/path/to/spectre python3 figures/analyse_spectre_hessian.py"""
import os, sys, csv
sys.path.insert(0, os.environ.get("SPECTRE_DIR", "../spectre"))
import numpy as np
import spectre_wilson as S
from scipy.optimize import minimize

def top(k):
    return np.linalg.eigvalsh(S.hessian(np.atleast_2d(np.array(k, float))))[0, -1]

N = 16
K = S.kgrid(N)
t = np.linalg.eigvalsh(S.hessian(K))[:, -1]
best = max((minimize(lambda k: -top(k), K[j], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=4000))
            for j in np.argsort(-t)[:6]), key=lambda r: -r.fun)
print(f"weights (W_TRI, W_MIX) = ({S.W_TRI}, {S.W_MIX}); grid N = {N}: max top eigenvalue {t.max():.6f}")
print(f"refined: {-best.fun:.12f} at k / pi = {np.round(best.x / np.pi, 4)}")
print("eigenvalues of M(k) at k = 1e-3 e_1, divided by k^2:", np.round(np.linalg.eigvalsh(S.hessian(np.array([[1e-3, 0, 0, 0]])))[0][:4] / 1e-6, 4))
with open("figures/spectre_hessian.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["W_TRI", "W_MIX", "grid_N", "grid_max_top", "refined_max_top", "k_over_pi"])
    w.writerow([S.W_TRI, S.W_MIX, N, t.max(), -best.fun, " ".join(f"{x:.4f}" for x in best.x / np.pi)])
