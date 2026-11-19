"""Dump the exact graph-directed system of the Spectre window (brentharts/spectre, boundary_certificate/spectre_window.py) as floats, so that
figures/analyse_spectre_window.py does not need the spectre repository (importing it takes about 70 s).
Z_p(S) = U over (i, S') in MAPS[(S, p)] of  O_{i,p}^T ( GA_{i,p} + Z_{p+1}(S') / lam ),  O orthogonal 2 x 2, GA a 2-vector, lam = sqrt(4 + sqrt 15).
Usage: SPECTRE_DIR=/path/to/spectre python3 figures/dump_spectre_window_maps.py     -> figures/spectre_window_maps.npz"""
import os, sys
sys.path.insert(0, os.path.join(os.environ.get("SPECTRE_DIR", "../spectre"), "boundary_certificate"))
import numpy as np
import spectre_window as W

O = np.zeros((8, 2, 2, 2)); GA = np.zeros((8, 2, 2))
for i in range(8):
    for p in (0, 1):
        O[i, p] = [[float(W.O[(i, p)][r, c]) for c in range(2)] for r in range(2)]
        GA[i, p] = [float(W.GA[(i, p)][r]) for r in range(2)]
maps = {}
for (S, p), lst in W.MAPS.items():
    maps[f"{S}{p}"] = np.array([[i, 0 if S2 == "O" else 1] for i, S2 in lst])
np.savez("figures/spectre_window_maps.npz", O=O, GA=GA, lam=float(W.LAM), rad=float(W.RAD), **maps)
print("saved; pieces per set:", {k: len(v) for k, v in maps.items()}, "lam", float(W.LAM), "radius bound", float(W.RAD))
