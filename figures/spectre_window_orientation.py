"""Internal-space image of the Spectre tile centres (brentharts/spectre, lift4d.py), coloured by tile orientation: is the orientation a function of
the internal position, i.e. does it partition the window?  Usage: SPECTRE_DIR=/path/to/spectre python3 figures/spectre_window_orientation.py [levels]"""
import os, sys
sys.path.insert(0, os.environ.get("SPECTRE_DIR", "../spectre"))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spectre_patch as sp, lift4d as L

levels = int(sys.argv[1]) if len(sys.argv) > 1 else 5
tiles, lifted, _ = L.lift_patch(levels)
P = sp.SpectrePatch(levels)
S = L.bounded_plane(4)[0]
X = np.array([np.mean(np.array(t, float), axis=0) for t in tiles]) @ S.T
k = np.round(np.degrees(np.arctan2(P.placements[:, 1, 0], P.placements[:, 0, 0])) / 30.0).astype(int) % 12
# overlap measure at several cell sizes h: among cells with at least 5 tile centres, the fraction holding more than one orientation class
# (a partition with a boundary of dimension about 1.1 would give a fraction that falls quickly with h; below h = 0.05 there are too few points)
print(f"levels {levels}: {len(X)} tiles, orientation counts {np.bincount(k, minlength=12).tolist()}, Mystic {int(P.shape.sum())}")
for h in (0.2, 0.1, 0.05, 0.025):
    cell = np.floor(X / h).astype(int)
    keys = cell[:, 0] * 100003 + cell[:, 1]
    order = np.argsort(keys)
    groups = np.split(k[order], np.flatnonzero(np.diff(keys[order])) + 1)
    sel = [g for g in groups if len(g) >= 5]
    n = np.array([len(set(g.tolist())) for g in sel])
    print(f"h = {h}: {len(sel)} cells, mixed {np.mean(n > 1):.3f}, mean classes per cell {n.mean():.2f}, max {n.max()}")
fig, ax = plt.subplots(1, 2, figsize=(14, 7))
for a, (lab, t) in zip(ax, ((k, "orientation (12 classes)"), (k % 6, "orientation mod 6"))):
    a.scatter(X[:, 0], X[:, 1], c=lab, s=1.2, cmap="tab20"); a.set_title(t); a.set_aspect("equal")
plt.tight_layout(); plt.savefig("figures/spectre_window_orientation.png", dpi=70)
