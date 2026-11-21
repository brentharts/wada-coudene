"""Rasterising the four window pieces Z_p(S) of the Spectre (S = O or G, parity p = 0 or 1) from their exact graph-directed system
(brentharts/spectre, boundary_certificate/spectre_window.py; maps dumped by figures/dump_spectre_window_maps.py into figures/spectre_window_maps.npz).

    Z_p(S) = U over (i, S') in MAPS[(S, p)] of  O_{i,p}^T ( GA_{i,p} + Z_{p+1}(S') / lam ),   lam = sqrt(4 + sqrt 15).

``points(S, p, depth)`` returns the images of the origin under all compositions of ``depth`` maps (every one is within rad * lam^-depth of
the attractor, with rad = 4.10034 the radius bound of the source) together with the index of the first map used (the child piece).
``label_image`` rasterises them by nearest point.  Nothing here is the window itself, only a finite approximation of it."""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

_DATA = None


def _data():
    global _DATA
    if _DATA is None:
        d = np.load(Path(__file__).resolve().parent.parent / "figures" / "spectre_window_maps.npz")
        _DATA = {k: d[k] for k in d.files}
    return _DATA


def _level(S, p, depth, cache):
    key = (S, p, depth)
    if key in cache:
        return cache[key]
    d = _data()
    lam = float(d["lam"])
    if depth == 0:
        res = (np.zeros((1, 2)), np.zeros(1, dtype=int))
    else:
        pts, kid = [], []
        for k, (i, s2) in enumerate(d[f"{S}{p}"]):
            sub, _ = _level("OG"[s2], 1 - p, depth - 1, cache)
            img = (d["GA"][i, p] + sub / lam) @ d["O"][i, p]        # O^T (GA + z / lam), as row vectors
            pts.append(img); kid.append(np.full(len(img), k))
        res = (np.concatenate(pts), np.concatenate(kid))
    cache[key] = res
    return res


def points(S, p, depth):
    """(n, 2) points and the child index (position in MAPS[(S, p)]) of each."""
    return _level(S, p, depth, {})


def count(S, p, depth):
    c = {("O", 0): 1, ("O", 1): 1, ("G", 0): 1, ("G", 1): 1}
    d = _data()
    for _ in range(depth):
        c = {(S, p): sum(c[("OG"[s2], 1 - p)] for _, s2 in d[f"{S}{p}"]) for S in "OG" for p in (0, 1)}
    return c[(S, p)]


def label_image(pts, lab, n, box=None, reach=1.5):
    """n x n image (index [i, j] = (x_i, y_j)) of the label of the nearest point, -1 where the nearest point is farther than ``reach`` pixels.
    Returns (image, box, pixel size)."""
    if box is None:
        lo, hi = pts.min(0), pts.max(0)
        c = 0.5 * (lo + hi); h = 0.5 * float((hi - lo).max()) * 1.02
        box = (c[0] - h, c[0] + h, c[1] - h, c[1] + h)
    x0, x1, y0, y1 = box
    px = (x1 - x0) / n
    xs = x0 + (np.arange(n) + 0.5) * px; ys = y0 + (np.arange(n) + 0.5) * px
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    tree = cKDTree(pts)
    dist, idx = tree.query(np.column_stack([X.ravel(), Y.ravel()]), distance_upper_bound=reach * px)
    out = np.full(X.size, -1, dtype=int)
    ok = np.isfinite(dist)
    out[ok] = lab[idx[ok]]
    return out.reshape(n, n), box, px
