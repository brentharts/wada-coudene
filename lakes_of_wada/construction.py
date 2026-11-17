"""Wada-style constructions on a raster: canals dug into an island.

Brouwer/Yoneyama/Wada: start with an island and K lakes (here: K sectors of the surrounding ocean). At
stage k dig a canal from lake (k mod K) into the remaining land so that every
land point ends up within eps_k of that lake's water, while the land never
touches the water of another lake (a gap is kept) and shrinks to a connected
set.  With eps_k -> 0 every remaining land point is a limit of all K lakes, so
the land becomes their common boundary.  On a pixel grid we stop at
``min_eps``.

Output labels: 0 = remaining land, 1..K = lakes.
"""
import numpy as np
from scipy import ndimage as ndi

from .pipeline import Picture


def _disk_offsets(w):
    r = np.arange(-w, w + 1)
    dy, dx = np.meshgrid(r, r, indexing="ij")
    m = dx * dx + dy * dy <= w * w
    return dy[m], dx[m]


def wada_canals(n=384, lakes=3, stages=15, eps0=None, shrink=0.35, min_eps=1.0,
                island=0.45, layout="ponds", sep_frac=0.3, width_frac=1 / 20) -> Picture:
    """Raster Wada construction; returns a :class:`Picture` (no dynamics).

    Stage k digs lake (k mod K) from its water through the remaining land until
    every reachable land pixel is within about ``eps_k = eps0 * shrink**k`` of
    it.  Canals are geodesic (they wind through land corridors), have width
    ``eps_k * width_frac``, and stay ``eps_k * sep_frac`` away from every other
    lake.  The thin canals and fast shrink (0.35) are what keep the land
    corridors left by earlier stages wide enough for the next lake to pass; with
    a slower shrink, one lake gets sealed in.

    Resolution limit: the construction is exact only down to the scale at which
    every lake has been dug at least once (about eps0 * shrink**(K-1)), so the
    Wada profile (see ``wadatest``) saturates at a radius of roughly n/8 pixels
    for the defaults, and improves, in units of n, as n grows.  Stages stop when
    eps < ``min_eps`` pixels.
    """
    yy, xx = np.mgrid[0:n, 0:n]
    cy = cx = (n - 1) / 2
    R = island * n
    rad = np.hypot(yy - cy, xx - cx)
    label = np.zeros((n, n), dtype=int)
    # the surrounding ocean is split into `lakes` equal sectors, one per lake, so
    # every lake touches the island and has the same opportunity to dig
    sector = np.floor(((np.arctan2(yy - cy, xx - cx) + np.pi) / (2 * np.pi)) * lakes).astype(int)
    sector = np.clip(sector, 0, lakes - 1) + 1
    if layout == "sectors":
        label[rad > R] = sector[rad > R]
    else:                                     # classical: sea outside, ponds inside
        label[rad > R] = 1
        for k in range(2, lakes + 1):
            ang = 2 * np.pi * (k - 2) / max(lakes - 1, 1) + 0.3
            py, px = cy + 0.35 * R * np.sin(ang), cx + 0.35 * R * np.cos(ang)
            label[np.hypot(yy - py, xx - px) <= max(2.0, 0.04 * n)] = k
    eps = eps0 if eps0 is not None else n / 4
    done = 0
    for k in range(stages):
        if eps < min_eps:
            break
        w = int(eps * width_frac)
        _dig(label, k % lakes + 1, lakes, eps, w, gap=1, max_targets=10 ** 9, sep_frac=sep_frac)
        eps *= shrink
        done += 1
    return Picture(label.astype(float), label, np.zeros_like(label),
                   dict(kind="canal", n=n, lakes=lakes, stages=done, eps0=eps0,
                        shrink=shrink, min_eps=min_eps, island=island, layout=layout))


_NB8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def _geodesic(water, allowed):
    """BFS steps (8-connected) from ``water`` through ``allowed`` pixels; inf if unreachable."""
    G = np.full(water.shape, np.inf)
    G[water] = 0
    visited = water.copy()
    cur = water
    st = np.ones((3, 3), bool)
    k = 0
    while cur.any():
        k += 1
        new = ndi.binary_dilation(cur, st) & allowed & ~visited
        G[new] = k
        visited |= new
        cur = new
    return G


def _dig(label, j, lakes, eps, w, gap, max_targets=10 ** 9, sep_frac=0.25):
    """Dig a canal tree from lake j so every reachable land pixel is within ~eps of it.
    Canals wind through the land (geodesic paths) and stay ``gap`` pixels away
    from every other lake.  Digs at most ``max_targets`` canals; returns True
    when lake j has nothing left to do at this eps."""
    n = label.shape[0]
    water = label == j
    if not water.any():
        return True
    others = (label > 0) & (label != j)
    # keep a scale-proportional separation so corridors stay open for the other lakes
    forbidden = ndi.distance_transform_edt(~others) <= max(gap, eps * sep_frac) + 1.0
    allowed = (label == 0) & ~forbidden
    G = _geodesic(water, allowed)
    dist = ndi.distance_transform_edt(~water)
    skip = np.zeros(label.shape, bool)
    dy, dx = _disk_offsets(w)
    for _ in range(max_targets):
        cand = np.where(allowed & np.isfinite(G) & ~skip, dist, 0.0)
        t = np.unravel_index(np.argmax(cand), cand.shape)
        if cand[t] <= eps:
            return True
        # walk down the geodesic gradient from t to the water
        path = [t]
        cur = t
        while G[cur] > 1:
            best = None
            for oy, ox in _NB8:
                y, x = cur[0] + oy, cur[1] + ox
                if 0 <= y < n and 0 <= x < n and G[y, x] == G[cur] - 1:
                    best = (y, x)
                    break
            if best is None:
                break
            path.append(best)
            cur = best
        drop = int(eps / 2)                      # stay within eps/2 of the target
        path = path[drop:] if len(path) > drop else path[-1:]
        py = np.array([p[0] for p in path]); px = np.array([p[1] for p in path])
        by = (py[:, None] + dy[None, :]).ravel(); bx = (px[:, None] + dx[None, :]).ravel()
        ok = (by >= 0) & (by < n) & (bx >= 0) & (bx < n)
        by, bx = by[ok], bx[ok]
        sel = (label[by, bx] == 0) & ~forbidden[by, bx]
        by, bx = by[sel], bx[sel]
        label[by, bx] = j
        allowed[by, bx] = False
        if by.size:
            m = int(eps) + 2                     # local distance-field update
            y0, y1 = max(by.min() - m, 0), min(by.max() + m + 1, n)
            x0, x1 = max(bx.min() - m, 0), min(bx.max() + m + 1, n)
            patch = np.ones((y1 - y0, x1 - x0), bool)
            patch[by - y0, bx - x0] = False
            dist[y0:y1, x0:x1] = np.minimum(dist[y0:y1, x0:x1],
                                            ndi.distance_transform_edt(patch))
            dist[by, bx] = 0.0
            G[by, bx] = 0.0
        if dist[t] > eps:
            skip[t] = True
    return False
