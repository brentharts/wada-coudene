"""Maverick-point detector (after Martí-Pete, Rempe & Waterman, arXiv:2108.10256).

For a wandering domain U, a boundary point z is *maverick* if some subsequence
f^{n_k}(z) converges to a point w that is not a limit of f^{n_k}(U); equivalently,
for w in U the limsup_n dist(f^n(z), f^n(w)) is positive.

Numerically, for every boundary pixel p and every basin L we pick the nearest
*interior* pixel q_L of basin L (at least ``margin`` pixels deep, with an orbit that
settled before the comparison window) and track

    div_L(p) = max over the last ``tail`` of ``n_iter`` iterates of
               dist(f^n(p), f^n(q_L)),

a finite-window proxy for the limsup.  div_L(p) ~ 0 means p shadows basin L;
``score(p) = min_L div_L(p)`` is large only if p tracks *none* of the adjacent
basins, i.e. p looks maverick with respect to all of them.

Resolution: if a basin has no pixel ``margin`` deep (it is interleaved with the
others at pixel scale) any settled pixel of that basin is used as the reference.

Caveats: (1) the maps are chaotic, so double-precision orbits are only
pseudo-orbits; keep ``n_iter`` moderate (tens of steps) and treat the result as
a heat map, not a proof.  (2) with a finite window the quantity is only a
proxy for the limsup, and a finite grid can only approximate "boundary".
"""
import numpy as np
from scipy import ndimage as ndi

from .wadatest import boundary_mask


def _orbit_tail(f, p, n_iter, tail):
    p = np.array(p)
    out = []
    for k in range(n_iter):
        p = f(p)
        if k >= n_iter - tail:
            out.append(p.copy())
    return np.stack(out, axis=1)               # (m, tail, ...)


def orbit_divergence(pic, basins=None, n_iter=60, tail=20, margin=2, region=None):
    """Return ``field`` of shape ``(n0, n1, B)`` and the list of basin labels.

    ``pic`` must carry ``f``, ``points`` and ``distance`` (torus / plane /
    entire pictures do).  ``region`` is a boolean mask of pixels to evaluate
    (default: the boundary mask); other pixels are NaN.
    """
    if pic.f is None:
        raise ValueError("this picture has no dynamics attached")
    label, pts = pic.label, pic.points
    n0, n1 = label.shape
    if basins is None:
        basins = [int(b) for b in np.unique(label[label > 0])]
    mask = boundary_mask(label, width=margin) if region is None else region
    sel = np.flatnonzero(mask.ravel())
    flat_pts = pts.reshape((n0 * n1,) + pts.shape[2:])

    partner = {}
    for L in basins:
        # the reference orbit must have settled before the comparison window starts
        settled = (label == L) & (pic.steps <= n_iter - tail)
        interior = ndi.binary_erosion(label == L, iterations=margin) & settled
        if not interior.any():
            # basin too finely interleaved at this resolution for a deep interior
            # (itself a Wada-like sign): fall back to any settled pixel
            interior = settled
        if not interior.any():
            continue
        _, idx = ndi.distance_transform_edt(~interior, return_indices=True)
        flat = np.ravel_multi_index((idx[0].ravel(), idx[1].ravel()), (n0, n1))
        partner[L] = flat[sel]
    used = sorted(partner)

    need = np.unique(np.concatenate([sel] + [partner[L] for L in used])) if used else sel
    tails = _orbit_tail(pic.f, flat_pts[need], n_iter, tail)
    pos = {g: i for i, g in enumerate(need)}
    row = np.vectorize(pos.get)(sel)
    field = np.full((n0 * n1, len(basins)), np.nan)
    for b, L in enumerate(basins):
        if L not in partner:
            continue
        prow = np.vectorize(pos.get)(partner[L])
        a, c = tails[row], tails[prow]
        d = np.stack([pic.distance(a[:, k], c[:, k]) for k in range(tail)], axis=1)
        field[sel, b] = d.max(axis=1)
    return field.reshape(n0, n1, len(basins)), basins


def maverick_score(field):
    """min over basins of the divergence (NaN where nothing was evaluated)."""
    with np.errstate(all="ignore"):
        s = np.nanmin(np.where(np.isnan(field), np.inf, field), axis=-1)
    s[np.isinf(s)] = np.nan
    return s


def maverick_mask(score, delta=0.1):
    """Candidate maverick pixels: score above ``delta`` in the metric's units."""
    return np.nan_to_num(score, nan=0.0) > delta
