"""Numerical tests of the Wada property on a label array of any dimension.

A point of a common boundary of K >= 3 domains has *every* one of the K
domains arbitrarily close.  Two complementary tests are provided:

* the **Wada profile** (:func:`wada_profile`): for each boundary pixel, which
  basins lie within r pixels?  Reports the fraction of boundary pixels that see
  all K basins as a function of r.  For a true Wada boundary this is ~1 once r
  exceeds the resolution; for an ordinary partition it stays far below 1.
  :func:`basin_profile` gives the same per basin, to expose *partially* Wada
  sets (some basins share a boundary, others do not).
* the **merging method** (:func:`merging_test`, Daza-Wagemakers-Sanjuán): merge
  all basins but one into a single colour; the set is Wada iff the resulting
  boundary is the same for every basin.  Compared by Hausdorff distance.

Label convention: ``>0`` = basins 1..K, ``0`` = undecided (counted as boundary),
``<0`` = invalid / outside the domain (ignored).  Arrays may have any number of
dimensions; distances are in pixels unless ``sampling`` (per-axis voxel size)
is given.
"""
import numpy as np
from scipy import ndimage as ndi


def basin_labels(label, min_fraction=1e-3):
    """Decided labels covering at least ``min_fraction`` of the valid points."""
    labs, counts = np.unique(label[label > 0], return_counts=True)
    nvalid = max(int(np.sum(label >= 0)), 1)
    return labs[counts >= min_fraction * nvalid]


def boundary_mask(label, width=1):
    """Valid points that are undecided or have a different valid label within
    ``width`` pixels (Chebyshev, any dimension).  Invalid points (label < 0)
    are never boundary and do not count as neighbours."""
    size = 2 * width + 1
    valid = label >= 0
    big = np.iinfo(np.int64).max
    lab = label.astype(np.int64)
    mx = ndi.maximum_filter(np.where(valid, lab, -1), size=size, mode="nearest")
    mn = ndi.minimum_filter(np.where(valid, lab, big), size=size, mode="nearest")
    return valid & ((mx != mn) | (label == 0))


def _distance_maps(label, labs, sampling=None):
    return {L: ndi.distance_transform_edt(label != L, sampling=sampling) for L in labs}


def wada_map(label, r, min_fraction=1e-3, sampling=None):
    """Fraction (0..1) of the K basins lying within ``r`` of each point."""
    labs = basin_labels(label, min_fraction)
    dist = _distance_maps(label, labs, sampling)
    cnt = sum((dist[L] <= r).astype(int) for L in labs)
    return cnt / max(len(labs), 1)


def wada_profile(label, radii=(1, 2, 4, 8), boundary=None, min_fraction=1e-3,
                 sampling=None):
    """``{"K", "boundary_pixels", "profile": {r: fraction of boundary points
    seeing all K basins within r}}``.  ``boundary`` defaults to
    :func:`boundary_mask`."""
    labs = basin_labels(label, min_fraction)
    dist = _distance_maps(label, labs, sampling)
    b = boundary_mask(label) if boundary is None else boundary
    out = {}
    for r in radii:
        cnt = sum((dist[L][b] <= r).astype(int) for L in labs)
        out[r] = float(np.mean(cnt == len(labs))) if b.any() else float("nan")
    return {"K": len(labs), "boundary_pixels": int(b.sum()), "profile": out}


def basin_profile(label, radii=(1, 2, 4, 8), min_fraction=1e-3, sampling=None):
    """Per-basin Wada profile: for each basin L, the fraction of *its* boundary
    points (points of L with a different valid label within 1 pixel) that see
    all K basins within r.  Values near 1 for every basin = Wada; a mixture of
    high and low = partially Wada."""
    labs = basin_labels(label, min_fraction)
    dist = _distance_maps(label, labs, sampling)
    bm = boundary_mask(label)
    out = {}
    for L in labs:
        b = bm & (label == L)
        row = {}
        for r in radii:
            cnt = sum((dist[M][b] <= r).astype(int) for M in labs)
            row[r] = float(np.mean(cnt == len(labs))) if b.any() else float("nan")
        out[int(L)] = row
    return out


# --- merging method -----------------------------------------------------------

def slim_boundary(binary, valid=None):
    """Points of a two-colour array that have a differently coloured neighbour
    (3-neighbourhood, any dimension).  Invalid points (``valid`` False) are never
    boundary and are ignored as neighbours."""
    b = binary.astype(np.int8)
    if valid is None:
        valid = np.ones(b.shape, bool)
    mx = ndi.maximum_filter(np.where(valid, b, -1), size=3, mode="nearest")
    mn = ndi.minimum_filter(np.where(valid, b, 2), size=3, mode="nearest")
    return valid & (mx != mn)


def _hausdorff(a, b, sampling=None):
    da = ndi.distance_transform_edt(~a, sampling=sampling)   # distance to a
    db = ndi.distance_transform_edt(~b, sampling=sampling)
    return float(max(db[a].max(), da[b].max()))


def merging_test(label, min_fraction=1e-3, sampling=None):
    """Merging method.  For each basin, merge everything else (all other basins,
    undecided points) into one colour and extract the slim boundary; return the
    pairwise Hausdorff distances between these boundaries.

    Wada  <=>  all the boundaries coincide, i.e. ``max_d`` is small compared with
    the array size (``max_d_normalized`` is ``max_d`` over the array diagonal) and
    comparable to ``min_d`` (``ratio = (max_d - min_d) / min_d``).  The decision
    is left to the caller: results are only reliable up to the resolution.
    Invalid points (label < 0) are excluded.
    """
    labs = basin_labels(label, min_fraction)
    valid = label >= 0
    bnds = {}
    for L in labs:
        bnds[int(L)] = slim_boundary(label == L, valid)
    pairs = {}
    labs_i = [int(L) for L in labs]
    for i, a in enumerate(labs_i):
        for b in labs_i[i + 1:]:
            if bnds[a].any() and bnds[b].any():
                pairs[(a, b)] = _hausdorff(bnds[a], bnds[b], sampling)
    ds = list(pairs.values())
    diag = float(np.linalg.norm(np.array(label.shape) * (1 if sampling is None else np.asarray(sampling))))
    if not ds:
        return {"K": len(labs_i), "pairs": pairs, "max_d": float("nan"),
                "min_d": float("nan"), "ratio": float("nan"),
                "max_d_normalized": float("nan"), "diagonal": diag}
    mx, mn = max(ds), min(ds)
    return {"K": len(labs_i), "pairs": pairs, "max_d": mx, "min_d": mn,
            "ratio": (mx - mn) / mn if mn > 0 else float("inf"),
            "max_d_normalized": mx / diag, "diagonal": diag}
