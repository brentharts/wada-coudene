"""Dynamics on the complex plane: entire functions and Newton maps.

A map is any vectorized callable ``f(z) -> z'`` on complex arrays.  The
classifier labels each starting point by what its orbit does:

* label 1: escapes (``|z| > escape_radius`` or overflows)
* label 2, 3, ...: converges to an attracting cycle (numbered in a fixed order)
* label 0: undecided within the iteration budget (typically Julia set / chaos)
"""
import numpy as np
from numpy.polynomial import polynomial as P

from .classify import Result


# --- map families -------------------------------------------------------------

def polynomial(coeffs):
    """z -> sum coeffs[k] z^k  (coeffs lowest degree first)."""
    c = np.asarray(coeffs, dtype=complex)
    return lambda z: P.polyval(z, c)


def newton(coeffs=None, degree=3):
    """Newton map z - p/p' of p (default z^degree - 1).  The three basins of
    z^3 - 1 are the classic example of a common boundary (Wada property)."""
    if coeffs is None:
        coeffs = np.r_[-1, np.zeros(degree - 1), 1]
    c = np.asarray(coeffs, dtype=complex)
    dc = P.polyder(c)

    def f(z):
        with np.errstate(all="ignore"):
            return z - P.polyval(z, c) / P.polyval(z, dc)
    return f


def exponential(lam=0.3):
    """z -> lam * exp(z)."""
    def f(z):
        with np.errstate(all="ignore"):
            return lam * np.exp(z)
    return f


def sine(lam=1.0):
    """z -> lam * sin(z)."""
    def f(z):
        with np.errstate(all="ignore"):
            return lam * np.sin(z)
    return f


FAMILIES = {"newton": newton, "poly": polynomial, "exp": exponential, "sin": sine}


# --- metric -------------------------------------------------------------------

def chordal(a, b):
    """Chordal (Riemann sphere) distance; finite even when an orbit escapes."""
    with np.errstate(all="ignore"):
        d = 2 * np.abs(a - b) / np.sqrt((1 + np.abs(a) ** 2) * (1 + np.abs(b) ** 2))
    return np.where(np.isfinite(d), d, 0.0)


# --- classifier ---------------------------------------------------------------

def classify_attractors(f, z0, max_iter=200, escape_radius=1e4, tol=1e-9,
                        max_period=3, key_digits=5) -> Result:
    """Label orbits of ``f`` by escape / attracting cycle (see module doc)."""
    z = np.array(z0, dtype=complex)
    shape = z.shape
    z = z.ravel()
    n = z.size
    label = np.zeros(n, dtype=int)
    steps = np.zeros(n, dtype=int)
    cyc = {}                      # index -> representative cycle point
    active = np.arange(n)
    for k in range(max_iter):
        if active.size == 0:
            break
        zn = f(z[active])
        zn = np.where(np.isfinite(zn), zn, np.inf)
        z[active] = zn
        steps[active] = k + 1
        esc = ~(np.abs(zn) < escape_radius)
        label[active[esc]] = 1
        alive = active[~esc]
        # cycle detection every few steps, once the orbit has settled a bit
        if k >= 4 and k % 3 == 0 and alive.size:
            w = z[alive].copy()
            hit = np.zeros(alive.size, bool)
            rep = np.zeros(alive.size, complex)
            cur = w.copy()
            for p in range(1, max_period + 1):
                cur = f(cur)
                close = (np.abs(cur - w) < tol) & ~hit
                if close.any():
                    # representative: lexicographically smallest of the cycle
                    pts = [w[close]]
                    c2 = w[close]
                    for _ in range(p - 1):
                        c2 = f(c2)
                        pts.append(c2)
                    st = np.stack(pts)                       # (p, m)
                    key = np.round(st.real, key_digits) * 1e6 + np.round(st.imag, key_digits)
                    rep[close] = st[key.argmin(axis=0), np.arange(st.shape[1])]
                    hit |= close
            for idx, r in zip(alive[hit], rep[hit]):
                cyc[idx] = r
            alive = alive[~hit]
        active = alive
    if cyc:
        idxs = np.array(list(cyc.keys()))
        reps = np.array([cyc[i] for i in idxs])
        keys = np.stack([np.round(reps.real, key_digits - 1),
                         np.round(reps.imag, key_digits - 1)], axis=1)
        order = np.lexsort((keys[:, 1], keys[:, 0]))
        uniq, inv = np.unique(keys[order], axis=0, return_inverse=True)
        inv = inv.ravel()
        lab = np.empty(len(idxs), int)
        lab[order] = inv + 2
        label[idxs] = lab
    return Result(label.reshape(shape), steps.reshape(shape), z.reshape(shape))
