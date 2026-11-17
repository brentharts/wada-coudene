"""Inverse of the Weierstrass P function for the square (lemniscatic) lattice.

The original script used Sage's ``elliptic_f(arcsin(z**-0.5), -1)``.  Here we
use the Carlson symmetric form instead, which scipy supports for complex
arguments:  F(phi | m=-1) = sin(phi) * RF(cos^2 phi, 1 + sin^2 phi, 1)
with sin(phi) = z**-0.5.  No Sage required.
"""
import numpy as np
from scipy.special import elliprf

# Gamma(3/4) / (2 sqrt(pi) Gamma(5/4)): normalizes the period lattice to Z[i]
INV_PERIOD = 0.38137988175090659403

_TINY = 1e-300


def inv_p(z):
    """Return w with P(w) = z (principal branch), vectorized over ``z``.

    Real z in (0, 1) lies on a branch cut of RF; we resolve it as the limit
    from Im z -> 0+ (this matches mpmath / Sage's behaviour there).
    """
    z = np.array(z, dtype=complex)  # copy
    cut = (z.imag == 0) & (z.real > 0) & (z.real < 1)
    z[cut] += 1j * _TINY
    t = z ** -0.5
    t2 = t * t
    return INV_PERIOD * t * elliprf(1 - t2, 1 + t2, 1 + 0j)
