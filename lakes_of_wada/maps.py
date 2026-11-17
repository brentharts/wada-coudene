"""Torus maps: a hyperbolic toral automorphism plus localized perturbations.

Points are float arrays of shape ``(..., 2)`` in torus coordinates (period 1).
A map is any callable ``f(p) -> p'`` with the same shape; ``Anosov`` below is
the family used by Coudene, but anything with that signature plugs into the
classifiers and views.
"""
from dataclasses import dataclass, field
from typing import Callable, Sequence
import numpy as np

PHI = (1 + 5 ** 0.5) / 2
#: cat map; hyperbolic with eigenvalues phi^2, phi^-2
CAT = np.array([[2.0, 1.0], [1.0, 1.0]])
#: projection onto the unstable direction of CAT (symmetric: eigenvector (phi,1))
UNSTABLE_PROJ = np.array([[PHI ** 2, PHI], [PHI, 1.0]]) / (1 + PHI ** 2)


def wrap(p):
    """Reduce to the fundamental domain [-1/2, 1/2)^2."""
    return (p + 0.5) % 1.0 - 0.5


def torus_distance(a, b):
    """Euclidean distance on the flat torus (shape (..., 2) -> (...))."""
    return np.linalg.norm(wrap(a - b), axis=-1)


# --- bump profiles: t -> weight in [0, 1], compactly supported on |t| <= 1 ----

def bump_product(u):
    """(1-x^2)^2 (1-y^2)^2 on the unit square; u has shape (..., 2)."""
    t = np.clip(1 - u ** 2, 0, None) ** 2
    return t[..., 0] * t[..., 1]


def bump_radial(u):
    """(1-r^2)^2 on the unit disk."""
    r2 = np.sum(u ** 2, axis=-1)
    return np.clip(1 - r2, 0, None) ** 2


PROFILES = {"product": bump_product, "radial": bump_radial}


@dataclass
class Perturbation:
    """Subtract ``amplitude * profile(offset/spread) * (B @ offset)`` where
    ``offset`` is the displacement from the nearest point of the lattice
    ``(spacing * Z)^2``.  spacing=1 -> one centre (origin); spacing=0.5 -> the
    four 2-torsion points (0,0), (0,1/2), (1/2,0), (1/2,1/2).
    """
    amplitude: float
    spread: float
    spacing: float = 1.0
    profile: str = "product"
    # Optionally restrict to centres other than the origin class, e.g. the
    # plane picture perturbs the origin and the other 2-torsion points separately.
    exclude_origin: bool = False

    def offset(self, p):
        return p - self.spacing * np.floor(p / self.spacing + 0.5)

    def __call__(self, p):
        off = self.offset(p)
        w = PROFILES[self.profile](off / self.spread)
        if self.exclude_origin:
            # centres whose nearest lattice point is the origin class (integers)
            near_int = np.all(np.abs(p - np.round(p)) <= 0.25, axis=-1)
            w = np.where(near_int, 0.0, w)
        return self.amplitude * w[..., None] * (off @ UNSTABLE_PROJ.T)


@dataclass
class Anosov:
    """f(p) = wrap(A p - sum_i perturbation_i(p))."""
    perturbations: Sequence[Callable] = field(default_factory=list)
    matrix: np.ndarray = field(default_factory=lambda: CAT.copy())

    def __call__(self, p):
        p = wrap(np.asarray(p, dtype=float))   # make the map genuinely periodic
        out = p @ self.matrix.T
        for pert in self.perturbations:
            out = out - pert(p)
        return wrap(out)


# --- presets matching the original scripts ------------------------------------

def torus_map(amplitude=2.2, spread=0.23, islands=4):
    """wada-torus.py: product bump at 1 or 4 centres."""
    if islands not in (1, 4):
        raise ValueError("islands must be 1 or 4")
    spacing = 1.0 if islands == 1 else 0.5
    return Anosov([Perturbation(amplitude, spread, spacing, "product")])


def plane_map(amplitude1=2.4, spread1=0.5, amplitude2=2.3, spread2=0.25):
    """wada-plane.py: radial bump at the origin plus radial bumps at the other
    three 2-torsion points."""
    return Anosov([
        Perturbation(amplitude1, spread1, 1.0, "radial"),
        Perturbation(amplitude2, spread2, 0.5, "radial", exclude_origin=True),
    ])
