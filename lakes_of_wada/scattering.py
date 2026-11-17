"""Continuous-time open systems: exit basins of Hamiltonian / ODE scattering.

A :class:`System` bundles an autonomous vector field ``rhs(state)`` (vectorized
over rows of an ``(m, d)`` array) with an ``exit_label(state)`` rule returning
0 while a trajectory is still inside and k >= 1 once it has left through exit
k.  :func:`integrate_exits` advances many trajectories at once with RK4 and
freezes each as soon as it exits, giving the exit basin of every initial
condition.  The exit basins of open Hamiltonian systems (Henon-Heiles,
binary black holes, ...) are the classic physical Wada examples.

Accuracy note: fixed-step RK4 is used for speed.  Trajectories in these systems
are chaotic, so individual pixels near basin boundaries are step-size dependent
(as with any numerical integrator); the large-scale basin structure is not.
"""
from dataclasses import dataclass
from typing import Callable, Optional
import numpy as np

from .pipeline import Picture


@dataclass
class System:
    name: str
    rhs: Callable           # (m, d) -> (m, d)
    exit_label: Callable    # (m, d) -> (m,) int, 0 = still inside
    energy: Optional[Callable] = None   # (m, d) -> (m,), for diagnostics


def integrate_exits(system: System, y0, dt=0.05, t_max=60.0):
    """Integrate each row of ``y0`` until it exits or ``t_max`` is reached.

    Returns ``(label, steps, final)``: exit index (0 = never exited by t_max),
    number of RK4 steps taken, and the final state.
    """
    y = np.array(y0, dtype=float)
    m = len(y)
    label = np.zeros(m, dtype=int)
    steps = np.zeros(m, dtype=int)
    f = system.rhs
    active = np.arange(m)
    for k in range(int(round(t_max / dt))):
        if active.size == 0:
            break
        ya = y[active]
        k1 = f(ya)
        k2 = f(ya + 0.5 * dt * k1)
        k3 = f(ya + 0.5 * dt * k2)
        k4 = f(ya + dt * k3)
        ya = ya + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        y[active] = ya
        steps[active] = k + 1
        finite = np.all(np.isfinite(ya), axis=1)
        lab = np.where(finite, system.exit_label(np.where(finite[:, None], ya, 0.0)), 0)
        done = (lab > 0) | ~finite          # blown-up trajectories stay undecided (label 0)
        label[active[lab > 0]] = lab[lab > 0]
        active = active[~done]
    return label, steps, y


# --- Henon-Heiles ---------------------------------------------------------------
# H = (px^2 + py^2)/2 + (x^2 + y^2)/2 + x^2 y - y^3/3,  state = (x, y, px, py).
# For E > 1/6 the three channels around theta = 90, 210, 330 degrees are open.

HH_ESCAPE_ENERGY = 1.0 / 6.0


def _hh_potential(x, y):
    return 0.5 * (x * x + y * y) + x * x * y - y ** 3 / 3.0


def henon_heiles(radius=5.0) -> System:
    """Henon-Heiles with exit when r > ``radius`` (at which point the potential
    is far below any reachable energy, so the trajectory cannot return)."""
    def rhs(s):
        x, y, px, py = s[:, 0], s[:, 1], s[:, 2], s[:, 3]
        return np.stack([px, py, -(x + 2 * x * y), -(y + x * x - y * y)], axis=1)

    def exit_label(s):
        x, y = s[:, 0], s[:, 1]
        out = x * x + y * y > radius * radius
        theta = np.arctan2(y, x)                       # channels centred at 90, 210, 330 deg
        sector = np.floor(((theta - np.pi / 6) % (2 * np.pi)) / (2 * np.pi / 3)).astype(int) + 1
        return np.where(out, np.clip(sector, 1, 3), 0)

    def energy(s):
        return 0.5 * (s[:, 2] ** 2 + s[:, 3] ** 2) + _hh_potential(s[:, 0], s[:, 1])

    return System("henon-heiles", rhs, exit_label, energy)


def henon_heiles_section(n, energy, y_range=(-0.6, 1.0), py_range=(-0.75, 0.75)):
    """Initial conditions on the surface x = 0, px >= 0 at the given energy.

    Returns ``(y0, valid)`` with ``y0`` of shape ``(n, n, 4)`` indexed
    ``[row, col]``: columns run over y, rows over p_y (top row = largest p_y).
    ``valid`` marks points compatible with the energy.
    """
    y = np.linspace(*y_range, n)
    py = np.linspace(*py_range, n)[::-1]
    Y, PY = np.meshgrid(y, py)
    px2 = 2 * (energy - _hh_potential(0.0, Y)) - PY ** 2
    valid = px2 >= 0
    s = np.stack([np.zeros_like(Y), Y, np.sqrt(np.where(valid, px2, 0.0)), PY], axis=-1)
    return s, valid


def henon_heiles_volume(n, energy, x_range=(-1.2, 1.2), y_range=(-0.7, 1.0),
                        py_range=(-0.75, 0.75)):
    """The full 3D energy surface (half with px >= 0), parametrised by
    ``(x, p_y, y)``: array shape ``(n, n, n, 4)`` indexed ``[ix, row(py), col(y)]``,
    so ``volume[n // 2]`` is the x = 0 section above."""
    x = np.linspace(*x_range, n)
    y = np.linspace(*y_range, n)
    py = np.linspace(*py_range, n)[::-1]
    X, PY, Y = np.meshgrid(x, py, y, indexing="ij")
    px2 = 2 * (energy - _hh_potential(X, Y)) - PY ** 2
    valid = px2 >= 0
    s = np.stack([X, Y, np.sqrt(np.where(valid, px2, 0.0)), PY], axis=-1)
    return s, valid


def energy_window(energy, y_max=1.0):
    """A natural (y_range, py_range) window for the x = 0 section at this energy:
    from the left turning point of V(0, y) to ``y_max`` (the saddle), and
    |p_y| <= sqrt(2 E)."""
    from scipy.optimize import brentq
    v = lambda y: _hh_potential(0.0, y)
    y_min = brentq(lambda y: v(y) - energy, -10.0, 0.0)
    return (float(y_min), float(y_max)), (-float(np.sqrt(2 * energy)), float(np.sqrt(2 * energy)))


def uncertainty_fraction(energy, eps=(1e-4, 3e-4, 1e-3, 3e-3, 1e-2), n_samples=20000,
                         dt=0.05, t_max=200.0, radius=5.0, seed=0, **ranges):
    """Uncertainty exponent of the Henon-Heiles exit basins on the x = 0 section.

    Draw ``n_samples`` random initial conditions (y, p_y) at this energy, perturb
    y by ``eps`` (p_x is recomputed so the energy is unchanged), integrate both and
    record how often the two exits differ.  f(eps) ~ eps**alpha: alpha = 1 for a
    smooth boundary, alpha < 1 for a fractal one (boundary dimension 2 - alpha on
    the section).  Pairs in which either trajectory had not exited by ``t_max``
    are dropped.  Returns ``{"alpha", "dimension", "f", "n_pairs"}``.
    """
    rng = np.random.default_rng(seed)
    yr, pr = (ranges.get("y_range"), ranges.get("py_range"))
    if yr is None or pr is None:
        yr, pr = energy_window(energy)
    y = rng.uniform(*yr, n_samples)
    py = rng.uniform(*pr, n_samples)

    def state(yv):
        px2 = 2 * (energy - _hh_potential(0.0, yv)) - py ** 2
        ok = px2 >= 0
        return np.stack([np.zeros_like(yv), yv, np.sqrt(np.where(ok, px2, 0.0)), py], axis=1), ok

    sysm = henon_heiles(radius)
    s0, ok0 = state(y)
    lab0, _, _ = integrate_exits(sysm, s0[ok0], dt, t_max)
    base = np.zeros(n_samples, dtype=int)
    base[ok0] = lab0
    f, npairs = {}, {}
    for e in eps:
        s1, ok1 = state(y + e)
        ok = ok0 & ok1
        lab1, _, _ = integrate_exits(sysm, s1[ok], dt, t_max)
        a, b = base[ok], lab1
        good = (a > 0) & (b > 0)
        f[e] = float(np.mean(a[good] != b[good])) if good.any() else float("nan")
        npairs[e] = int(good.sum())
    pts = [(np.log(e), np.log(v)) for e, v in f.items() if v > 0 and v < 0.5]
    alpha = float(np.polyfit(*np.array(pts).T, 1)[0]) if len(pts) >= 2 else float("nan")
    return {"alpha": alpha, "dimension": 2 - alpha, "f": f, "n_pairs": npairs}


def scattering_picture(system="henon-heiles", n=300, energy=0.25, dim=2, dt=0.05,
                       t_max=60.0, radius=5.0, **ranges) -> Picture:
    """Exit-basin picture.  ``dim=2``: the x = 0 section in (y, p_y);
    ``dim=3``: the whole energy surface.  ``label`` is 1..3 for the exit used,
    0 if the trajectory had not exited by ``t_max``, and -1 for points not
    allowed at this energy."""
    if system != "henon-heiles":
        raise ValueError(f"unknown system {system!r}")
    if energy <= HH_ESCAPE_ENERGY:
        raise ValueError("Henon-Heiles is closed (no exits) for E <= 1/6")
    sysm = henon_heiles(radius)
    if dim == 2:
        y0, valid = henon_heiles_section(n, energy, **ranges)
    elif dim == 3:
        y0, valid = henon_heiles_volume(n, energy, **ranges)
    else:
        raise ValueError("dim must be 2 or 3")
    label = np.full(valid.shape, -1, dtype=int)
    steps = np.zeros(valid.shape, dtype=int)
    lab, st, _ = integrate_exits(sysm, y0[valid], dt, t_max)
    label[valid] = lab
    steps[valid] = st
    return Picture(label.astype(float), label, steps,
                   dict(kind="scatter", system=system, n=n, energy=energy, dim=dim,
                        dt=dt, t_max=t_max, radius=radius))
