"""Exit basins of a particle with quantised-inertia (MiHsC-type) dynamics.

Toy for the connection to McCulloch and Gine, Adv. Stud. Theor. Phys. 14 (2020) 1-8
(see QINERTIA.md).  Units: c = 1, GM of each attractor = 1, spacing of the
three attractors = 1.  A test particle moves in the field g(x) of three softened
point masses at the corners of an equilateral triangle, but its acceleration is
not g: it follows an inertia law m_i(a).

Laws (``law=``), with g = |F|/m the Newtonian acceleration and a0 = 2 c^2 / Theta:

* ``"newton"``  a = g.  Control.
* ``"mond"``    the paper's Eq. 14, m_i = m / (1 + a0/a), solved for a given g:
                a = (g + sqrt(g^2 + 4 g a0)) / 2.  Single valued and smooth, a -> sqrt(a0 g)
                for g << a0 (the deep-MoND limit).
* ``"qi1"``     the paper's Eq. 1, m_i = m (1 - a0/|a|), on its positive-mass branch
                |a| > a0: F = m_i a gives |a| = g + a0.  The branch |a| < a0 (negative m_i,
                a antiparallel to F) is NOT used; the paper does not treat it.  g = 0 would give
                a jump from 0 to a0, so g/sqrt(g^2 + soft^2) replaces the unit step.

The cosmic horizon is at R_H = c^2 / a0 = Theta / 2 (``a0 = 1 / R_H`` here): the Rindler horizon
distance c^2 / a equals R_H exactly at a = a0.  The exit rule is: the particle has reached R_H,
labelled by which of three 120-degree sectors (centred on the attractors) it is in.

The modified-inertia law a = nu(g) g is not derived from a Lagrangian, so the Newtonian
energy is NOT conserved when a0 > 0; only ``law="newton"`` conserves it (a test).
"""
import numpy as np

from .scattering import System, integrate_exits
from .pipeline import Picture

ANGLES = np.deg2rad([90.0, 210.0, 330.0])
SIDE = 1.0
CENTERS = (SIDE / np.sqrt(3.0)) * np.stack([np.cos(ANGLES), np.sin(ANGLES)], axis=1)


def field(x, y, soft=0.1):
    """Newtonian acceleration (gx, gy) from three unit masses, Plummer-softened."""
    gx = np.zeros_like(x)
    gy = np.zeros_like(y)
    for cx, cy in CENTERS:
        dx, dy = cx - x, cy - y
        d3 = (dx * dx + dy * dy + soft * soft) ** 1.5
        gx += dx / d3
        gy += dy / d3
    return gx, gy


def potential(x, y, soft=0.1):
    out = np.zeros_like(x)
    for cx, cy in CENTERS:
        out -= 1.0 / np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + soft * soft)
    return out


def inertia_factor(g, a0, law, tiny=1e-300, reg=1e-3):
    """nu(g) = a / g for the chosen law (so that a = nu * g_vector)."""
    if law == "newton" or a0 == 0.0:
        return np.ones_like(g)
    if law == "mond":
        return 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * a0 / np.maximum(g, tiny)))
    if law == "qi1":
        return 1.0 + a0 / np.sqrt(g * g + (reg * a0) ** 2)
    raise ValueError(f"unknown law {law!r}")


def sector_label(x, y):
    """1, 2, 3 for the sector centred on the attractor at 90, 210, 330 degrees."""
    th = np.arctan2(y, x)
    return 1 + np.floor(np.mod(th - np.deg2rad(30.0), 2 * np.pi) / (2 * np.pi / 3.0)).astype(int)


def system(R_H=20.0, law="mond", soft=0.1, a0=None) -> System:
    """The system; ``a0`` defaults to 1/R_H (c = 1), ``a0=0`` or ``law='newton'`` is the control."""
    a0 = 1.0 / R_H if a0 is None else a0

    def rhs(s):
        x, y, vx, vy = s[:, 0], s[:, 1], s[:, 2], s[:, 3]
        gx, gy = field(x, y, soft)
        nu = inertia_factor(np.hypot(gx, gy), a0, law)
        return np.stack([vx, vy, nu * gx, nu * gy], axis=1)

    def exit_label(s):
        out = np.hypot(s[:, 0], s[:, 1]) > R_H
        return np.where(out, sector_label(s[:, 0], s[:, 1]), 0)

    def energy(s):
        return 0.5 * (s[:, 2] ** 2 + s[:, 3] ** 2) + potential(s[:, 0], s[:, 1], soft)

    return System(f"quantised-inertia-{law}", rhs, exit_label, energy)


def launch(pos, speed, direction_deg=0.0, mode="tangential", energy=None, soft=0.1):
    """Velocities for starting positions ``pos`` (m, 2).  ``tangential``: speed * (-y, x)/r, which
    commutes with the 120-degree rotation of the three attractors, so the picture has an exact
    cyclic symmetry.  ``beam``: the same fixed velocity (speed, direction_deg) everywhere."""
    if energy is not None:
        # fixed launch energy: |v| = sqrt(2 (E - Phi)); starts with E < Phi are invalid (NaN speed)
        k = 2.0 * (energy - potential(pos[:, 0], pos[:, 1], soft))
        speed = np.where(k >= 0, np.sqrt(np.maximum(k, 0.0)), np.nan)
    if mode == "tangential":
        r = np.maximum(np.hypot(pos[:, 0], pos[:, 1]), 1e-12)
        vx, vy = -speed * pos[:, 1] / r, speed * pos[:, 0] / r
    elif mode == "beam":
        ph = np.deg2rad(direction_deg)
        vx, vy = np.full(len(pos), speed * np.cos(ph)), np.full(len(pos), speed * np.sin(ph))
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return np.column_stack([pos, vx, vy])


def beam_states(n, half_width=2.0, speed=2.0, direction_deg=0.0, mode="tangential", energy=None, soft=0.1):
    """n x n starting positions in a square (array order [ix, iy]) with velocities from ``launch``."""
    xs = np.linspace(-half_width, half_width, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    return launch(np.column_stack([X.ravel(), Y.ravel()]), speed, direction_deg, mode, energy, soft)


def qi_basins(n=200, R_H=20.0, law="mond", soft=0.1, a0=None, dt=0.02, t_max=None,
              half_width=2.0, speed=2.0, direction_deg=0.0, mode="tangential", energy=None,
              return_final=False):
    """Exit labels on an n x n grid of launch points (array indexed [ix, iy]); 0 = no exit by t_max,
    -1 = invalid start (``energy`` given and below the potential there).  With ``energy`` the launch speed
    is sqrt(2 (E - Phi)) so that every start has the same Newtonian energy."""
    t_max = 6.0 * R_H / speed if t_max is None else t_max
    s0 = beam_states(n, half_width, speed, direction_deg, mode, energy, soft)
    valid = np.all(np.isfinite(s0), axis=1)
    lab = np.full(len(s0), -1, dtype=int)
    fin = np.full_like(s0, np.nan)
    l, _, f = integrate_exits(system(R_H, law, soft, a0), s0[valid], dt, t_max)
    lab[valid], fin[valid] = l, f
    lab = lab.reshape(n, n)
    return (lab, fin.reshape(n, n, 4)) if return_final else lab


def qi_picture(**kw) -> Picture:
    lab = qi_basins(**kw)
    return Picture(lab.astype(float), lab, np.zeros_like(lab), dict(kind="qinertia", **kw))


def uncertainty_fraction(eps=(1e-4, 3e-4, 1e-3, 3e-3, 1e-2), n_samples=4000, seed=0,
                         half_width=2.0, dt=0.02, **kw):
    """f(eps): fraction of random launch points whose exit sector changes when the launch point is
    moved by eps in a random direction (pairs with an undecided end are dropped).  Returns
    {"alpha", "dimension", "f", "n_pairs", "undecided"}."""
    rng = np.random.default_rng(seed)
    speed = kw.get("speed", 2.0)
    energy = kw.get("energy")
    R_H = kw.get("R_H", 20.0)
    t_max = kw.get("t_max") or 6.0 * R_H / speed
    sysm = system(R_H, kw.get("law", "mond"), kw.get("soft", 0.1), kw.get("a0"))
    p = rng.uniform(-half_width, half_width, (n_samples, 2))
    mode, ddeg = kw.get("mode", "tangential"), kw.get("direction_deg", 0.0)

    soft = kw.get("soft", 0.1)

    def state(pos):
        return launch(pos, speed, ddeg, mode, energy, soft)

    s_base = state(p)
    p = p[np.all(np.isfinite(s_base), axis=1)]
    n_samples = len(p)
    base, _, _ = integrate_exits(sysm, state(p), dt, t_max)
    f, npairs = {}, {}
    for e in eps:
        th = rng.uniform(0, 2 * np.pi, n_samples)
        s1 = state(p + e * np.column_stack([np.cos(th), np.sin(th)]))
        ok1 = np.all(np.isfinite(s1), axis=1)
        lab1 = np.zeros(n_samples, dtype=int)
        lab1[ok1], _, _ = integrate_exits(sysm, s1[ok1], dt, t_max)
        good = (base > 0) & (lab1 > 0)
        f[e] = float(np.mean(base[good] != lab1[good])) if good.any() else float("nan")
        npairs[e] = int(good.sum())
    pts = [(np.log(e), np.log(v)) for e, v in f.items() if 0 < v < 0.5]
    alpha = float(np.polyfit(*np.array(pts).T, 1)[0]) if len(pts) >= 2 else float("nan")
    return {"alpha": alpha, "dimension": 2 - alpha, "f": f, "n_pairs": npairs,
            "undecided": float(np.mean(base == 0))}


# --- the repository's own chaotic scatterers with the inertia law ----------------------------
# The three fixed softened attractors above turned out regular (smooth basins even for law="newton"):
# a gravity well has no barrier to trap a particle at positive energy.  Henon-Heiles and its
# harmonic-polynomial relatives do, so the law is applied to them instead.

def hh_system(degree=3, radius=5.0, law="mond", a0=None) -> System:
    """Henon-Heiles (degree 3) or V = r^2/2 + Im(z^n)/n with acceleration nu(g) * g, g the Newtonian
    acceleration.  ``a0`` defaults to 1/radius (c = 1: the exit circle is the 'cosmic horizon' c^2/a0).
    The Newtonian energy is conserved only for law='newton' or a0 = 0."""
    from .scattering import harmonic_polynomial
    base = harmonic_polynomial(degree, radius)
    a0 = 1.0 / radius if a0 is None else a0

    def rhs(s):
        d = base.rhs(s)
        gx, gy = d[:, 2], d[:, 3]
        nu = inertia_factor(np.hypot(gx, gy), a0, law)
        return np.stack([d[:, 0], d[:, 1], nu * gx, nu * gy], axis=1)

    return System(f"hh{degree}-qi-{law}", rhs, base.exit_label, base.energy)


def hh_labels(n=200, energy=0.25, degree=3, radius=5.0, law="mond", a0=None, dt=0.05, t_max=200.0):
    """Exit labels on the x = 0, p_x >= 0 section at launch energy ``energy`` (-1 = not allowed, 0 = no
    exit by t_max); same layout as ``scattering.harmonic_section``."""
    from .scattering import harmonic_section
    y0, valid = harmonic_section(degree, n, energy)
    lab = np.full(valid.shape, -1, dtype=int)
    lab[valid], _, _ = integrate_exits(hh_system(degree, radius, law, a0), y0[valid], dt, t_max)
    return lab


def hh_uncertainty(energy=0.25, degree=3, radius=5.0, law="mond", a0=None, eps=(1e-4, 3e-4, 1e-3, 3e-3, 1e-2),
                   n_samples=6000, dt=0.05, t_max=200.0, seed=0):
    """As ``scattering.uncertainty_fraction`` (y shifted by eps, p_x recomputed from the launch energy), for hh_system."""
    from .scattering import harmonic_window, _harm_potential
    rng = np.random.default_rng(seed)
    yr, pr = harmonic_window(degree, energy)
    y, py = rng.uniform(*yr, n_samples), rng.uniform(*pr, n_samples)
    sysm = hh_system(degree, radius, law, a0)

    def st(yv):
        px2 = 2 * (energy - _harm_potential(degree, 0.0, yv)) - py ** 2
        return np.stack([np.zeros_like(yv), yv, np.sqrt(np.maximum(px2, 0.0)), py], axis=1), px2 >= 0

    s0, ok0 = st(y)
    base = np.zeros(n_samples, dtype=int)
    base[ok0], _, _ = integrate_exits(sysm, s0[ok0], dt, t_max)
    f, npairs = {}, {}
    for e in eps:
        s1, ok1 = st(y + e)
        ok = ok0 & ok1
        lab1 = np.zeros(n_samples, dtype=int)
        lab1[ok], _, _ = integrate_exits(sysm, s1[ok], dt, t_max)
        good = ok & (base > 0) & (lab1 > 0)
        f[e] = float(np.mean(base[good] != lab1[good])) if good.any() else float("nan")
        npairs[e] = int(good.sum())
    pts = [(np.log(e), np.log(v)) for e, v in f.items() if 0 < v < 0.5]
    alpha = float(np.polyfit(*np.array(pts).T, 1)[0]) if len(pts) >= 2 else float("nan")
    return {"alpha": alpha, "dimension": 2 - alpha, "f": f, "n_pairs": npairs}


# --- Unruh-type noise ---------------------------------------------------------------------
def hh_noisy_labels(starts, law="mond", a0=None, D0=1e-4, mode="unruh", degree=3, radius=5.0, dt=0.05,
                    t_max=200.0, n_real=10, seed=0):
    """Exit-label probabilities with velocity noise added after every RK4 step.

    ``mode="unruh"``: velocity diffusion rate D_v = D0 * |a| (|a| the acceleration the law gives), the
    scaling of an Unruh bath, k_B T = hbar |a| / (2 pi c) (the paper's Eq. 9), with D0 absorbing hbar / (2 pi c m)
    and a coupling rate.  ``mode="const"``: D_v = D0 * 0.3 (0.3 is a typical |a| on the section), the repository's
    constant-temperature Langevin without friction.  D0 = 0 is the deterministic map.  Returns an array
    (len(starts), 4) of fractions of the ``n_real`` realisations that ended with label 0 (no exit), 1, 2, 3.
    The noise is a toy, with no friction and no derivation of D0 from hbar, c and m (see QINERTIA.md)."""
    rng = np.random.default_rng(seed)
    sysm = hh_system(degree, radius, law, a0)
    s0 = np.repeat(np.asarray(starts, float), n_real, axis=0)
    y, m = s0.copy(), len(s0)
    lab = np.zeros(m, dtype=int)
    active = np.arange(m)
    f = sysm.rhs
    for _ in range(int(round(t_max / dt))):
        if active.size == 0:
            break
        ya = y[active]
        k1 = f(ya); k2 = f(ya + 0.5 * dt * k1); k3 = f(ya + 0.5 * dt * k2); k4 = f(ya + dt * k3)
        ya = ya + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        if D0 > 0:
            acc = np.hypot(k1[:, 2], k1[:, 3])
            dv = (D0 * acc if mode == "unruh" else np.full(len(acc), D0 * 0.3))
            ya[:, 2:] += np.sqrt(2.0 * dv * dt)[:, None] * rng.standard_normal((len(ya), 2))
        y[active] = ya
        l = sysm.exit_label(ya)
        lab[active[l > 0]] = l[l > 0]
        active = active[l == 0]
    out = np.zeros((len(s0) // n_real, 4))
    lab = lab.reshape(-1, n_real)
    for k in range(4):
        out[:, k] = np.mean(lab == k, axis=1)
    return out
