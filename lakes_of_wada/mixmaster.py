"""Mixmaster (vacuum Bianchi IX) as an exit system: which axis collapses first?

The vacuum Bianchi IX universe has spatial slices that are three-spheres with three scale factors a, b, c
(a Berger sphere when two are equal, the unequal axis being the Hopf fibre).  With x_i = ln(a_i^2) and the BKL
time tau (dt = abc dtau), the field equations are, with q_i = exp(x_i) = a_i^2,

    x_i'' = (q_j - q_k)^2 - q_i^2        (i, j, k cyclic)                        [Cornish & Levin, eq. (1)]

and the Hamiltonian constraint is the conserved quantity

    I = sum_{i<j} x_i' x_j' - V = 0 ,    V = sum q_i^2 - 2 sum_{i<j} q_i q_j   [V is their eq. (2); the constraint was
                                                                                  found numerically, see the tests]

Between bounces V ~ 0 and the motion is a Kasner epoch: x_i' = S p_i with S = sum x_i' and sum p_i = sum p_i^2 = 1.
The Mixmaster dynamics never leaves; it bounces forever (the BKL map).  Cornish and Levin turned it into an *exit
system* by cutting pockets in the three corners of the potential in the anisotropy plane
(beta_+, beta_-) = ((1/6)(x_2 + x_3 - 2 x_1), (1/(2 sqrt 3))(x_2 - x_3)) at Omega = (1/6) sum x_i = 0, and found fractal
basin boundaries between the three outcomes ([gr-qc/9709037]).  **Their pocket geometry is not in the text available to
us, so the exit rule used here is our own**, defined in Kasner terms: a trajectory exits through axis i at the first
step with S < 0 (collapsing), |V| < kappa * sum x_j'^2 (curvature negligible: a Kasner epoch) and p_i > p_exit (axis i
carries almost all of the expansion rate, a Kasner epoch with large LK parameter u).  p_i -> 1 is the *pancake* limit
a_i ~ t, a_j, a_k ~ const: a Berger sphere whose collapsing axis i is the Hopf fibre.  So the three exits are the three
pancakes, one for each axis of SU(2).

Initial conditions.  A two-parameter family, the moment of maximum expansion: Omega = 0, S = 0 (volume stationary),
anisotropic velocity of magnitude fixed by the constraint, |beta'|^2 = -V/3 (needs V < 0), pointing radially outward from
the isotropic point.  The family is invariant under the cyclic permutation of axes (rotation by 120 degrees).

Accuracy: RK4 with a step h0 / (1 + max|x'|) chosen per trajectory.  Chaotic: individual pixels near basin boundaries
depend on the step and on the exit thresholds (tests report how many); the large-scale structure does not.
"""
import numpy as np

from .pipeline import Picture

_SQ3 = np.sqrt(3.0)


def potential(x):
    """V = sum q_i^2 - 2 sum_{i<j} q_i q_j with q = exp(x); x has shape (..., 3)."""
    q = np.exp(np.asarray(x, dtype=float))
    return np.sum(q ** 2, axis=-1) - 2.0 * (q[..., 0] * q[..., 1] + q[..., 1] * q[..., 2] + q[..., 0] * q[..., 2])


def force(x):
    q = np.exp(x)
    return np.stack([(q[..., (i + 1) % 3] - q[..., (i + 2) % 3]) ** 2 - q[..., i] ** 2 for i in range(3)], axis=-1)


def constraint(x, v):
    """I = sum_{i<j} v_i v_j - V; zero on physical solutions."""
    S = v.sum(axis=-1)
    return (S ** 2 - np.sum(v ** 2, axis=-1)) / 2.0 - potential(x)


def beta_to_x(bp, bm):
    """(beta_+, beta_-) at Omega = 0 to x = ln(a_i^2)."""
    return np.stack([-2.0 * bp, bp + _SQ3 * bm, bp - _SQ3 * bm], axis=-1)


def x_to_beta(x):
    return (x[..., 1] + x[..., 2] - 2.0 * x[..., 0]) / 6.0, (x[..., 1] - x[..., 2]) / (2.0 * _SQ3)


def turnaround_state(bp, bm):
    """Maximum-expansion initial state at (beta_+, beta_-): returns (x, v, V); v is meaningful only where V < 0."""
    bp, bm = np.asarray(bp, dtype=float), np.asarray(bm, dtype=float)
    x = beta_to_x(bp, bm)
    V = potential(x)
    r = np.hypot(bp, bm)
    r = np.where(r == 0, 1.0, r)
    sp = np.sqrt(np.maximum(-V, 0.0) / 3.0)
    dbp, dbm = sp * bp / r, sp * bm / r
    v = dbp[..., None] * np.array([-2.0, 1.0, 1.0]) + (_SQ3 * dbm)[..., None] * np.array([0.0, 1.0, -1.0])
    return x, v, V


def exit_basins(bp, bm, h0=0.01, nmax=30000, kappa=1e-3, p_exit=0.95, return_info=False, return_time=False, return_epochs=False):
    """Exit axis (1, 2, 3 = a, b, c) of the maximum-expansion state at each (beta_+, beta_-).

    Labels: 0 = no exit within ``nmax`` steps, -1 = not a valid initial state (V >= 0).  With ``return_info``
    also returns the largest relative constraint violation seen (sampled every 1000 steps); with ``return_time`` also the BKL
    time tau elapsed (the exit time where label > 0, the time reached at ``nmax`` steps where label = 0, nan for invalid states); with ``return_epochs`` also the number of Kasner epochs
    entered (|V| < kappa * sum x'^2, left again above 10 kappa), counted up to the exit or up to ``nmax`` steps (-1 for invalid
    states).  Returns label, then drift, time, epochs, as requested."""
    shape = np.shape(bp)
    x0, v0, V0 = turnaround_state(np.ravel(bp), np.ravel(bm))
    label = np.where(V0 < 0, 0, -1)
    act = np.where(V0 < 0)[0]
    x, v = x0[act].copy(), v0[act].copy()
    out = np.zeros(act.size, dtype=int)
    alive = np.arange(act.size)
    drift = 0.0
    T = np.zeros(act.size)
    nep = np.zeros(act.size, dtype=int)
    inK = np.zeros(act.size, dtype=bool)
    for it in range(nmax):
        if alive.size == 0:
            break
        xa, va = x[alive], v[alive]
        h = (h0 / (1.0 + np.max(np.abs(va), axis=1)))[:, None]
        k1x, k1v = va, force(xa)
        k2x, k2v = va + 0.5 * h * k1v, force(xa + 0.5 * h * k1x)
        k3x, k3v = va + 0.5 * h * k2v, force(xa + 0.5 * h * k2x)
        k4x, k4v = va + h * k3v, force(xa + h * k3x)
        xa = xa + h / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
        va = va + h / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        x[alive], v[alive] = xa, va
        T[alive] += h[:, 0]
        S = va.sum(axis=1)
        sq = np.sum(va ** 2, axis=1)
        Vc = potential(xa)
        p = va / np.where(S == 0, 1.0, S)[:, None]
        kap = np.abs(Vc) / np.maximum(sq, 1e-300)
        ex = (kap < kappa) & (p.max(axis=1) > p_exit) & (S < 0)
        if return_epochs:
            now = np.where(inK[alive], kap < 10 * kappa, kap < kappa)
            nep[alive] += (now & ~inK[alive]) & (S < 0)
            inK[alive] = now
        out[alive[ex]] = np.argmax(p[ex], axis=1) + 1
        if return_info and it % 1000 == 0:
            drift = max(drift, float(np.max(np.abs((S ** 2 - sq) / 2.0 - Vc) / np.maximum(np.abs(Vc) + sq, 1e-12))))
        alive = alive[~ex]
    label[act] = out
    label = label.reshape(shape)
    res = [label]
    if return_info:
        res.append(drift)
    if return_time:
        tau = np.full(label.size, np.nan)
        tau[act] = T
        res.append(tau.reshape(shape))
    if return_epochs:
        ne = np.full(label.size, -1)
        ne[act] = nep
        res.append(ne.reshape(shape))
    return res[0] if len(res) == 1 else tuple(res)


def mixmaster_picture(n=200, half=1.2, **kw) -> Picture:
    """Exit-basin picture on the window [-half, half]^2 of the (beta_+, beta_-) plane (n x n, indexed [i, j] = (beta_+,
    beta_-)).  Label 1..3 = exit axis, 0 = no exit, -1 = invalid initial state (V >= 0)."""
    ax = (np.arange(n) + 0.5) / n * 2 * half - half
    BP, BM = np.meshgrid(ax, ax, indexing="ij")
    lab = exit_basins(BP, BM, **kw)
    return Picture(lab.astype(float), lab, np.zeros(lab.shape, dtype=int),
                   dict(kind="mixmaster", n=n, half=half, **kw))
