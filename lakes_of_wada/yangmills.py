"""Homogeneous SU(2) Yang-Mills mechanics, H = (px^2 + py^2)/2 + x^2 y^2 / 2, as an open scattering system.

The potential has four valleys (the half axes +x, -x, +y, -y) of constant transverse width ~ 1/|x|.  A particle in a valley is not
free: with the adiabatic invariant J = E_y / |x| the motion along the valley feels the potential J |x|, so there is no true escape
to infinity, and the exit has to be defined.  Here it is *first passage to the circle of radius R*: the label is the half axis
(1: +x, 2: +y, 3: -x, 4: -y) on which the orbit crosses it.  For R above (8E)^(1/4) = 1.68 (E = 1) the circle lies outside the
region x^2 y^2 < 2E except in the valleys, so the label is unambiguous.  The label depends on R (see YANGMILLS.md).
The energy can be set to 1 by the scaling x -> s x, t -> t / s, E -> s^4 E.

Initial conditions: the section x = 0 (V = 0 there for every y), coordinates (y, p_y) with p_y^2 < 2E and p_x = +sqrt(2E - p_y^2).
"""
import numpy as np


def potential(x, y):
    return 0.5 * x * x * y * y


def energy(x, y, px, py):
    return 0.5 * (px * px + py * py) + potential(x, y)


def exit_basins(y, py, E=1.0, R=4.0, dt=0.004, t_max=100.0, return_time=False, return_info=False):
    """Exit half axis (1 +x, 2 +y, 3 -x, 4 -y) of the orbit started at (x, y, px, py) = (0, y, +sqrt(2E - py^2), py).

    Labels: 0 = no passage through radius R before ``t_max``, -1 = invalid (py^2 >= 2E).  Velocity Verlet.  With
    ``return_time`` also the passage time (nan where undecided or invalid); with ``return_info`` also the largest relative energy drift."""
    shape = np.shape(y)
    y0 = np.ravel(np.asarray(y, float))
    p0 = np.ravel(np.asarray(py, float))
    valid = p0 * p0 < 2 * E
    label = np.where(valid, 0, -1)
    idx = np.where(valid)[0]
    x = np.zeros(idx.size)
    yy = y0[idx].copy()
    py_ = p0[idx].copy()
    px = np.sqrt(2 * E - py_ ** 2)
    out = np.zeros(idx.size, dtype=int)
    tex = np.full(idx.size, np.nan)
    alive = np.arange(idx.size)
    drift = 0.0
    nsteps = int(t_max / dt)
    fx = -x * yy * yy
    fy = -x * x * yy
    for it in range(nsteps):
        if alive.size == 0:
            break
        a = alive
        px[a] += 0.5 * dt * fx[a]; py_[a] += 0.5 * dt * fy[a]
        x[a] += dt * px[a]; yy[a] += dt * py_[a]
        fx[a] = -x[a] * yy[a] ** 2; fy[a] = -x[a] ** 2 * yy[a]
        px[a] += 0.5 * dt * fx[a]; py_[a] += 0.5 * dt * fy[a]
        r2 = x[a] ** 2 + yy[a] ** 2
        ex = r2 > R * R
        if ex.any():
            e = a[ex]
            ax_x = np.abs(x[e]) > np.abs(yy[e])
            lab = np.where(ax_x, np.where(x[e] > 0, 1, 3), np.where(yy[e] > 0, 2, 4))
            out[e] = lab
            tex[e] = (it + 1) * dt
        if return_info and it % 500 == 0:
            drift = max(drift, float(np.max(np.abs(energy(x[a], yy[a], px[a], py_[a]) - E) / E)))
        alive = a[~ex]
    label[idx] = out
    res = [label.reshape(shape)]
    if return_time:
        t = np.full(label.size, np.nan)
        t[idx] = tex
        res.append(t.reshape(shape))
    if return_info:
        res.append(drift)
    return res[0] if len(res) == 1 else tuple(res)


def section_picture(n=200, half=2.0, E=1.0, **kw):
    """Labels on an n x n grid of (y, p_y) in [-half, half] x (-sqrt(2E), sqrt(2E)) (index [i, j] = (y_i, p_j))."""
    ys = (np.arange(n) + 0.5) / n * 2 * half - half
    pm = np.sqrt(2 * E)
    ps = (np.arange(n) + 0.5) / n * 2 * pm - pm
    Y, P = np.meshgrid(ys, ps, indexing="ij")
    return exit_basins(Y, P, E=E, **kw), ys, ps


# ------------------------------------------------------------------------------------------ gauge-Higgs toy with true final states

def higgs_potential(x, y, lam=1.0, v=1.0):
    """V = x^2 y^2 / 2 + (lam/4) (x^2 + y^2 - v^2)^2: the Yang-Mills valleys capped by a Higgs-like ring; four degenerate minima
    (+-v, 0), (0, +-v) with V = 0.  A toy chosen for this repository, not a derived truncation."""
    return 0.5 * x * x * y * y + 0.25 * lam * (x * x + y * y - v * v) ** 2


def higgs_force(x, y, lam=1.0, v=1.0):
    s = lam * (x * x + y * y - v * v)
    return -x * y * y - s * x, -x * x * y - s * y


def higgs_basins(x0, y0, lam=1.0, gamma=0.1, v=1.0, dt=0.01, t_max=400.0, tol=0.03, return_time=False):
    """Final well (1: +x, 2: +y, 3: -x, 4: -y) of the damped motion  q'' = -grad V - gamma q'  released from rest at (x0, y0).

    A state is settled when it is within ``tol`` of a minimum and its speed is below ``tol``; the label is that minimum.  Labels: 0 =
    not settled by ``t_max``.  RK4 with the friction term included."""
    shape = np.shape(x0)
    x = np.ravel(np.asarray(x0, float)).copy()
    y = np.ravel(np.asarray(y0, float)).copy()
    px = np.zeros_like(x)
    py = np.zeros_like(x)
    out = np.zeros(x.size, dtype=int)
    tex = np.full(x.size, np.nan)
    alive = np.arange(x.size)
    wells = np.array([[v, 0.0], [0.0, v], [-v, 0.0], [0.0, -v]])

    def acc(xa, ya, pxa, pya):
        fx, fy = higgs_force(xa, ya, lam, v)
        return fx - gamma * pxa, fy - gamma * pya

    for it in range(int(t_max / dt)):
        if alive.size == 0:
            break
        xa, ya, pxa, pya = x[alive], y[alive], px[alive], py[alive]
        k1x, k1y = pxa, pya; k1px, k1py = acc(xa, ya, pxa, pya)
        k2x, k2y = pxa + .5 * dt * k1px, pya + .5 * dt * k1py; k2px, k2py = acc(xa + .5 * dt * k1x, ya + .5 * dt * k1y, k2x, k2y)
        k3x, k3y = pxa + .5 * dt * k2px, pya + .5 * dt * k2py; k3px, k3py = acc(xa + .5 * dt * k2x, ya + .5 * dt * k2y, k3x, k3y)
        k4x, k4y = pxa + dt * k3px, pya + dt * k3py; k4px, k4py = acc(xa + dt * k3x, ya + dt * k3y, k4x, k4y)
        xa = xa + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x); ya = ya + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
        pxa = pxa + dt / 6 * (k1px + 2 * k2px + 2 * k3px + k4px); pya = pya + dt / 6 * (k1py + 2 * k2py + 2 * k3py + k4py)
        x[alive], y[alive], px[alive], py[alive] = xa, ya, pxa, pya
        if it % 10 == 0:
            d = np.sqrt((xa[:, None] - wells[:, 0]) ** 2 + (ya[:, None] - wells[:, 1]) ** 2)
            near = d.min(axis=1) < tol
            slow = pxa ** 2 + pya ** 2 < tol ** 2
            done = near & slow
            if done.any():
                out[alive[done]] = np.argmin(d[done], axis=1) + 1
                tex[alive[done]] = (it + 1) * dt
                alive = alive[~done]
    res = [out.reshape(shape)]
    if return_time:
        res.append(tex.reshape(shape))
    return res[0] if len(res) == 1 else tuple(res)


def higgs_picture(n=200, half=1.8, **kw):
    """Basins on an n x n grid of released-from-rest positions in [-half, half]^2 (index [i, j] = (x_i, y_j))."""
    a = (np.arange(n) + 0.5) / n * 2 * half - half
    X, Yg = np.meshgrid(a, a, indexing="ij")
    return higgs_basins(X, Yg, **kw), a
