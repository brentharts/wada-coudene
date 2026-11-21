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


def wilson_potential(x, y, eps):
    """Wilson-plaquette version of x^2 y^2 / 2 for constant SU(2) links U_1 = exp(-i eps x s_1 / 2), U_2 = exp(-i eps y s_2 / 2):
    V = (8 / eps^4) sin^2(eps x / 2) sin^2(eps y / 2) = (4 / eps^4) (1 - (1/2) Re Tr U_p), which is ||F||^2_HS for F = [A_1, A_2] in the limit.  Periodic, bounded by 8 / eps^4, and
    -> x^2 y^2 / 2 as eps -> 0 (the small-plaquette expansion, Lockwood 2026, section 15.4)."""
    return 8.0 / eps ** 4 * np.sin(0.5 * eps * x) ** 2 * np.sin(0.5 * eps * y) ** 2


def energy(x, y, px, py, eps=None):
    v = potential(x, y) if eps is None else wilson_potential(x, y, eps)
    return 0.5 * (px * px + py * py) + v


def _force(x, y, eps):
    if eps is None:
        return -x * y * y, -x * x * y
    return (-4.0 / eps ** 3 * np.sin(eps * x) * np.sin(0.5 * eps * y) ** 2,
            -4.0 / eps ** 3 * np.sin(eps * y) * np.sin(0.5 * eps * x) ** 2)


def exit_basins(y, py, E=1.0, R=4.0, dt=0.004, t_max=100.0, return_time=False, return_info=False, eps=None):
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
    fx, fy = _force(x, yy, eps)
    for it in range(nsteps):
        if alive.size == 0:
            break
        a = alive
        px[a] += 0.5 * dt * fx[a]; py_[a] += 0.5 * dt * fy[a]
        x[a] += dt * px[a]; yy[a] += dt * py_[a]
        fx[a], fy[a] = _force(x[a], yy[a], eps)
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
            drift = max(drift, float(np.max(np.abs(energy(x[a], yy[a], px[a], py_[a], eps) - E) / E)))
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


# ------------------------------------------------------------------------------------------ Langevin version of the Higgs toy

def higgs_langevin(x0, y0, T, gamma=0.1, lam=1.0, v=1.0, M=20, dt=0.01, t_max=None, r_well=0.3, e_settle=0.04, seed=0):
    """Probability of each final well for  q'' = -grad V - gamma q' + sqrt(2 gamma T) xi,  released from rest at (x0, y0).

    ``M`` noise realisations per start (BAOAB splitting, so the invariant measure is exp(-(V + p^2/2)/T) up to O(dt^2)).  A realisation
    is settled when it is within ``r_well`` of a minimum with mechanical energy below ``e_settle`` (the saddle between wells is near
    V = 0.084 for lam = 1, v = 1).  Returns an array of shape x0.shape + (5,): index 0 = not settled by ``t_max`` (default
    max(25 / gamma, 40)), indices 1..4 = wells +x, +y, -x, -y.  T = 0 is the deterministic limit."""
    rng = np.random.default_rng(seed)
    shape = np.shape(x0)
    x = np.repeat(np.ravel(np.asarray(x0, float)), M)
    y = np.repeat(np.ravel(np.asarray(y0, float)), M)
    px = np.zeros_like(x); py = np.zeros_like(x)
    nst = x.size // M
    out = np.zeros(x.size, dtype=int)
    wells = np.array([[v, 0.0], [0.0, v], [-v, 0.0], [0.0, -v]])
    if t_max is None:
        t_max = max(25.0 / gamma, 40.0)
    c = np.exp(-gamma * dt); s = np.sqrt(T * (1 - c * c))
    alive = np.arange(x.size)
    fx, fy = higgs_force(x, y, lam, v)
    for it in range(int(t_max / dt)):
        if alive.size == 0:
            break
        a = alive
        px[a] += 0.5 * dt * fx[a]; py[a] += 0.5 * dt * fy[a]
        x[a] += 0.5 * dt * px[a]; y[a] += 0.5 * dt * py[a]
        if T > 0:
            px[a] = c * px[a] + s * rng.standard_normal(a.size); py[a] = c * py[a] + s * rng.standard_normal(a.size)
        else:
            px[a] *= c; py[a] *= c
        x[a] += 0.5 * dt * px[a]; y[a] += 0.5 * dt * py[a]
        fx[a], fy[a] = higgs_force(x[a], y[a], lam, v)
        px[a] += 0.5 * dt * fx[a]; py[a] += 0.5 * dt * fy[a]
        if it % 10 == 0:
            xa, ya = x[a], y[a]
            d = np.sqrt((xa[:, None] - wells[:, 0]) ** 2 + (ya[:, None] - wells[:, 1]) ** 2)
            em = 0.5 * (px[a] ** 2 + py[a] ** 2) + higgs_potential(xa, ya, lam, v)
            done = (d.min(axis=1) < r_well) & (em < e_settle)
            if done.any():
                out[a[done]] = np.argmin(d[done], axis=1) + 1
                alive = a[~done]
    cnt = np.zeros((nst, 5))
    lab = out.reshape(nst, M)
    for k in range(5):
        cnt[:, k] = (lab == k).mean(axis=1)
    return cnt.reshape(shape + (5,))


# ------------------------------------------------------------------------------------------ three-variable model, z free

def potential3(x, y, z):
    return 0.5 * (x * x * y * y + y * y * z * z + z * z * x * x)


def exit_basins3(y, py, z=0.0, pz=0.0, E=1.0, R=2.0, dt=0.004, t_max=60.0):
    """First-passage exit of H = p^2/2 + (x^2 y^2 + y^2 z^2 + z^2 x^2)/2 from x = 0, started at (0, y, z) with momentum (px, py, pz),
    px = +sqrt(2E - py^2 - pz^2 - V) (V = y^2 z^2 / 2).  Labels: 1 +x, 2 +y, 3 -x, 4 -y, 5 +z, 6 -z (the axis with the largest
    |coordinate| when the orbit crosses radius R); 0 = undecided, -1 = invalid (px^2 < 0).  Velocity Verlet.  z = pz = 0 is an invariant
    subspace, in which this is ``exit_basins``."""
    shape = np.broadcast(np.asarray(y), np.asarray(py), np.asarray(z), np.asarray(pz)).shape
    Y, PY, Z, PZ = [np.ravel(np.broadcast_to(np.asarray(a, float), shape)).copy() for a in (y, py, z, pz)]
    px2 = 2 * E - PY ** 2 - PZ ** 2 - 0.5 * Y * Y * Z * Z
    valid = px2 > 0
    label = np.where(valid, 0, -1)
    idx = np.where(valid)[0]
    X = np.zeros(idx.size); yy = Y[idx]; zz = Z[idx]
    px = np.sqrt(px2[idx]); py_ = PY[idx].copy(); pz_ = PZ[idx].copy()
    out = np.zeros(idx.size, dtype=int)
    alive = np.arange(idx.size)

    def force(a):
        x_, y_, z_ = X[a], yy[a], zz[a]
        return -x_ * (y_ * y_ + z_ * z_), -y_ * (x_ * x_ + z_ * z_), -z_ * (x_ * x_ + y_ * y_)
    f = [np.zeros(idx.size) for _ in range(3)]
    f[0][:], f[1][:], f[2][:] = force(alive)
    for it in range(int(t_max / dt)):
        if alive.size == 0:
            break
        a = alive
        px[a] += 0.5 * dt * f[0][a]; py_[a] += 0.5 * dt * f[1][a]; pz_[a] += 0.5 * dt * f[2][a]
        X[a] += dt * px[a]; yy[a] += dt * py_[a]; zz[a] += dt * pz_[a]
        f[0][a], f[1][a], f[2][a] = force(a)
        px[a] += 0.5 * dt * f[0][a]; py_[a] += 0.5 * dt * f[1][a]; pz_[a] += 0.5 * dt * f[2][a]
        ex = X[a] ** 2 + yy[a] ** 2 + zz[a] ** 2 > R * R
        if ex.any():
            e = a[ex]
            c = np.stack([X[e], yy[e], zz[e]])
            k = np.argmax(np.abs(c), axis=0)
            sgn = np.take_along_axis(c, k[None], 0)[0] > 0
            out[e] = np.array([[3, 1], [4, 2], [6, 5]])[k, sgn.astype(int)]
        alive = a[~ex]
    label[idx] = out
    return label.reshape(shape)


def exit_labels3(pos, mom, R=3.0, dt=0.005, t_max=150.0):
    """First-passage exit of H = p^2/2 + (x^2 y^2 + y^2 z^2 + z^2 x^2)/2 from arbitrary starts.  ``pos`` and ``mom`` are (n, 3) arrays.
    Labels as in ``exit_basins3`` (1 +x, 2 +y, 3 -x, 4 -y, 5 +z, 6 -z, the axis with the largest |coordinate| at radius R; 0 = undecided).
    Velocity Verlet.  The flow is invariant under the octahedral group (permutations and sign changes of x, y, z), and ``exit_labels3``
    is equivariant under it up to rounding, which chaos can amplify."""
    P = np.array(pos, float).reshape(-1, 3).copy()
    Q = np.array(mom, float).reshape(-1, 3).copy()
    n = P.shape[0]
    out = np.zeros(n, dtype=int)
    alive = np.arange(n)

    def force(p):
        x, y, z = p[:, 0], p[:, 1], p[:, 2]
        return np.stack([-x * (y * y + z * z), -y * (x * x + z * z), -z * (x * x + y * y)], axis=1)
    F = force(P)
    for _ in range(int(t_max / dt)):
        if alive.size == 0:
            break
        a = alive
        Q[a] += 0.5 * dt * F[a]
        P[a] += dt * Q[a]
        F[a] = force(P[a])
        Q[a] += 0.5 * dt * F[a]
        ex = np.einsum("ij,ij->i", P[a], P[a]) > R * R
        if ex.any():
            e = a[ex]
            k = np.argmax(np.abs(P[e]), axis=1)
            sgn = P[e][np.arange(e.size), k] > 0
            out[e] = np.array([[3, 1], [4, 2], [6, 5]])[k, sgn.astype(int)]
        alive = a[~ex]
    return out


# the plane x + y + z = 0 with the orthonormal frame (E1, E2); the six half axes project to the vertices of a regular hexagon
E1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
E2 = np.array([1.0, 1.0, -2.0]) / np.sqrt(6.0)
# angle (degrees) of the projection of each label's half axis in the frame (E1, E2): +x 30, ... computed, see HEX_ORDER
HEX_ORDER = None   # filled below


def _hex_order():
    ang = {}
    for lab, v in ((1, (1, 0, 0)), (2, (0, 1, 0)), (3, (-1, 0, 0)), (4, (0, -1, 0)), (5, (0, 0, 1)), (6, (0, 0, -1))):
        v = np.array(v, float)
        ang[lab] = float(np.degrees(np.arctan2(v @ E2, v @ E1)) % 360.0)
    return sorted(ang, key=ang.get), ang


HEX_ORDER, HEX_ANGLES = _hex_order()


def plane_start(u, v):
    """Points (u, v) of the plane x + y + z = 0 as (n, 3) positions."""
    u = np.asarray(u, float).ravel(); v = np.asarray(v, float).ravel()
    return u[:, None] * E1[None, :] + v[:, None] * E2[None, :]


def exit_plane_rest(u, v, R=3.0, dt=0.005, t_max=150.0):
    """Labels of the orbits released from rest at the points (u, v) of the plane x + y + z = 0."""
    P = plane_start(u, v)
    return exit_labels3(P, np.zeros_like(P), R=R, dt=dt, t_max=t_max).reshape(np.shape(u))
