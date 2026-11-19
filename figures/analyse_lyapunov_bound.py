"""Finite-time Lyapunov exponents of the section-x=0 orbits (E = 1) for the continuum potential x^2 y^2 / 2 and the Wilson-plaquette
potential of spacing eps, against the Hessian bound.  For V_eps = (1/eps^2) H(theta, phi) the eigenvalues of H lie in [-4, 4] (grid 1500^2,
the minimum -4 is attained exactly at theta = phi = pi), so the tangent flow  d'' = -Hess V d  grows at most like exp(2 t / eps).
Usage: python3 figures/analyse_lyapunov_bound.py [n_states] [T]"""
import sys, csv
sys.path.insert(0, ".")
import numpy as np

def hess(x, y, eps):
    if eps is None:
        return y * y, 2 * x * y, x * x
    return (4 / eps ** 2 * np.cos(eps * x) * np.sin(eps * y / 2) ** 2,
            2 / eps ** 2 * np.sin(eps * x) * np.sin(eps * y),
            4 / eps ** 2 * np.cos(eps * y) * np.sin(eps * x / 2) ** 2)

def force(x, y, eps):
    if eps is None:
        return -x * y * y, -x * x * y
    return (-4 / eps ** 3 * np.sin(eps * x) * np.sin(eps * y / 2) ** 2, -4 / eps ** 3 * np.sin(eps * y) * np.sin(eps * x / 2) ** 2)

def lyap(eps, n, T, dt=0.002, E=1.0, seed=0, ymax=2.0, renorm=100):
    rng = np.random.default_rng(seed)
    y = rng.uniform(-ymax, ymax, n); py = rng.uniform(-np.sqrt(2 * E), np.sqrt(2 * E), n)
    s = np.stack([np.zeros(n), y, np.sqrt(2 * E - py ** 2), py])        # x, y, px, py
    d = rng.normal(size=(4, n)); d /= np.linalg.norm(d, axis=0)           # tangent (dx, dy, dpx, dpy)
    acc = np.zeros(n); inst = np.zeros(n)

    def rhs(s, d):
        fx, fy = force(s[0], s[1], eps); a, b, c = hess(s[0], s[1], eps)
        return (np.stack([s[2], s[3], fx, fy]), np.stack([d[2], d[3], -(a * d[0] + b * d[1]), -(b * d[0] + c * d[1])]))
    for k in range(int(T / dt)):
        k1s, k1d = rhs(s, d); k2s, k2d = rhs(s + .5 * dt * k1s, d + .5 * dt * k1d)
        k3s, k3d = rhs(s + .5 * dt * k2s, d + .5 * dt * k2d); k4s, k4d = rhs(s + dt * k3s, d + dt * k3d)
        s = s + dt / 6 * (k1s + 2 * k2s + 2 * k3s + k4s); d = d + dt / 6 * (k1d + 2 * k2d + 2 * k3d + k4d)
        if k % 50 == 0:
            a, b, c = hess(s[0], s[1], eps); lm = (a + c) / 2 - np.sqrt(((a - c) / 2) ** 2 + b * b)
            inst = np.maximum(inst, np.sqrt(np.maximum(-lm, 0)))
        if k % renorm == renorm - 1:
            nr = np.linalg.norm(d, axis=0); acc += np.log(nr); d /= nr
    return acc / T, inst


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    T = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
    rows = []
    for eps in (None, 0.4, 0.2, 0.1):
        lam, inst = lyap(eps, n, T)
        bound = "inf" if eps is None else 2 / eps
        print(f"eps={eps}: Lyapunov exponent mean {lam.mean():.3f}, median {np.median(lam):.3f}, 99th pct {np.percentile(lam,99):.3f}, max {lam.max():.3f};"
              f" orbit-maximum of sqrt(-lambda_min(Hess)): mean {inst.mean():.3f}, max {inst.max():.3f}; global bound {bound}", flush=True)
        rows.append((eps, lam.mean(), np.median(lam), np.percentile(lam, 99), lam.max(), inst.mean(), inst.max(), bound))
    with open("figures/lyapunov_bound.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["eps", "lyap_mean", "lyap_median", "lyap_p99", "lyap_max", "orbit_sqrt_neg_lmin_mean", "orbit_sqrt_neg_lmin_max", "global_bound_2_over_eps"])
        w.writerows(rows)
