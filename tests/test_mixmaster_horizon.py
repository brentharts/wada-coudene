import numpy as np
from lakes_of_wada import mixmaster_horizon as MH
from lakes_of_wada.mixmaster import turnaround_state


def kasner_start(u, scale=-12.0):
    """Pure-Kasner start: tiny scale factors (curvature negligible for a long BKL time), S = -1, exponents from u."""
    p = np.array([-u / (1 + u + u * u), (1 + u) / (1 + u + u * u), u * (1 + u) / (1 + u + u * u)])
    assert abs(p.sum() - 1) < 1e-12 and abs((p ** 2).sum() - 1) < 1e-12
    x0 = np.full((1, 3), scale)
    v0 = (-1.0 * p)[None, :]
    return x0, v0, p


def test_kasner_horizon_is_one_over_one_minus_p():
    x0, v0, p = kasner_start(2.0)
    orb = MH.horizon_orbits(None, None, x0=x0, v0=v0, tau_budget=10.0, tau_cap=10.0, h0=0.02, thin=5)[0]
    # early in the run the curvature is still negligible and Kasner continues to the singularity within the tail estimate
    R = orb["R"][3]
    assert np.allclose(R, 1.0 / (1.0 - p), rtol=2e-3)


def test_reverse_sum_matches_direct_quadrature():
    # a_i = e^{x_i/2}; for a Kasner orbit H_i(tau) = integral_tau^inf sqrt(q_j q_k) dtau has a closed form
    x0, v0, p = kasner_start(1.3)
    orb = MH.horizon_orbits(None, None, x0=x0, v0=v0, tau_budget=8.0, tau_cap=8.0, h0=0.02, thin=5)[0]
    S = -1.0
    e = np.exp(0.5 * orb["x"])
    g = np.stack([e[:, 1] * e[:, 2], e[:, 2] * e[:, 0], e[:, 0] * e[:, 1]], axis=1)
    exact = g / (0.5 * abs(S) * (1.0 - p))
    assert np.allclose(orb["H"][2], exact[2], rtol=2e-3)


def test_real_orbit_has_epochs_and_finite_horizons():
    rng = np.random.default_rng(3)
    for bp, bm in rng.uniform(-0.8, 0.8, (30, 2)):
        x, v, V = turnaround_state(bp, bm)
        if V < 0:
            break
    orbs = MH.horizon_orbits(np.array([bp]), np.array([bm]), tau_budget=40.0, tau_cap=80.0)
    assert len(orbs) == 1
    o = orbs[0]
    assert np.all(np.isfinite(o["R"])) and np.all(o["R"] > 0)
    assert len(MH.epochs(o)) >= 2
