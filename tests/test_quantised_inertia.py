import numpy as np
from lakes_of_wada import quantised_inertia as Q
from lakes_of_wada.scattering import integrate_exits


def test_laws_have_the_stated_limits():
    g = np.array([1e-9, 1e-3, 1.0, 1e3])
    a0 = 1e-3
    mond = Q.inertia_factor(g, a0, "mond") * g
    # Eq. 14: m_i = m/(1 + a0/a)  <=>  g = a / (1 + a0/a)  (F = m_i a)
    assert np.allclose(mond * mond / (mond + a0), g, rtol=1e-9)
    assert np.allclose(mond[:1], np.sqrt(a0 * g[:1]), rtol=1e-3)          # deep MoND limit
    assert np.allclose(Q.inertia_factor(g, a0, "newton"), 1.0)
    qi1 = Q.inertia_factor(g[2:], a0, "qi1") * g[2:]
    assert np.allclose(qi1, g[2:] + a0, rtol=1e-6)                         # |a| = g + a0 (Eq. 1, positive branch)
    big = Q.inertia_factor(g[-1:], a0, "mond") * g[-1:]
    assert np.allclose(big, g[-1:], rtol=1e-5)                              # Newton recovered at large g


def test_eq1_branch_satisfies_the_force_law():
    g = np.array([0.3, 2.0])
    a0 = 0.05
    a = Q.inertia_factor(g, a0, "qi1") * g
    assert np.allclose(g, (1 - a0 / a) * a, rtol=1e-5)                      # F = m_i a with m_i = m(1 - a0/|a|)


def test_newton_conserves_energy():
    s = Q.system(R_H=10.0, law="newton")
    st = Q.beam_states(20, speed=3.0)
    _, _, fin = integrate_exits(s, st, dt=0.01, t_max=3.0)
    e0, e1 = s.energy(st), s.energy(fin)
    assert np.max(np.abs(e1 - e0)) < 5e-3


def test_cyclic_symmetry_of_the_picture_for_tangential_launch():
    rot = np.array([[np.cos(2 * np.pi / 3), -np.sin(2 * np.pi / 3)], [np.sin(2 * np.pi / 3), np.cos(2 * np.pi / 3)]])
    rng = np.random.default_rng(0)
    p = rng.uniform(-2, 2, (300, 2))
    for law in ("newton", "mond", "qi1"):
        s = Q.system(R_H=8.0, law=law)
        a, _, _ = integrate_exits(s, Q.launch(p, 3.0), dt=0.02, t_max=100)
        b, _, _ = integrate_exits(s, Q.launch(p @ rot.T, 3.0), dt=0.02, t_max=100)
        ok = (a > 0) & (b > 0)
        assert ok.sum() > 150
        # rotating the start by +120 degrees moves the exit sector label 1 -> 2 -> 3 -> 1
        assert np.mean(((a[ok] % 3) + 1) == b[ok]) > 0.95


def test_horizon_is_where_the_rindler_horizon_equals_the_cosmic_one():
    R_H = 7.0
    a0 = 1.0 / R_H                       # c = 1
    assert abs(1.0 / a0 - R_H) < 1e-12   # Rindler horizon distance c^2/a at a = a0 is R_H = Theta/2


def test_hh_with_zero_a0_is_the_repository_henon_heiles():
    from lakes_of_wada.scattering import henon_heiles_section, henon_heiles
    y0, valid = henon_heiles_section(40, 0.25)
    ref, _, _ = integrate_exits(henon_heiles(5.0), y0[valid], 0.05, 100.0)
    new, _, _ = integrate_exits(Q.hh_system(3, 5.0, "mond", a0=0.0), y0[valid], 0.05, 100.0)
    assert np.array_equal(ref, new)


def test_hh_inertia_law_changes_labels_and_breaks_energy_conservation():
    y0, valid = __import__("lakes_of_wada.scattering", fromlist=["x"]).henon_heiles_section(30, 0.25)
    s = y0[valid]
    for law, conserves in (("newton", True), ("mond", False)):
        sysm = Q.hh_system(3, 5.0, law, a0=0.05)
        _, _, fin = integrate_exits(sysm, s, 0.02, 5.0)
        drift = np.max(np.abs(sysm.energy(fin) - sysm.energy(s)))
        assert (drift < 1e-4) == conserves


def test_noise_free_limit_is_the_deterministic_map_and_noise_changes_labels():
    from lakes_of_wada.scattering import henon_heiles_section
    y0, valid = henon_heiles_section(12, 0.25)
    s = y0[valid]
    ref, _, _ = integrate_exits(Q.hh_system(3, 5.0, "mond", 0.05), s, 0.05, 100.0)
    p0 = Q.hh_noisy_labels(s, "mond", 0.05, D0=0.0, n_real=2, t_max=100.0)
    assert np.array_equal(np.argmax(p0, axis=1), ref) and np.all(p0.max(axis=1) == 1.0)
    p1 = Q.hh_noisy_labels(s, "mond", 0.05, D0=1e-2, n_real=6, t_max=100.0)
    assert np.mean(p1.max(axis=1) < 1.0) > 0.05
