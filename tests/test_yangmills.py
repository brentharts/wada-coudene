import numpy as np
from lakes_of_wada import yangmills as Y


def _sample(n=400, seed=0):
    rng = np.random.default_rng(seed)
    return rng.uniform(-2, 2, n), rng.uniform(-1.4, 1.4, n)


def test_invalid_states_are_marked():
    lab = Y.exit_basins(np.array([0.0, 0.5]), np.array([1.5, 0.0]), R=2.0)
    assert lab[0] == -1 and lab[1] > 0


def test_symmetry_y_p_to_minus_y_minus_p_swaps_the_y_channels():
    y, p = _sample()
    a = Y.exit_basins(y, p, R=3.0, t_max=60)
    b = Y.exit_basins(-y, -p, R=3.0, t_max=60)
    swap = {0: 0, 1: 1, 2: 4, 3: 3, 4: 2}
    assert np.array_equal(np.vectorize(swap.get)(b), a)


def test_energy_is_conserved():
    y, p = _sample(200, 1)
    _, drift = Y.exit_basins(y, p, R=3.0, t_max=60, return_info=True)
    assert drift < 1e-3


def test_all_four_channels_are_used_and_labels_converge_at_small_radius():
    y, p = _sample(1500, 2)
    a = Y.exit_basins(y, p, R=2.0, dt=0.002, t_max=60)
    b = Y.exit_basins(y, p, R=2.0, dt=0.001, t_max=60)
    m = (a > 0) & (b > 0)
    assert m.mean() > 0.99
    assert np.mean(a[m] != b[m]) < 0.01
    for k in (1, 2, 3, 4):
        assert np.mean(a == k) > 0.1


def test_higgs_wells_are_degenerate_minima():
    for x, y in ((1, 0), (0, 1), (-1, 0), (0, -1)):
        assert abs(Y.higgs_potential(x, y)) < 1e-14
        fx, fy = Y.higgs_force(x, y)
        assert abs(fx) < 1e-14 and abs(fy) < 1e-14
    assert Y.higgs_potential(0.0, 0.0) > 0 and Y.higgs_potential(0.7, 0.7) > 0


def test_higgs_basins_have_the_symmetry_of_the_square_and_settle():
    rng = np.random.default_rng(3)
    x, y = rng.uniform(-1.8, 1.8, 300), rng.uniform(-1.8, 1.8, 300)
    kw = dict(gamma=0.3, t_max=120)
    a = Y.higgs_basins(x, y, **kw)
    assert np.mean(a == 0) < 0.03
    sw = {0: 0, 1: 2, 2: 1, 3: 4, 4: 3}            # (x, y) -> (y, x)
    b = Y.higgs_basins(y, x, **kw)
    assert np.array_equal(np.vectorize(sw.get)(b), a)
    fl = {0: 0, 1: 3, 2: 2, 3: 1, 4: 4}            # x -> -x
    c = Y.higgs_basins(-x, y, **kw)
    assert np.array_equal(np.vectorize(fl.get)(c), a)


def test_strong_friction_gives_smooth_basins():
    n = 60
    lab, _ = Y.higgs_picture(n=n, gamma=1.0, t_max=60)
    ok = (lab[:-1] > 0) & (lab[1:] > 0)
    assert np.mean(lab[:-1][ok] == lab[1:][ok]) > 0.95


# ------------------------------------------------------------------------------------------ Wilson-plaquette version (Lockwood 2026, section 15.4)

def _su2(axis, t):
    s = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)][axis]
    return np.cos(t / 2) * np.eye(2) - 1j * np.sin(t / 2) * s


def test_wilson_potential_is_the_plaquette():
    rng = np.random.default_rng(3)
    for _ in range(200):
        eps = rng.uniform(0.05, 1.5)
        x, y = rng.uniform(-4, 4, 2)
        U1, U2 = _su2(0, eps * x), _su2(1, eps * y)
        Up = U1 @ U2 @ np.linalg.inv(U1) @ np.linalg.inv(U2)
        v = 4.0 / eps ** 4 * (1 - 0.5 * np.trace(Up).real)
        assert abs(v - Y.wilson_potential(x, y, eps)) < 1e-9 * (1 + v)


def test_wilson_potential_limit_force_and_periodicity():
    x, y = np.array([0.7, -1.3, 2.0]), np.array([0.4, 0.9, -0.2])
    assert np.max(abs(Y.wilson_potential(x, y, 0.01) - Y.potential(x, y))) < 1e-4
    eps, h = 0.7, 1e-6
    fx, fy = Y._force(x, y, eps)
    nx = -(Y.wilson_potential(x + h, y, eps) - Y.wilson_potential(x - h, y, eps)) / (2 * h)
    ny = -(Y.wilson_potential(x, y + h, eps) - Y.wilson_potential(x, y - h, eps)) / (2 * h)
    assert np.max(abs(fx - nx)) < 1e-6 and np.max(abs(fy - ny)) < 1e-6
    assert np.allclose(Y.wilson_potential(x + 2 * np.pi / eps, y, eps), Y.wilson_potential(x, y, eps))
    assert np.max(Y.wilson_potential(np.linspace(-20, 20, 401), 0.5 * np.linspace(-20, 20, 401), eps)) <= 8 / eps ** 4 + 1e-12


def test_wilson_basins_approach_continuum_and_keep_symmetry():
    n = 40
    base, ys, ps = Y.section_picture(n=n, R=2.0, dt=0.004, t_max=60.0)
    lab, _, _ = Y.section_picture(n=n, R=2.0, dt=0.004, t_max=60.0, eps=0.02)
    ok = (base > 0) & (lab > 0)
    assert np.mean(base[ok] != lab[ok]) < 0.05
    sw = {2: 4, 4: 2, 1: 1, 3: 3}
    flipped = lab[::-1, ::-1]
    assert np.mean(np.vectorize(sw.get)(lab[lab > 0]) != flipped[lab > 0]) < 0.02


# ------------------------------------------------------------------------------------------ Hessian bounds (Romanov 2026, section 1.2 analogue)

def _hess_wilson(x, y, eps):
    return (4 / eps ** 2 * np.cos(eps * x) * np.sin(eps * y / 2) ** 2, 2 / eps ** 2 * np.sin(eps * x) * np.sin(eps * y),
            4 / eps ** 2 * np.cos(eps * y) * np.sin(eps * x / 2) ** 2)


def test_continuum_negative_curvature_is_bounded_by_xy():
    rng = np.random.default_rng(0)
    x, y = rng.normal(size=(2, 100000)) * 3
    a, b, c = y * y, 2 * x * y, x * x
    lmin = (a + c) / 2 - np.sqrt(((a - c) / 2) ** 2 + b * b)
    assert np.all(-lmin <= np.abs(x * y) + 1e-9)          # hence -lambda_min <= sqrt(2E) on x^2 y^2 / 2 <= E


def test_wilson_hessian_formulas_and_bound():
    eps, h = 0.6, 1e-4
    x, y = 0.9, -1.7
    V = lambda u, v: Y.wilson_potential(u, v, eps)
    vxx = (V(x + h, y) - 2 * V(x, y) + V(x - h, y)) / h ** 2
    vyy = (V(x, y + h) - 2 * V(x, y) + V(x, y - h)) / h ** 2
    vxy = (V(x + h, y + h) - V(x + h, y - h) - V(x - h, y + h) + V(x - h, y - h)) / (4 * h * h)
    a, b, c = _hess_wilson(x, y, eps)
    assert abs(a - vxx) < 1e-3 and abs(c - vyy) < 1e-3 and abs(b - vxy) < 1e-3
    t = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    T, P = np.meshgrid(t, t, indexing="ij")
    a, b, c = _hess_wilson(T, P, 1.0)
    d = np.sqrt(((a - c) / 2) ** 2 + b * b)
    assert np.max((a + c) / 2 + d) <= 4 + 1e-9 and np.min((a + c) / 2 - d) >= -4 - 1e-9
    a0, b0, c0 = _hess_wilson(np.pi, np.pi, 1.0)
    assert abs(a0 + 4) < 1e-12 and abs(c0 + 4) < 1e-12 and abs(b0) < 1e-12      # the minimum -4 is attained at (pi, pi)


# ------------------------------------------------------------------------------------------ Langevin Higgs toy and the three-variable model

def test_langevin_probabilities_and_zero_noise_limit():
    x0 = np.array([0.5, -0.6, 1.2, 0.1]); y0 = np.array([0.4, 0.2, -1.1, 0.9])
    p = Y.higgs_langevin(x0, y0, T=1e-3, gamma=1.0, M=8, seed=1)
    assert p.shape == (4, 5) and np.allclose(p.sum(axis=-1), 1.0)
    det = Y.higgs_basins(x0, y0, gamma=1.0)
    p0 = Y.higgs_langevin(x0, y0, T=0.0, gamma=1.0, M=1)
    assert np.all(np.argmax(p0, axis=-1) == det)                      # T = 0 reproduces the deterministic final well
    inner = np.array([0, 1, 3])                                          # the point (1.2, -1.1) is near a basin boundary
    assert np.all(p[inner, det[inner]] > 0.9)                          # strong friction, small noise: same well


def test_langevin_square_symmetry():
    x0 = np.array([0.5, 0.3]); y0 = np.array([0.2, -0.4])
    a = Y.higgs_langevin(x0, y0, T=0.0, gamma=0.3, M=1)
    b = Y.higgs_langevin(-x0, y0, T=0.0, gamma=0.3, M=1)               # x -> -x exchanges wells 1 and 3
    sw = np.array([0, 3, 2, 1, 4])
    assert np.all(np.argmax(a, axis=-1) == sw[np.argmax(b, axis=-1)])


def test_three_variable_reduces_to_two_and_is_symmetric():
    n = 30
    base, ys, ps = Y.section_picture(n=n, R=2.0, dt=0.004, t_max=60.0)
    Yg, P = np.meshgrid(ys, ps, indexing="ij")
    lab3 = Y.exit_basins3(Yg, P, 0.0, 0.0, R=2.0, dt=0.004, t_max=60.0)
    ok = base > 0
    assert np.mean(lab3[ok] == base[ok]) > 0.999                       # z = pz = 0 is invariant
    assert np.all(lab3[~ok] == -1)
    a = Y.exit_basins3(Yg, P, 0.05, 0.02, R=2.0, dt=0.004, t_max=40.0)
    b = Y.exit_basins3(Yg, P, -0.05, -0.02, R=2.0, dt=0.004, t_max=40.0)
    swz = {0: 0, -1: -1, 1: 1, 2: 2, 3: 3, 4: 4, 5: 6, 6: 5}
    assert np.mean(np.vectorize(swz.get)(a) == b) > 0.999            # (z, pz) -> -(z, pz) exchanges +z and -z


def test_exit_labels3_matches_exit_basins3_and_is_hexagon_equivariant():
    """General 3D exit function = exit_basins3 on x = 0 starts; the plane x+y+z=0 map is equivariant under the 60-degree rotation (labels permute along HEX_ORDER)."""
    rng = np.random.default_rng(1)
    y = rng.uniform(-1, 1, 60); z = rng.uniform(-0.2, 0.2, 60); py = rng.uniform(-0.5, 0.5, 60); pz = rng.uniform(-0.2, 0.2, 60)
    px = np.sqrt(2.0 - py ** 2 - pz ** 2 - 0.5 * y * y * z * z)
    a = Y.exit_basins3(y, py, z, pz, R=2.0, dt=0.004, t_max=60)
    b = Y.exit_labels3(np.stack([np.zeros(60), y, z], 1), np.stack([px, py, pz], 1), R=2.0, dt=0.004, t_max=60)
    assert np.mean(a == b) > 0.99
    assert sorted(Y.HEX_ORDER) == [1, 2, 3, 4, 5, 6]
    assert np.allclose(sorted(Y.HEX_ANGLES.values()), [30, 90, 150, 210, 270, 330])
    th = rng.uniform(0, 2 * np.pi, 200)
    l0 = Y.exit_plane_rest(np.cos(th), np.sin(th), R=1.15, dt=0.005, t_max=100)
    l1 = Y.exit_plane_rest(np.cos(th + np.pi / 3), np.sin(th + np.pi / 3), R=1.15, dt=0.005, t_max=100)
    nxt = {Y.HEX_ORDER[i]: Y.HEX_ORDER[(i + 1) % 6] for i in range(6)}
    dec = (l0 > 0) & (l1 > 0)
    assert dec.mean() > 0.95 and np.mean(np.array([nxt[v] for v in l0[dec]]) == l1[dec]) > 0.98
