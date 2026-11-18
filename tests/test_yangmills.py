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
