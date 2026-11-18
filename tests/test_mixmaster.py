import numpy as np
from scipy.integrate import solve_ivp

from lakes_of_wada import mixmaster as M, wadatest as W


def test_constraint_is_conserved_by_the_field_equations():
    # I = sum_{i<j} x_i' x_j' - V is conserved: it is the Hamiltonian constraint (found numerically, see module doc)
    rng = np.random.default_rng(0)
    for _ in range(3):
        y0 = np.r_[rng.normal(size=3) * 0.3, rng.normal(size=3)]
        rhs = lambda t, y: np.r_[y[3:], M.force(y[:3])]
        sol = solve_ivp(rhs, [0, 3], y0, rtol=1e-11, atol=1e-12, dense_output=True)
        vals = [M.constraint(sol.sol(t)[:3], sol.sol(t)[3:]) for t in np.linspace(0, 3, 7)]
        assert np.ptp(vals) < 1e-8


def test_turnaround_states_satisfy_the_constraint():
    bp, bm = np.array([0.1, -0.2, 0.3]), np.array([0.05, 0.1, -0.1])
    x, v, V = M.turnaround_state(bp, bm)
    assert np.all(V < 0) and np.allclose(v.sum(axis=-1), 0, atol=1e-14)
    assert np.allclose(M.constraint(x, v), 0, atol=1e-12)


def test_three_fold_symmetry_of_the_exit_basins():
    rng = np.random.default_rng(3)
    pts = rng.uniform(-0.8, 0.8, size=(300, 2))
    x, _, V = M.turnaround_state(pts[:, 0], pts[:, 1])
    pts, x = pts[V < 0], x[V < 0]
    l0 = M.exit_basins(pts[:, 0], pts[:, 1])
    bp, bm = M.x_to_beta(np.roll(x, 1, axis=-1))        # cyclic permutation of the axes
    l1 = M.exit_basins(bp, bm)
    perm = {1: 2, 2: 3, 3: 1}
    m = (l0 > 0) & (l1 > 0)
    assert m.sum() > 50 and all(perm[a] == b for a, b in zip(l0[m], l1[m]))


def test_biaxial_line_exits_through_its_own_axis_and_stays_biaxial():
    # beta_- = 0 is b = c (a Berger sphere, axis a = the Hopf fibre): the Taub solutions exit through axis a
    bp = np.linspace(-0.4, 1.0, 25)
    lab = M.exit_basins(bp, np.zeros_like(bp))
    assert set(lab[lab != -1].tolist()) == {1}


def test_exit_basins_are_three_and_wada_at_pixel_scale():
    pic = M.mixmaster_picture(n=80, half=1.2)
    prof = W.wada_profile(pic.label, (4, 8, 16))
    assert prof["K"] == 3
    assert prof["profile"][16] > 0.99
    # the same picture with other exit thresholds
    for kw in (dict(p_exit=0.9), dict(p_exit=0.99, nmax=15000)):
        p = W.wada_profile(M.mixmaster_picture(n=80, half=1.2, **kw).label, (8, 16))
        assert p["K"] == 3 and p["profile"][16] > 0.98


def test_time_and_epoch_outputs_do_not_change_the_labels():
    bp = np.array([0.1, -0.2, 0.3, 2.0, 5.0])
    bm = np.array([0.05, 0.1, -0.1, 0.0, 0.0])
    lab = M.exit_basins(bp, bm)
    lab2, drift, tau, ne = M.exit_basins(bp, bm, return_info=True, return_time=True, return_epochs=True)
    assert np.array_equal(lab, lab2) and drift < 1e-3
    ok = lab > 0
    assert np.all(tau[ok] > 0) and np.all(ne[ok] >= 1)          # an exit is a Kasner epoch, so at least one was entered
    assert np.all(np.isnan(tau[lab < 0])) and np.all(ne[lab < 0] == -1)


def test_labels_inside_the_smooth_wedges_do_not_depend_on_the_step():
    bp = np.array([-0.4, 0.9, 0.7, 0.3, 0.5])
    bm = np.array([0.0, 0.0, 0.0, -0.35, 0.0])
    for h0 in (0.01, 0.005):
        assert M.exit_basins(bp, bm, h0=h0).tolist() == [1, 1, 1, 3, 1]
