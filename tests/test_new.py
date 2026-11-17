import numpy as np
import pytest
import lakes_of_wada as low
from lakes_of_wada import wadatest as W, maverick as M, entire


def sectors(n=200, k=3):
    yy, xx = np.mgrid[0:n, 0:n]
    ang = np.arctan2(yy - n / 2, xx - n / 2)
    return np.clip(((ang + np.pi) / (2 * np.pi) * k).astype(int), 0, k - 1) + 1


def test_newton_has_three_root_basins_and_is_wada():
    pic = low.entire_picture("newton", n=240, radius=1.5, degree=3)
    big = W.basin_labels(pic.label, min_fraction=0.01)
    assert len(big) == 3
    prof = W.wada_profile(pic.label, (4, 8))
    assert prof["K"] == 3
    assert prof["profile"][4] > 0.95 and prof["profile"][8] > 0.99


def test_ordinary_partition_is_not_wada():
    prof = W.wada_profile(sectors(), (2, 4))
    assert prof["profile"][2] < 0.1 and prof["profile"][4] < 0.2


def test_newton_degree_k_basins():
    for k in (3, 4, 5):
        pic = low.entire_picture("newton", n=160, radius=1.5, degree=k)
        # ignore the sliver of pixels that hit a pole exactly (grid artefact)
        assert len(W.basin_labels(pic.label, min_fraction=0.01)) == k


def test_exp_has_escape_and_attractor():
    pic = low.entire_picture("exp", n=120, radius=4, lam=0.3)
    assert {1, 2} <= set(np.unique(pic.label))


def test_chordal_finite_for_huge_values():
    d = entire.chordal(np.array([1e300 + 0j]), np.array([-1e300 + 0j]))
    assert np.isfinite(d).all() and d[0] <= 2


def test_maverick_zero_inside_basins_large_on_chaotic_pixels():
    pic = low.torus_picture(100)
    fld, basins = M.orbit_divergence(pic, n_iter=40, tail=15)
    s = M.maverick_score(fld)
    # pixels that did converge shadow their own basin's interior orbit
    # (only pixels that settle well inside the iteration window; longer orbits
    # are still wandering at step 40)
    decided = (pic.label > 0) & ~np.isnan(s) & (pic.steps <= 20)
    assert decided.any() and np.nanmax(s[decided]) < 1e-2
    # undecided (never-captured) pixels do not shadow any basin
    und = (pic.label == 0) & ~np.isnan(s)
    if und.any():
        assert np.nanmedian(s[und]) > 0.1


def test_maverick_requires_dynamics():
    pic = low.construction.wada_canals(n=64, stages=3)
    with pytest.raises(ValueError):
        M.orbit_divergence(pic)


def test_canal_construction_invariants():
    from scipy import ndimage as ndi
    c = low.construction.wada_canals(n=96, stages=6)
    assert set(np.unique(c.label)) == {0, 1, 2, 3}
    # lakes never touch each other directly (land keeps them apart)
    for a in (1, 2, 3):
        d = ndi.distance_transform_edt(c.label != a)
        for b in (1, 2, 3):
            if a != b:
                assert d[c.label == b].min() >= 2
    # every lake actually digs into the island
    assert min(np.sum(c.label == k) for k in (1, 2, 3)) > 1.5 * np.sum(low.construction.wada_canals(n=96, stages=0).label == 3)


def test_canal_construction_is_wada_at_its_resolution():
    """At radius n/4 nearly every land pixel sees all three lakes, far more than
    an ordinary 3-sector partition does at the same radius."""
    n = 192
    c = low.construction.wada_canals(n=n)
    r = n / 4
    got = W.wada_profile(c.label, (r,))["profile"][r]
    control = W.wada_profile(sectors(n), (r,))["profile"][r]
    assert got > 0.9 and got > control + 0.4


# --- n-dimensional Wada tests, merging method, scattering ----------------------

def test_wada_tools_work_in_3d():
    nw = low.entire_picture("newton", n=120, radius=1.5, degree=3).label
    vol = np.repeat(nw[None], 12, axis=0)                # extrusion of a Wada set is Wada
    assert W.wada_profile(vol, (4,))["profile"][4] > 0.95
    sec3 = np.repeat(sectors(120)[None], 12, axis=0)
    assert W.wada_profile(sec3, (4,))["profile"][4] < 0.1


def test_invalid_points_are_ignored():
    lab = low.entire_picture("newton", n=120, radius=1.5, degree=3).label.copy()
    lab[:, :30] = -1                                     # blank a strip: not a boundary
    b = W.boundary_mask(lab)
    assert not b[:, :30].any()
    assert W.wada_profile(lab, (4,))["profile"][4] > 0.9   # edge pixels cannot see removed basins


def test_merging_method_separates_wada_from_ordinary():
    nw = low.entire_picture("newton", n=200, radius=1.5, degree=3).label
    good = W.merging_test(nw)
    bad = W.merging_test(sectors(200))
    assert good["max_d_normalized"] < 0.03
    assert bad["max_d_normalized"] > 0.2


def test_basin_profile_flags_partial_wada():
    # basins 1 and 2 share a long straight boundary; basin 3 touches it only at one point
    n = 120
    lab = np.ones((n, n), int)
    lab[:, n // 2:] = 2
    lab[:12, n // 2 - 6: n // 2 + 6] = 3
    prof = W.basin_profile(lab, (3,), min_fraction=1e-3)
    assert prof[1][3] < 0.2 and prof[2][3] < 0.2


def test_hh_integrator_conserves_energy_and_has_no_exits_below_escape_energy():
    S = low.scattering
    sysm = S.henon_heiles()
    y0, valid = S.henon_heiles_section(6, 0.12)
    y = y0[valid][:30]
    lab, _, yf = S.integrate_exits(sysm, y, dt=0.01, t_max=20)
    assert (lab == 0).all()
    assert np.abs(sysm.energy(yf) - sysm.energy(y)).max() < 1e-7


def test_hh_matches_scipy_reference():
    from scipy.integrate import solve_ivp
    S = low.scattering
    sysm, probe = S.henon_heiles(), S.henon_heiles(4.99)
    y0, valid = S.henon_heiles_section(30, 0.25)
    pts = y0[valid]
    lab, st, _ = S.integrate_exits(sysm, pts, dt=0.02, t_max=30)
    quick = np.flatnonzero((lab > 0) & (st * 0.02 < 8))[::6][:25]
    assert len(quick) > 10
    for i in quick:
        ev = lambda t, s: s[0] ** 2 + s[1] ** 2 - 25
        ev.terminal = True
        r = solve_ivp(lambda t, s: sysm.rhs(s[None])[0], (0, 30), pts[i], method="DOP853",
                      rtol=1e-10, atol=1e-11, events=ev)
        assert int(probe.exit_label(r.y[:, -1][None])[0]) == lab[i]


def test_hh_exit_basins_are_wada_2d_and_3d():
    S = low.scattering
    p2 = S.scattering_picture(n=240, energy=0.25)
    assert set(np.unique(p2.label)) >= {-1, 1, 2, 3}
    assert W.wada_profile(p2.label, (4,))["profile"][4] > 0.97
    assert W.merging_test(p2.label)["max_d_normalized"] < 0.03
    p3 = S.scattering_picture(n=40, energy=0.25, dim=3)
    assert p3.label.shape == (40, 40, 40)
    assert W.wada_profile(p3.label, (3,))["profile"][3] > 0.9


def test_scattering_rejects_closed_energy():
    with pytest.raises(ValueError):
        low.scattering.scattering_picture(n=10, energy=0.1)


# --- lake labelling of the torus / plane pictures --------------------------------

def test_torus_and_plane_labels_are_the_lakes_not_the_first_hit_point():
    """(0,1/2), (1/2,1/2), (1/2,0) form a period-3 cycle, so the lakes are the
    cyclic colour classes (basins of the cycle points under f^3), not the index of
    the point an orbit happens to hit first (which is dust at every scale)."""
    t = low.torus_picture(160)
    assert np.array_equal(t.label, t.color.astype(int))
    p = low.plane_picture(160)
    assert np.array_equal(p.label, p.color.astype(int))


def test_torus_lakes_are_wada_and_not_dust():
    t = low.torus_picture(300)
    prof = W.wada_profile(t.label, (4,))
    assert prof["K"] == 4 and prof["profile"][4] > 0.98
    assert W.merging_test(t.label)["max_d_normalized"] < 0.03
    # each lake has genuinely solid interior (not interleaved at pixel scale)
    from scipy import ndimage as ndi
    for L in (1, 2, 3, 4):
        assert ndi.binary_erosion(t.label == L, iterations=3).sum() > 500


def test_maverick_high_scores_are_long_transients():
    p = low.torus_picture(200)
    fld, _ = M.orbit_divergence(p, n_iter=60, tail=20)
    s = M.maverick_score(fld)
    hi = M.maverick_mask(s, 0.1)
    # a pixel can only look 'maverick' if its orbit had not settled before the window
    assert (p.steps[hi] > 60 - 20).all()


# --- energy sweep -----------------------------------------------------------------

def test_energy_window_contains_the_allowed_region():
    S = low.scattering
    for E in (0.2, 1.0, 5.0):
        (y0, y1), (p0, p1) = S.energy_window(E)
        assert y0 < 0 < y1 and p0 == -p1
        # the left turning point sits exactly on the energy surface
        assert abs(S._hh_potential(0.0, y0) - E) < 1e-8


def test_uncertainty_exponent_distinguishes_fractal_from_smooth():
    S = low.scattering
    fractal = S.uncertainty_fraction(0.25, n_samples=3000)
    smoother = S.uncertainty_fraction(4.0, n_samples=3000)
    assert fractal["alpha"] < 0.4 < smoother["alpha"] < 1.0
    assert abs(fractal["dimension"] - (2 - fractal["alpha"])) < 1e-12
    # f(eps) grows with eps for a boundary of positive dimension
    f = fractal["f"]
    assert f[1e-4] < f[1e-2]


def test_energy_sweep_columns_and_trend_near_threshold():
    rows, pics = low.sweep.hh_energy_sweep([0.19, 0.3], n=60, alpha_samples=1500,
                                           t_max=60, keep=(0.3,))
    assert [r["energy"] for r in rows] == [0.19, 0.3]
    assert {"undecided_frac", "wada_r3", "min_basin_r3", "merging_pct_diag", "alpha"} <= set(rows[0])
    # just above the escape energy far more trajectories are still trapped
    assert rows[0]["undecided_frac"] > rows[1]["undecided_frac"]
    assert 0.3 in pics and pics[0.3].label.shape == (60, 60)


def test_sweep_csv_roundtrip(tmp_path):
    rows, _ = low.sweep.hh_energy_sweep([0.5], n=40, alpha_samples=800, t_max=40)
    path = tmp_path / "s.csv"
    low.sweep.write_csv(rows, path)
    import csv
    got = list(csv.DictReader(open(path)))
    assert float(got[0]["energy"]) == 0.5 and "alpha" in got[0]


# --- harmonic-polynomial family ---------------------------------------------------

def test_harmonic_degree3_is_henon_heiles():
    sc = low.scattering
    h, g = sc.henon_heiles(), sc.harmonic_polynomial(3)
    rng = np.random.default_rng(3)
    s = rng.normal(size=(500, 4)) * 3
    assert np.abs(h.rhs(s) - g.rhs(s)).max() < 1e-10
    assert np.abs(h.energy(s) - g.energy(s)).max() < 1e-10
    assert (h.exit_label(s * 3) == g.exit_label(s * 3)).all()
    assert sc.harmonic_saddle_energy(3) == pytest.approx(sc.HH_ESCAPE_ENERGY)


@pytest.mark.parametrize("degree", [4, 5, 6, 8])
def test_harmonic_channels_sit_at_the_saddles(degree):
    sc = low.scattering
    sysm = sc.harmonic_polynomial(degree)
    for k in range(1, degree + 1):
        th = (-np.pi / 2 + 2 * np.pi * k) / degree
        # the saddle of height 1/2 - 1/n is at r = 1 along theta_k, and the potential
        # falls off beyond it
        x, y = np.cos(th), np.sin(th)
        v = lambda r: sc._harm_potential(degree, r * x, r * y)
        assert v(1.0) == pytest.approx(sc.harmonic_saddle_energy(degree))
        assert v(1.2) < v(1.0) and v(0.8) < v(1.0)
        assert sysm.exit_label(np.array([[6 * x, 6 * y, 0, 0]]))[0] == k
    assert sysm.exit_label(np.array([[0.1, 0.1, 0, 0]]))[0] == 0   # inside


@pytest.mark.parametrize("degree", [4, 5])
def test_harmonic_conserves_energy_and_uses_all_channels(degree):
    sc = low.scattering
    E = 1.5 * sc.harmonic_saddle_energy(degree)
    y0, valid = sc.harmonic_section(degree, 60, E)
    sysm = sc.harmonic_polynomial(degree)
    s = y0[valid]
    assert np.abs(sysm.energy(s) - E).max() < 1e-12
    lab, _, _ = sc.integrate_exits(sysm, s, dt=0.02, t_max=60)
    assert set(np.unique(lab)) >= set(range(1, degree + 1))
    # RK4 energy drift over a short, still-bounded stretch of the motion (the final
    # plunge down a steep channel is not a fair test of a fixed step)
    lab2, _, fin = sc.integrate_exits(sysm, s, dt=0.01, t_max=1.0)
    inside = lab2 == 0
    assert inside.sum() > 100
    assert np.abs(sysm.energy(fin[inside]) - E).max() < 1e-5


def test_harmonic_picture_guards():
    sc = low.scattering
    with pytest.raises(ValueError):
        sc.scattering_picture("harmonic", n=10, energy=0.2, degree=4)   # below 1/4
    with pytest.raises(ValueError):
        sc.scattering_picture("harmonic", n=10, energy=0.5, degree=4, dim=3)
    with pytest.raises(ValueError):
        sc.harmonic_polynomial(2)


def test_degree_sweep_smoke():
    rows, _ = low.sweep.harmonic_degree_sweep([3, 4], n=90, alpha_samples=600, t_max=100)
    assert [r["degree"] for r in rows] == [3, 4]
    assert rows[0]["energy"] == pytest.approx(0.25)
    for r in rows:
        assert r["n_basins_seen"] == r["degree"]
        assert 0 < r["alpha"] <= 1.2
        assert 0.0 <= r["wada_r3_exits"] <= 1.0


# --- mixing band and the needed-radius diagnostic ---------------------------------

def test_wada_radius_separates_wada_from_non_wada_by_scaling():
    sl = low.sweep.radius_slope
    sizes = (120, 240, 480)
    sect = [W.wada_radius(sectors(n))["radius"][90] for n in sizes]
    assert sl(sizes, sect) > -0.1                  # stalls at a fixed physical length
    newt = [W.wada_radius(low.entire_picture("newton", n=n, radius=1.5, degree=3).label,
                          min_fraction=0.01)["radius"][90] for n in sizes]
    assert sl(sizes, newt) < -0.8                  # shrinks like 1/n


def test_wada_radius_matches_profile():
    L = low.entire_picture("newton", n=200, radius=1.5, degree=3).label
    r = W.wada_radius(L, (50,), min_fraction=0.01)["radius"][50]
    prof = W.wada_profile(L, (r,), min_fraction=0.01)["profile"][r]
    assert prof >= 0.5 - 1e-9


def test_mixing_band_is_a_subset_of_the_boundary():
    L = sectors(200)
    b = W.boundary_mask(L)
    band = W.mixing_band(L, 8, 0.25)
    assert (band & ~b).sum() == 0
    d = W.boundary_density(L, 8)
    assert d.min() >= 0 and d.max() <= 1
    out = W.band_wada_profile(L, (1, 3), 8, 0.25)
    assert 0.0 <= out["band_fraction"] <= 1.0
    assert out["boundary_pixels"] == int(band.sum())
    # the three-sector control has only a thin boundary: no dense mixing core
    dense = W.mixing_band(L, 8, 0.6)
    assert dense.sum() == 0


def test_harmonic_radius_scaling_smoke():
    rows = low.sweep.harmonic_radius_scaling([3], sizes=(80, 160, 320), t_max=60)
    assert [r["n"] for r in rows] == [80, 160, 320]
    assert rows[0]["slope_p90"] < -0.5
    assert 0 <= rows[0]["rest0.6_pixels"] <= rows[0]["boundary_pixels"]
