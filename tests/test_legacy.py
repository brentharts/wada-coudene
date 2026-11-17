"""Compare the vectorized module against scalar copies of the legacy logic."""
import numpy as np
import pytest
from scipy.special import elliprf
import lakes_of_wada as low
from lakes_of_wada.elliptic import inv_p

phi = 1.6180339887498948482
A = np.array([[2, 1], [1, 1]])
B = (1.0 / (1 + phi * phi)) * np.array([[phi * phi, phi], [phi, 1]])
wrap = lambda p: ((p + 0.5) % 1) - 0.5


def legacy_torus_color(p0, amp, spr, islands):
    k = lambda t: 0.0 if abs(t) > 1 else (1 - t * t) ** 2
    ctr = (lambda p: p - np.floor(p + 0.5)) if islands == 1 else \
          (lambda p: p - 0.5 * np.floor(2.0 * p + 0.5))
    f = lambda p: wrap(A.dot(p) - amp * k(ctr(p)[0] / spr) * k(ctr(p)[1] / spr) * B.dot(ctr(p)))
    eps = 1e-4

    def near(p):
        x, y = p
        if abs(x) < eps:
            return 1 if abs(y) < eps else (2 if abs(abs(y) - .5) < eps else 0)
        if abs(abs(x) - .5) < eps:
            return 3 if abs(y) < eps else (4 if abs(abs(y) - .5) < eps else 0)
        return 0
    p, N, nr = p0, 0, 0
    while not nr:
        if N > 500:
            break
        p = f(p); N += 1
        nr = near(p) if islands == 4 else (1 if abs(p[0]) < eps and abs(p[1]) < eps else 0)
    if islands == 1:
        x, y = p
        return (0 if x == 0 and y == 0 else 0.5 * (1 + x / np.sqrt(x * x + y * y))), N
    return (nr if nr <= 1 else (nr + N) % 3 + 2), N


@pytest.mark.parametrize("islands", [1, 4])
def test_torus_matches_legacy(islands):
    """The map is chaotic (error grows ~2.6x per step), so float round-off makes
    long orbits diverge from the scalar code; short orbits must agree exactly."""
    n = 24
    pic = low.torus_picture(n, islands=islands)
    pts = low.views.torus_view(n)
    checked = 0
    for i in range(n):
        for j in range(n):
            col, N = legacy_torus_color(pts[i, j], 2.2, 0.23, islands)
            if N <= 25:
                assert pic.steps[i, j] == N
                if islands == 4:   # (1-island colour is the angle of a sub-eps residual)
                    assert pic.color[i, j] == col
                checked += 1
    assert checked > 100


def test_inv_p_values():
    # reference values from mpmath: invper * ellipf(asin(z**-0.5), -1)
    assert abs(inv_p(2.0) - 0.27724255490843436) < 1e-12
    assert abs(inv_p(1.0) - 0.5) < 1e-12
    # real z in (0,1) is a branch cut: resolved as Im z -> 0+
    assert abs(inv_p(0.5) - (0.5 - 0.22275744509156561j)) < 1e-12
    assert abs(inv_p(0.5) - inv_p(0.5 + 1e-9j)) < 1e-6
    assert np.all(np.isfinite(inv_p(np.array([1e-3, -3 + 1j, 4j]))))


def test_plane_runs():
    pic = low.plane_picture(20)
    assert pic.color.shape == (20, 20)
    assert set(np.unique(pic.label)) <= {0, 1, 2, 3, 4}
