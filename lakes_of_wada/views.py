"""Views: how pixels map to points on the torus."""
import numpy as np
from .elliptic import inv_p


def torus_view(n):
    """Legacy torus sampling.  Returns points of shape (n, n, 2) indexed [i, j]
    with x = -1/2 + i/(n-1), y = 1/2 - j/(n-1)."""
    i = np.arange(n) / (n - 1)
    x = -0.5 + i[:, None] + 0 * i[None, :]
    y = 0.5 - i[None, :] + 0 * i[:, None]
    return np.stack([x, y], axis=-1)


def plane_view(n, radius=4.0, theta=1.7):
    """Legacy plane sampling: pixels are points z of the complex plane
    (window of half-size ``radius``, rotated by ``theta``), lifted to the torus
    through the inverse Weierstrass P function.  Returns (points[n,n,2],
    z[n,n]); z == 0 is flagged in ``zero``.
    """
    s = -1 + 2 * np.arange(n) / (n - 1)
    z = radius * np.exp(1j * theta) * (s[:, None] + 1j * s[None, :])
    zero = z == 0
    w = inv_p(np.where(zero, 1.0, z))
    pts = np.stack([w.real, w.imag], axis=-1)
    return pts, zero


def complex_window(n, center=0j, radius=2.0):
    """Square window of the complex plane, standard image orientation:
    ``z[row, col]`` with x increasing along columns and y increasing upward."""
    s = -1 + 2 * np.arange(n) / (n - 1)
    return center + radius * (s[None, :] + 1j * (-s)[:, None])
