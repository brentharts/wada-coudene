"""Classify points by where their orbits go (basins of attracting fixed points)."""
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple
import numpy as np
from .maps import wrap

TORSION4 = ((0.0, 0.0), (0.0, 0.5), (0.5, 0.0), (0.5, 0.5))


@dataclass
class Result:
    label: np.ndarray    # int: 0 = undecided, k = reached targets[k-1]
    steps: np.ndarray    # iterations performed
    final: np.ndarray    # final positions (..., 2)


def orbit_to_targets(f, p0, targets: Sequence[Tuple[float, float]] = TORSION4,
                     eps=1e-4, max_iter=500) -> Result:
    """Iterate ``f`` on every point until it lands within ``eps`` (per
    coordinate, periodically) of one of ``targets``, or the budget runs out.
    Matches the legacy loop: up to ``max_iter + 1`` applications of f.
    """
    p = np.array(p0, dtype=float)
    shape = p.shape[:-1]
    p = p.reshape(-1, 2)
    n = len(p)
    label = np.zeros(n, dtype=int)
    steps = np.zeros(n, dtype=int)
    active = np.arange(n)
    tg = np.asarray(targets, dtype=float)
    for k in range(max_iter + 1):
        if active.size == 0:
            break
        p[active] = f(p[active])
        steps[active] = k + 1
        d = np.abs(wrap(p[active][:, None, :] - tg[None, :, :])) < eps  # (m, T, 2)
        hit = d.all(axis=-1)
        any_hit = hit.any(axis=1)
        label[active[any_hit]] = hit[any_hit].argmax(axis=1) + 1
        active = active[~any_hit]
    return Result(label.reshape(shape), steps.reshape(shape), p.reshape(shape + (2,)))


# --- colouring rules: Result -> float array -----------------------------------

def color_cyclic(res: Result):
    """Legacy: basin 0/1 get their own value; other basins cycle with orbit
    length mod 3 (this is what makes the Wada boundaries visible)."""
    lab = res.label
    return np.where(lab <= 1, lab, (lab + res.steps) % 3 + 2).astype(float)


def color_angle(res: Result):
    """Legacy one-island torus picture: colour by the direction of the final
    position, 0.5 * (1 + cos(angle))."""
    x, y = res.final[..., 0], res.final[..., 1]
    r = np.hypot(x, y)
    return np.where(r == 0, 0.0, 0.5 * (1 + x / np.where(r == 0, 1, r)))
