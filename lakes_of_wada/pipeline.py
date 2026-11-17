"""High-level entry points that glue maps, views and classifiers together."""
from dataclasses import dataclass
import numpy as np
from . import maps, classify, views, entire


@dataclass
class Picture:
    color: np.ndarray   # float, indexed [i, j] like the legacy scripts
    label: np.ndarray
    steps: np.ndarray
    params: dict
    # optional: lets analysis tools (maverick) re-iterate the dynamics
    f: object = None
    points: object = None
    distance: object = None


def torus_picture(n=100, amplitude=2.2, spread=0.23, islands=4, max_iter=500) -> Picture:
    f = maps.torus_map(amplitude, spread, islands)
    pts = views.torus_view(n)
    if islands == 1:
        res = classify.orbit_to_targets(f, pts, ((0.0, 0.0),), max_iter=max_iter)
        color = classify.color_angle(res)
    else:
        res = classify.orbit_to_targets(f, pts, max_iter=max_iter)
        color = classify.color_cyclic(res)
        # (0, 1/2), (1/2, 1/2), (1/2, 0) form a period-3 *cycle*, not three fixed
        # points: the lakes are the basins of the three cycle points under f^3, i.e.
        # the cyclic colour class, not "which point was hit first"
        res.label = color.astype(int)
    return Picture(color, res.label, res.steps,
                   dict(kind="torus", n=n, amplitude=amplitude, spread=spread, islands=islands),
                   f=f, points=pts, distance=maps.torus_distance)


def plane_picture(n=100, amplitude1=2.4, spread1=0.5, amplitude2=2.3, spread2=0.25,
                  radius=4.0, theta=1.7, max_iter=500) -> Picture:
    f = maps.plane_map(amplitude1, spread1, amplitude2, spread2)
    pts, zero = views.plane_view(n, radius, theta)
    res = classify.orbit_to_targets(f, pts, max_iter=max_iter)
    color = classify.color_cyclic(res)
    color[zero] = 1.0
    res.label = color.astype(int)       # lakes = cyclic colour classes (see torus_picture)
    return Picture(color, res.label, res.steps,
                   dict(kind="plane", n=n, amplitude1=amplitude1, spread1=spread1,
                        amplitude2=amplitude2, spread2=spread2, radius=radius, theta=theta),
                   f=f, points=pts, distance=maps.torus_distance)


def entire_picture(family="newton", n=400, center=0j, radius=2.0, max_iter=200,
                   escape_radius=1e4, **fparams) -> Picture:
    """Basin picture for a plane map from ``entire.FAMILIES`` (or pass a
    callable as ``family``)."""
    f = entire.FAMILIES[family](**fparams) if isinstance(family, str) else family
    z = views.complex_window(n, center, radius)
    res = entire.classify_attractors(
        f, z, max_iter=max_iter, escape_radius=escape_radius)
    lab = res.label
    color = np.where(lab > 0, lab + 0.9 * np.exp(-res.steps / 15.0), 0.0)
    name = family if isinstance(family, str) else getattr(family, "__name__", "custom")
    return Picture(color, lab, res.steps,
                   dict(kind="entire", family=name, n=n, center=center, radius=radius,
                        max_iter=max_iter, **fparams),
                   f=f, points=z, distance=entire.chordal)
