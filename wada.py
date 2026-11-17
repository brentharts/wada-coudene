#!/usr/bin/env python3
"""Draw dynamical Lakes of Wada (method of Coudene) using lakes_of_wada.

Examples:
  wada.py torus -n 400 -o torus.png
  wada.py plane -n 400 --data plane.npz        # also/only save raw arrays
  wada.py plane -n 400 -o - --data plane.npz   # data only, no image
  wada.py entire --family newton --degree 3 -n 600 --wada      # Wada profile
  wada.py torus -n 300 --maverick mav.png                      # orbit-divergence map
  wada.py canal -n 256 --wada                                  # canal construction
  wada.py scatter -n 500 --energy 0.25 --wada --merging        # Henon-Heiles exit basins (2D section)
  wada.py scatter -n 64 --dim 3 --data hh3d.npz --wada         # full 3D energy surface
"""
import argparse
import sys
import time

import numpy as np

import lakes_of_wada as low


def add_output_args(p):
    p.add_argument("-n", "--size", type=int, default=100, help="image size in pixels")
    p.add_argument("--max-iter", type=int, default=500, help="iteration budget per point")
    p.add_argument("-c", "--colormap", help="matplotlib colormap name")
    p.add_argument("-o", "--output", help="image filename ('-' for no image; "
                   "default: parameter/time based name)")
    p.add_argument("--data", metavar="FILE.npz",
                   help="save raw arrays (color, label, steps, params) as .npz")
    p.add_argument("--wada", nargs="?", const="1,2,4,8,16", metavar="RADII",
                   help="print the Wada profile (fraction of boundary pixels seeing all "
                        "basins within r pixels); optional comma-separated radii")
    p.add_argument("--merging", action="store_true",
                   help="also run the merging-method Wada test (Hausdorff distances)")
    p.add_argument("--maverick", metavar="FILE.png",
                   help="save a maverick-point orbit-divergence map (needs dynamics)")
    p.add_argument("--maverick-iters", type=int, default=60)
    p.add_argument("--maverick-tail", type=int, default=20)


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("torus", help="picture on the 2-torus")
    t.add_argument("-a", "--amplitude", type=float, default=2.2)
    t.add_argument("-s", "--spread", type=float, default=0.23)
    t.add_argument("-i", "--islands", type=int, choices=(1, 4), default=4)
    add_output_args(t)

    p = sub.add_parser("plane", help="picture in the plane (via inverse Weierstrass P)")
    p.add_argument("--amplitude1", type=float, default=2.4)
    p.add_argument("--spread1", type=float, default=0.5)
    p.add_argument("--amplitude2", type=float, default=2.3)
    p.add_argument("--spread2", type=float, default=0.25)
    p.add_argument("-r", "--radius", type=float, default=4.0)
    p.add_argument("-t", "--theta", type=float, default=1.7, help="rotation (radians)")
    add_output_args(p)

    e = sub.add_parser("entire", help="basins of a plane map (Newton, exp, sin, polynomial)")
    e.add_argument("--family", choices=sorted(low.entire.FAMILIES), default="newton")
    e.add_argument("--degree", type=int, default=3, help="newton: p(z)=z^degree-1")
    e.add_argument("--coeffs", type=complex, nargs="+", help="newton/poly coefficients, lowest first")
    e.add_argument("--lam", type=float, default=None, help="exp/sin parameter")
    e.add_argument("--center", type=complex, default=0j)
    e.add_argument("-r", "--radius", type=float, default=2.0)
    e.add_argument("--escape-radius", type=float, default=1e4)
    add_output_args(e)

    sc = sub.add_parser("scatter", help="exit basins of an open Hamiltonian (Henon-Heiles)")
    sc.add_argument("--system", choices=("henon-heiles",), default="henon-heiles")
    sc.add_argument("--energy", type=float, default=0.25, help="energy (> 1/6 for exits)")
    sc.add_argument("--dim", type=int, choices=(2, 3), default=2,
                    help="2: section x=0 in (y, p_y); 3: whole energy surface")
    sc.add_argument("--dt", type=float, default=0.05)
    sc.add_argument("--t-max", type=float, default=60.0)
    add_output_args(sc)

    c = sub.add_parser("canal", help="Wada canal construction on a raster (no dynamics)")
    c.add_argument("--lakes", type=int, default=3)
    c.add_argument("--stages", type=int, default=15)
    c.add_argument("--layout", choices=("ponds", "sectors"), default="ponds")
    add_output_args(c)
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.cmd == "torus":
        pic = low.torus_picture(a.size, a.amplitude, a.spread, a.islands, a.max_iter)
        default_cmap = "plasma"
        stem = f"wada-torus-i{a.islands}-n{a.size}"
    elif a.cmd == "entire":
        kw = {}
        if a.family == "newton":
            kw = dict(coeffs=a.coeffs, degree=a.degree)
        elif a.family == "poly":
            kw = dict(coeffs=a.coeffs)
        elif a.lam is not None:
            kw = dict(lam=a.lam)
        pic = low.entire_picture(a.family, a.size, a.center, a.radius, a.max_iter,
                                 a.escape_radius, **kw)
        default_cmap = "viridis"
        stem = f"wada-entire-{a.family}-n{a.size}"
    elif a.cmd == "scatter":
        pic = low.scattering.scattering_picture(a.system, a.size, a.energy, a.dim,
                                                a.dt, a.t_max)
        default_cmap = "Set2"
        stem = f"wada-scatter-{a.system}-E{a.energy}-d{a.dim}-n{a.size}"
    elif a.cmd == "canal":
        from lakes_of_wada.construction import wada_canals
        pic = wada_canals(a.size, a.lakes, a.stages, layout=a.layout)
        default_cmap = "Set2"
        stem = f"wada-canal-n{a.size}"
    else:
        pic = low.plane_picture(a.size, a.amplitude1, a.spread1, a.amplitude2,
                                a.spread2, a.radius, a.theta, a.max_iter)
        default_cmap = "magma"
        stem = f"wada-plane-n{a.size}-r{a.radius}-t{a.theta}"

    if a.wada:
        radii = tuple(float(x) for x in a.wada.split(","))
        res = low.wadatest.wada_profile(pic.label, radii)
        print(f"Wada profile: K={res['K']} basins, {res['boundary_pixels']} boundary pixels")
        for r, v in res["profile"].items():
            print(f"  r = {r:g} px: {v:.3f} of boundary pixels see all basins")
    if a.merging:
        m = low.wadatest.merging_test(pic.label)
        print(f"Merging test: K={m['K']}, max Hausdorff distance {m['max_d']:.2f} px "
              f"(min {m['min_d']:.2f}); {100 * m['max_d_normalized']:.2f}% of the array diagonal")
    if a.maverick and pic.f is None:
        print("--maverick needs a map with attached dynamics; skipped for this picture")
    elif a.maverick:
        fld, basins = low.maverick.orbit_divergence(pic, n_iter=a.maverick_iters,
                                                    tail=a.maverick_tail)
        score = low.maverick.maverick_score(fld)
        low.save_scalar(score, a.maverick, "magma")
        print(f"Saved maverick map: {a.maverick} (basins {basins}; "
              f"{low.maverick.maverick_mask(score).sum()} candidate pixels)")
    if a.data:
        np.savez_compressed(a.data, color=pic.color, label=pic.label,
                            steps=pic.steps, params=np.array(repr(pic.params)))
        print(f"Saved data: {a.data}")
    if a.output != "-":
        out = a.output or f"{stem}-{int(time.time())}.png"
        shown = pic
        if pic.label.ndim == 3:            # 3D volume: save the central slice
            k = pic.label.shape[0] // 2
            shown = low.Picture(pic.color[k], pic.label[k], pic.steps[k], pic.params)
            print(f"3D volume: image is the central slice (index {k})")
        low.save_image(shown, out, a.colormap or default_cmap)
        print(f"Saved image: {out}")


if __name__ == "__main__":
    sys.exit(main())
