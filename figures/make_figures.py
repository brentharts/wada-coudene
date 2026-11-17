#!/usr/bin/env python3
"""Regenerate every figure in figures/ (and the table behind the profile chart).

    python3 figures/make_figures.py            # all figures
    python3 figures/make_figures.py torus newton   # just some

Figures are written next to this script.  Sizes/parameters are fixed here so the
README's numbers can be reproduced.
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import lakes_of_wada as low                      # noqa: E402
from lakes_of_wada import wadatest as W          # noqa: E402

RADII = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96)

# categorical slots 1-5 of the validated reference palette (light mode), fixed order
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
SERIES = {  # key: (label, colour, marker)
    "newton": ("Newton map of z³ − 1 (plane)", "#2a78d6", "o"),
    "henon": ("Hénon–Heiles exit basins (E = 0.25)", "#eb6834", "s"),
    "plane": ("Plane picture (℘ lift of perturbed cat map)", "#1baf7a", "^"),
    "torus": ("Torus picture, 4 islands", "#e87ba4", "D"),
    "canal": ("Canal construction (n = 512)", "#eda100", "v"),
    "sectors": ("Control: ordinary 3-sector partition", "#8a8984", "x"),
}

_cache = {}


def picture(key):
    """Compute (and memoise) the pictures used by several figures."""
    if key in _cache:
        return _cache[key]
    if key == "torus":
        p = low.torus_picture(700)
    elif key == "plane":
        p = low.plane_picture(800)
    elif key == "newton":
        p = low.entire_picture("newton", n=800, radius=1.5, degree=3)
    elif key == "henon":
        p = low.scattering.scattering_picture(n=600, energy=0.25)
    elif key == "canal":
        p = low.construction.wada_canals(n=512)
    elif key == "sectors":
        n = 512
        yy, xx = np.mgrid[0:n, 0:n]
        ang = np.arctan2(yy - n / 2, xx - n / 2)
        lab = np.clip(((ang + np.pi) / (2 * np.pi) * 3).astype(int), 0, 2) + 1
        p = low.Picture(lab.astype(float), lab, np.zeros_like(lab), dict(kind="control"))
    _cache[key] = p
    return p


def fig_image(key, cmap, fname, crop=False):
    p = picture(key)
    if crop:   # trim the empty outer region to the picture's bounding box
        occupied = p.label > (1 if p.params.get("kind") == "plane" else 0)
        rows, cols = np.flatnonzero(occupied.any(axis=1)), np.flatnonzero(occupied.any(axis=0))
        pad = 24
        r0, r1 = max(rows[0] - pad, 0), rows[-1] + pad
        c0, c1 = max(cols[0] - pad, 0), cols[-1] + pad
        p = low.Picture(p.color[r0:r1, c0:c1], p.label[r0:r1, c0:c1],
                        p.steps[r0:r1, c0:c1], p.params)
    low.save_image(p, HERE / fname, cmap)
    print("wrote", fname)


def fig_maverick():
    p = picture("torus")
    fld, _ = low.maverick.orbit_divergence(p, n_iter=60, tail=20)
    score = low.maverick.maverick_score(fld)
    fig, axs = plt.subplots(1, 2, figsize=(10, 5.2), dpi=130, facecolor="white")
    img = np.ma.masked_where(p.label <= 0, p.color)
    cm = plt.get_cmap("plasma").copy(); cm.set_bad("white")
    axs[0].imshow(img, cmap=cm, interpolation="nearest")
    axs[0].set_title("The four lakes (basin of the origin + three cycle-point basins)", color=INK, fontsize=11)
    cm2 = plt.get_cmap("magma_r").copy(); cm2.set_bad("white")
    im = axs[1].imshow(np.ma.masked_invalid(score), cmap=cm2, interpolation="nearest",
                       vmin=0, vmax=0.71)
    axs[1].set_title("Orbit divergence (dark: orbit still not shadowing its lake)",
                     color=INK, fontsize=11)
    for a in axs:
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, ax=axs[1], fraction=0.046, pad=0.02)
    cb.set_label("max torus distance over the last 20 of 60 iterates", color=INK2, fontsize=9)
    cb.ax.tick_params(colors=INK2, labelsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "maverick.png", facecolor="white")
    plt.close(fig)
    print("wrote maverick.png")


def fig_profiles():
    rows, curves = [], {}
    for key, (label, _, _) in SERIES.items():
        p = picture(key)
        n = p.label.shape[-1]
        prof = W.wada_profile(p.label, RADII)["profile"]
        curves[key] = (n, prof)
        for r in RADII:
            rows.append(dict(series=label, n=n, radius_px=r,
                             radius_pct_of_width=round(100 * r / n, 3),
                             fraction_seeing_all_basins=round(prof[r], 4)))
    mrows = []
    for key, (label, _, _) in SERIES.items():
        m = W.merging_test(picture(key).label)
        mrows.append(dict(series=label, n=picture(key).label.shape[-1], K=m["K"],
                          max_hausdorff_px=round(m["max_d"], 2), min_hausdorff_px=round(m["min_d"], 2),
                          max_over_diagonal_pct=round(100 * m["max_d_normalized"], 3)))
    with open(HERE / "wada_merging.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(mrows[0]))
        w.writeheader(); w.writerows(mrows)
    with open(HERE / "wada_profiles.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    fig, ax = plt.subplots(figsize=(8.6, 6.3), dpi=150, facecolor="white")
    ax.set_facecolor("white")
    for key, (label, col, mk) in SERIES.items():
        n, prof = curves[key]
        x = [100 * r / n for r in RADII]
        y = [prof[r] for r in RADII]
        ls = (0, (4, 3)) if key == "sectors" else "-"
        ax.plot(x, y, color=col, lw=2, ls=ls, marker=mk, ms=5.5, mfc=col if key != "sectors" else col,
                mec="white", mew=0.8 if key != "sectors" else 1.4, label=label, clip_on=False)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 20)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([0.1, 0.3, 1, 3, 10])
    ax.set_xticklabels(["0.1%", "0.3%", "1%", "3%", "10%"])
    ax.set_xlabel("neighbourhood radius r, as a percentage of the image width", color=INK2, fontsize=10)
    ax.set_ylabel("fraction of boundary pixels that see all basins within r", color=INK2, fontsize=10)
    ax.set_title("Wada profiles: how quickly do all basins meet at every boundary point?",
                 color=INK, fontsize=12, loc="left")
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=INK2, labelsize=9)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False,
              fontsize=9, labelcolor=INK, columnspacing=1.6, handlelength=2.6)
    fig.tight_layout()
    fig.savefig(HERE / "wada_profiles.png", facecolor="white")
    plt.close(fig)
    print("wrote wada_profiles.png + wada_profiles.csv")
    return curves


JOBS = {
    "torus": lambda: fig_image("torus", "plasma", "torus.png"),
    "plane": lambda: fig_image("plane", "turbo", "plane.png", crop=True),
    "newton": lambda: fig_image("newton", "viridis", "newton.png"),
    "henon": lambda: fig_image("henon", "Set2", "henon_heiles.png"),
    "canal": lambda: fig_image("canal", "Set2", "canals.png"),
    "maverick": fig_maverick,
    "profiles": fig_profiles,
}

if __name__ == "__main__":
    for name in (sys.argv[1:] or JOBS):
        JOBS[name]()
