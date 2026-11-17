"""Parameter sweeps: how do the Wada diagnostics change with the energy?

For each energy of the Henon-Heiles system (x = 0 section, window chosen from the
energy) we record

* undecided fraction (trajectories not out by ``t_max``) and exit fractions,
* the Wada profile at ``r`` pixels (fraction of boundary pixels that see all
  three exits) and the *smallest per-basin* value (a partially Wada basin
  drags this down),
* the merging-method distance (percentage of the diagonal),
* the uncertainty exponent alpha from perturbing initial conditions
  (alpha ~ 1: smooth boundary; alpha << 1: fractal),
* optionally the same profile numbers at a finer resolution ``n_fine``.  A
  boundary that is Wada keeps (or improves) its fraction as the grid is refined;
  one that is not (basins meeting only along thin curves) loses it.
"""
import csv
import time
import numpy as np

from . import scattering, wadatest


def hh_energy_sweep(energies, n=500, n_fine=None, r=3, t_max=None, alpha_samples=20000,
                    dt=0.05, keep=(), progress=None, fine_min_energy=0.25):
    """Return ``(rows, pictures)``: a list of dicts (one per energy) and a dict
    ``{energy: Picture}`` for the energies listed in ``keep``."""
    rows, pics = [], {}
    for E in energies:
        t0 = time.time()
        tm = t_max if t_max is not None else (200.0 if E < 0.3 else 100.0)
        yr, pr = scattering.energy_window(E)

        def run(nn):
            return scattering.scattering_picture(n=nn, energy=E, dt=dt, t_max=tm,
                                                 y_range=yr, py_range=pr)

        p = run(n)
        L = p.label
        valid = L >= 0
        row = {"energy": E, "n": n, "t_max": tm,
               "undecided_frac": float(np.mean(L[valid] == 0))}
        for k in (1, 2, 3):
            row[f"exit{k}_frac"] = float(np.mean(L[valid] == k))
        if len(wadatest.basin_labels(L)) == 3:
            prof = wadatest.wada_profile(L, (r, 2 * r))["profile"]
            bp = wadatest.basin_profile(L, (r,))
            m = wadatest.merging_test(L)
            row.update({f"wada_r{r}": prof[r], f"wada_r{2 * r}": prof[2 * r],
                        f"min_basin_r{r}": min(v[r] for v in bp.values()),
                        "merging_pct_diag": 100 * m["max_d_normalized"]})
        else:
            row.update({f"wada_r{r}": np.nan, f"wada_r{2 * r}": np.nan,
                        f"min_basin_r{r}": np.nan, "merging_pct_diag": np.nan})
        if E in keep:
            pics[E] = p
        if n_fine and E >= fine_min_energy:
            Lf = run(n_fine).label
            if len(wadatest.basin_labels(Lf)) == 3:
                bpf = wadatest.basin_profile(Lf, (r,))
                row[f"wada_r{r}_fine"] = wadatest.wada_profile(Lf, (r,))["profile"][r]
                row[f"min_basin_r{r}_fine"] = min(v[r] for v in bpf.values())
            else:
                row[f"wada_r{r}_fine"] = row[f"min_basin_r{r}_fine"] = np.nan
            row["n_fine"] = n_fine
        elif n_fine:
            row[f"wada_r{r}_fine"] = row[f"min_basin_r{r}_fine"] = np.nan
            row["n_fine"] = n_fine
        u = scattering.uncertainty_fraction(E, n_samples=alpha_samples, dt=dt, t_max=tm)
        row["alpha"], row["dimension"] = u["alpha"], u["dimension"]
        row["seconds"] = time.time() - t0
        rows.append(row)
        if progress:
            progress(row)
    return rows, pics


def write_csv(rows, path):
    keys = list(rows[0])
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()})
