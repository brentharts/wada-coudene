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


def harmonic_degree_sweep(degrees, ratio=1.5, n=500, n_fine=None, r=3, t_max=100.0,
                          alpha_samples=20000, dt=0.05, keep=(), progress=None):
    """Wada diagnostics of the harmonic-polynomial family V = r^2/2 + Im(z^n)/n
    (``degree`` = n exit channels; n = 3 is Henon-Heiles) at a *common relative
    energy* ``E = ratio * E_saddle(n)``, where ``E_saddle = 1/2 - 1/n`` is the
    energy at which the channels open.  ``ratio = 1.5`` gives E = 0.25 for n = 3.

    Row keys: degree, energy, undecided_frac, share_min/share_max (smallest and
    largest exit basin), wada_r{r}, wada_r{2r}, min_basin_r{r} (all computed over
    all ``degree`` basins, nan if some basin is missing from the section),
    n_basins_seen, the ``*_exits`` variants (trapped points ignored instead of
    counted as boundary), optional ``*_exits_fine`` values at ``n_fine``, alpha,
    dimension.
    """
    rows, pics = [], {}
    for d in degrees:
        t0 = time.time()
        E = ratio * scattering.harmonic_saddle_energy(d)

        def run(nn):
            return scattering.scattering_picture("harmonic", n=nn, energy=E, dt=dt, t_max=t_max,
                                                 degree=d)

        def wada_numbers(L, radii):
            if len(wadatest.basin_labels(L)) != d:
                return None
            prof = wadatest.wada_profile(L, radii)["profile"]
            bp = wadatest.basin_profile(L, (radii[0],))
            return prof, min(v[radii[0]] for v in bp.values())

        p = run(n)
        L = p.label
        valid = L >= 0
        shares = [float(np.mean(L[valid] == k)) for k in range(1, d + 1)]
        row = {"degree": d, "energy": E, "n": n, "t_max": t_max,
               "undecided_frac": float(np.mean(L[valid] == 0)),
               "share_min": min(shares), "share_max": max(shares),
               "n_basins_seen": len(wadatest.basin_labels(L))}
        w = wada_numbers(L, (r, 2 * r))
        row[f"wada_r{r}"] = w[0][r] if w else np.nan
        row[f"wada_r{2 * r}"] = w[0][2 * r] if w else np.nan
        row[f"min_basin_r{r}"] = w[1] if w else np.nan
        # Same numbers with the trapped (regular-island) points ignored, so the
        # question is only "are the boundaries *between exit basins* Wada?"
        Lx = L.copy()
        Lx[Lx == 0] = -1
        wx = wada_numbers(Lx, (r, 2 * r))
        row[f"wada_r{r}_exits"] = wx[0][r] if wx else np.nan
        row[f"min_basin_r{r}_exits"] = wx[1] if wx else np.nan
        if n_fine:
            Lf = run(n_fine).label
            Lfx = Lf.copy()
            Lfx[Lfx == 0] = -1
            wf = wada_numbers(Lfx, (r,))
            row[f"wada_r{r}_exits_fine"] = wf[0][r] if wf else np.nan
            row[f"min_basin_r{r}_exits_fine"] = wf[1] if wf else np.nan
            row["n_fine"] = n_fine
        u = scattering.uncertainty_fraction(E, n_samples=alpha_samples, dt=dt, t_max=t_max,
                                            degree=d)
        row["alpha"], row["dimension"] = u["alpha"], u["dimension"]
        if d in keep:
            pics[d] = p
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
