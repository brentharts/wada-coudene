"""The four Spectre window pieces Z_p(S) (S = O, G; p = 0, 1) rasterised from their exact graph-directed system (lakes_of_wada/spectre_window.py,
data figures/spectre_window_maps.npz from brentharts/spectre boundary_certificate/spectre_window.py).

Part 1 (calibration of the dimension estimators).  The boundary of each piece against the value quoted by the source: 'seam dimension'
log(5 + 2 sqrt 6) / log(4 + sqrt 15) = 1.1110 (a proved lower bound for the boundary dimension of the window; box counts of the source give 1.10 to 1.13).
  (a) box counts: number of boundary pixels of the in/out image at n = 100, 200, 400, 800, log-log slope;
  (b) the uncertainty-fraction form: fraction of the piece's area within eps of its boundary (n = 1600), alpha = log-log slope, D = 2 - alpha,
      fitted in three windows of eps (0.008-0.04, 0.016-0.08, 0.03-0.15) because the answer depends on the window.
Part 2 (Wada-type test).  The 13 children of Z_0(O) and Z_1(O), labelled by child; for each pixel on a seam between two children, the number of
different children within physical radius rho (n = 800 and 1600), against rho.
Usage: python3 figures/analyse_spectre_window.py [depth]      (default 7; about 10 minutes)"""
import csv, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import spectre_window as SW     # noqa: E402

DEPTH = int(sys.argv[1]) if len(sys.argv) > 1 else 7
KNOWN = np.log(5 + 2 * np.sqrt(6)) / np.log(4 + np.sqrt(15))
rows1, rows2 = [], []


def boundary_pixels(inside):
    pad = np.pad(inside, 1, constant_values=False)
    nb = pad[:-2, 1:-1] & pad[2:, 1:-1] & pad[1:-1, :-2] & pad[1:-1, 2:]
    return inside & ~nb


print(f"depth {DEPTH}; known seam dimension {KNOWN:.4f}")
for S, p in (("O", 0), ("O", 1), ("G", 0), ("G", 1)):
    t0 = time.time()
    pts, kid = SW.points(S, p, DEPTH)
    # (a) box counts, a common box for all n
    img, box, px = SW.label_image(pts, kid, 800)
    nb, ns = [], (100, 200, 400, 800)
    for n in ns:
        im, _, _ = SW.label_image(pts, kid, n, box=box)
        nb.append(int(boundary_pixels(im >= 0).sum()))
    sl = float(np.polyfit(np.log(ns), np.log(nb), 1)[0])
    loc = np.diff(np.log(nb)) / np.diff(np.log(ns))
    # (b) uncertainty fraction from the distance to the boundary, n = 1600
    im, box16, px16 = SW.label_image(pts, kid, 1600, box=box)
    ins = im >= 0
    d_in = ndi.distance_transform_edt(ins) * px16
    Dw = []
    for lo, hi in ((0.008, 0.04), (0.016, 0.08), (0.03, 0.15)):               # three fit windows, each 8 values of eps
        eps = np.geomspace(lo, hi, 8)
        f = np.array([np.mean(d_in[ins] < e) for e in eps])
        Dw.append(2 - float(np.polyfit(np.log(eps), np.log(f), 1)[0]))
    area = ins.sum() * px16 ** 2
    print(f"Z_{p}({S}): {len(pts)} points, area {area:.3f}; (a) boundary-pixel slope {sl:.3f} (local {np.round(loc, 3)}); "
          f"(b) D from the uncertainty fraction in three eps windows {np.round(Dw, 3)}   [{time.time() - t0:.0f}s]", flush=True)
    rows1.append((S, p, DEPTH, len(pts), area, sl, *Dw, KNOWN))
    # Part 2 for the sets with more than two children
    if S == "O":
        for n in (800, 1600):
            im, _, pxn = SW.label_image(pts, kid, n, box=box)
            labs = [k for k in range(13) if (im == k).any()]
            dist = np.stack([ndi.distance_transform_edt(im != k) * pxn for k in labs])        # distance to child k (0 inside it)
            pad = np.pad(im, 1, constant_values=-1)
            seam = (im >= 0) & ((pad[:-2, 1:-1] != im) & (pad[:-2, 1:-1] >= 0) | (pad[2:, 1:-1] != im) & (pad[2:, 1:-1] >= 0)
                                | (pad[1:-1, :-2] != im) & (pad[1:-1, :-2] >= 0) | (pad[1:-1, 2:] != im) & (pad[1:-1, 2:] >= 0))
            for rho in (0.004, 0.008, 0.016, 0.032, 0.064):
                if rho < 2 * pxn:
                    continue
                cnt = (dist[:, seam] <= rho).sum(axis=0)
                rows2.append((S, p, n, rho, int(seam.sum()), float(np.mean(cnt == 2)), float(np.mean(cnt == 3)), float(np.mean(cnt >= 4)), float(cnt.mean())))
                print(f"   Z_{p}(O) n={n} rho={rho}: {int(seam.sum())} seam pixels; exactly 2 children {np.mean(cnt == 2):.3f}, 3 {np.mean(cnt == 3):.3f}, "
                      f">=4 {np.mean(cnt >= 4):.3f}, mean {cnt.mean():.2f}", flush=True)
with open(HERE / "spectre_window_analysis.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["shape", "parity", "depth", "points", "area", "box_slope", "D_eps_0.008_0.04", "D_eps_0.016_0.08", "D_eps_0.03_0.15", "known_seam_dim"]); w.writerows(rows1)
    w.writerow([]); w.writerow(["shape", "parity", "n", "rho", "seam_pixels", "frac_2", "frac_3", "frac_ge4", "mean_children"]); w.writerows(rows2)
