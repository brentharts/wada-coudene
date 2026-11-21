"""Section 8 of YANGMILLS.md ('Langevin noise = displacement of the start, flip fraction ~ delta^alpha') applied to a set of known boundary dimension: the Spectre
window piece Z_0(O) (D = 1.1110 for its seams, source), displaced in internal space by a Gaussian of standard deviation delta per coordinate (a random
phason displacement).  f_exit(delta) = probability that a uniform point of the window leaves it; f_child(delta) = probability that its child label (13 children)
changes or it leaves.  Expected f ~ delta^(2 - D) = delta^0.889 if the boundary of the child partition has D = 1.1110 at these scales.
Usage: python3 figures/analyse_phason_noise.py [n] -> figures/phason_noise.csv"""
import sys, csv
sys.path.insert(0, ".")
import numpy as np
from scipy import ndimage as ndi
from lakes_of_wada import spectre_window as SW

n = int(sys.argv[1]) if len(sys.argv) > 1 else 2400
pts, kid = SW.points("O", 0, 7)
lo, hi = pts.min(0), pts.max(0); c = 0.5 * (lo + hi); h = 0.5 * float((hi - lo).max()) + 0.6
lab, box, px = SW.label_image(pts, kid, n, box=(c[0] - h, c[0] + h, c[1] - h, c[1] + h), reach=1.5)
inw = lab >= 0
area = inw.sum() * px * px
print(f"n={n}, pixel {px:.5f}, area of Z_0(O) {area:.4f}", flush=True)
D_SEAM = np.log(5 + 2 * np.sqrt(6)) / np.log(4 + np.sqrt(15))
deltas = np.array([0.006, 0.008, 0.012, 0.016, 0.024, 0.032, 0.048, 0.064, 0.096, 0.128, 0.192, 0.256])
fe, fc = [], []
for d in deltas:
    s = d / px
    stay = ndi.gaussian_filter(inw.astype(float), s, mode="constant", truncate=4.0)
    fe.append(float(1 - stay[inw].mean()))
    same = np.zeros(inw.shape)
    for k in range(13):
        m = (lab == k)
        same[m] = ndi.gaussian_filter(m.astype(float), s, mode="constant", truncate=4.0)[m]
    fc.append(float(1 - same[inw].mean()))
fe, fc = np.array(fe), np.array(fc)
print("delta      ", deltas.tolist()); print("f_exit     ", np.round(fe, 4).tolist()); print("f_child    ", np.round(fc, 4).tolist(), flush=True)
rows = [("delta", deltas.tolist()), ("f_exit", fe.tolist()), ("f_child", fc.tolist()), ("seam_dimension_source", float(D_SEAM)), ("2_minus_D", float(2 - D_SEAM))]
print(f"expected exponent 2 - D = {2 - D_SEAM:.4f}; in terms of T = delta^2: {(2 - D_SEAM) / 2:.4f}")
for name, f in (("exit", fe), ("child", fc)):
    loc = np.diff(np.log(f)) / np.diff(np.log(deltas))
    print(name, "local slopes:", np.round(loc, 3).tolist())
    for a, b in ((0, 5), (2, 7), (4, 9), (6, 11)):
        al = np.polyfit(np.log(deltas[a:b]), np.log(f[a:b]), 1)[0]
        print(f"   fit delta {deltas[a]:.3f}-{deltas[b-1]:.3f}: alpha = {al:.3f}  (D = {2 - al:.3f})")
        rows.append((f"alpha_{name}_{deltas[a]}_{deltas[b-1]}", float(al)))
with open("figures/phason_noise.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerows(rows)
