"""Exit map of the three-variable model H = p^2/2 + (x^2 y^2 + y^2 z^2 + z^2 x^2)/2 for orbits released from rest on the plane x + y + z = 0, the plane
perpendicular to the three-fold axis of the cube, in which the six exit half axes project to a regular hexagon (labels in the order Y.HEX_ORDER).
H is homogeneous of degree 4, so a start at (rho, theta) is a rescaled copy of the start at (1, theta): the exit axis at radius R rho depends on theta
and R only, and the map is a function of the angle theta.  This script computes it on the sector [30, 90] degrees (between the projections of +x and -z),
tests the 60-degree equivariance against the step-size floor, and measures (i) the uncertainty exponent of theta (fraction of labels changed by a shift
delta) and (ii) the six-channel Wada profile along theta.
Usage: python3 figures/analyse_hexagon_exit.py [n] [R] [t_max]  ->  figures/hexagon_exit_theta_R{R}_n{n}.npy, figures/hexagon_exit_analysis_R{R}.csv"""
import sys, csv, time
sys.path.insert(0, ".")
import numpy as np
from lakes_of_wada import yangmills as Y

n = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
R = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
TMAX = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
DT = 0.005
rows = []
def labels(theta, dt=DT):
    P = Y.plane_start(np.cos(theta), np.sin(theta))
    return Y.exit_labels3(P, np.zeros_like(P), R=R, dt=dt, t_max=TMAX)
# --- equivariance under the rotation by 60 degrees, against the step floor (same starts, dt and dt/2)
rng = np.random.default_rng(0)
th = rng.uniform(0, 2 * np.pi, 3000)
l0, l0h, l1 = labels(th), labels(th, DT / 2), labels(th + np.pi / 3)
nxt = {Y.HEX_ORDER[i]: Y.HEX_ORDER[(i + 1) % 6] for i in range(6)}
dec = (l0 > 0) & (l1 > 0)
eq = float(np.mean(np.array([nxt[a] for a in l0[dec]]) == l1[dec]))
d2 = (l0 > 0) & (l0h > 0)
floor = float(np.mean(l0[d2] == l0h[d2]))
print(f"R={R}: rotation by 60 deg permutes the labels along {Y.HEX_ORDER}: agree on {eq:.4f} of {int(dec.sum())} decided pairs; step floor (dt vs dt/2) agree on {floor:.4f}", flush=True)
rows.append(("equivariance_60deg_agree", eq, int(dec.sum()))); rows.append(("step_floor_agree", floor, int(d2.sum())))
# --- the map on the sector
t0 = time.time()
theta = np.radians(np.linspace(30, 90, n))
lab = labels(theta)
np.save(f"figures/hexagon_exit_theta_R{R}_n{n}.npy", lab)
print(f"{n} angles in {time.time() - t0:.0f}s; undecided {np.mean(lab == 0):.4f}; label fractions {[round(float(np.mean(lab == L)), 4) for L in range(1, 7)]}", flush=True)
rows.append(("undecided", float(np.mean(lab == 0)), n))
# --- uncertainty exponent: fraction of decided samples whose label differs from that at a shift of s samples (delta = 60 deg * s / (n - 1))
ok = lab > 0
dl, f = [], []
for s in (1, 2, 4, 8, 16, 32, 64, 128, 256):
    a, b = lab[:-s], lab[s:]
    m = (a > 0) & (b > 0)
    dl.append(np.radians(60.0) * s / (n - 1)); f.append(float(np.mean(a[m] != b[m])))
dl, f = np.array(dl), np.array(f)
print("delta (rad):", np.round(dl, 6).tolist()); print("f(delta):   ", np.round(f, 4).tolist(), flush=True)
loc = np.diff(np.log(f)) / np.diff(np.log(dl))
print("local slopes alpha between successive deltas:", np.round(loc, 3).tolist(), flush=True)
for lo, hi in ((0, 5), (2, 7), (4, 9)):
    al = np.polyfit(np.log(dl[lo:hi]), np.log(f[lo:hi]), 1)[0]
    print(f"  fit over delta {dl[lo]:.5f} to {dl[hi-1]:.5f}: alpha = {al:.3f}  (D = 1 - alpha = {1 - al:.3f})", flush=True)
    rows.append((f"alpha_fit_{lo}_{hi}", float(al), f"{dl[lo]:.6f}-{dl[hi-1]:.6f}"))
rows.append(("f_of_delta", f.tolist(), dl.tolist()))
# --- six-channel Wada profile along theta: boundary = sample whose neighbour has another decided label; channels within r samples
bnd = np.where((lab[:-1] > 0) & (lab[1:] > 0) & (lab[:-1] != lab[1:]))[0]
print("boundary samples:", bnd.size, "of", n, flush=True)
for r in (1, 2, 4, 8, 16, 32, 64):
    seen = np.array([len(set(lab[max(0, i - r):i + r + 2][lab[max(0, i - r):i + r + 2] > 0])) for i in bnd])
    h = [float(np.mean(seen == k)) for k in range(1, 7)]
    print(f"  r={r:3d} samples ({np.degrees(60.0 / (n - 1)) * 0 + 60.0 * r / (n - 1):.4f} deg): channels seen 1..6:", [round(x, 3) for x in h], " all six:", round(h[5], 3), flush=True)
    rows.append((f"channels_within_{r}", h, bnd.size))
with open(f"figures/hexagon_exit_analysis_R{R}.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["quantity", "value", "extra"]); w.writerows(rows)
