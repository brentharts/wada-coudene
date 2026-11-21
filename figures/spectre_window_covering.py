"""Covering of internal space by the windows of the tiles that overlap a given Spectre tile, and the Wada-type test across window boundaries.

Setting.  A tile of window class a (shape S_a, frame A_a) placed at lattice vector u is present in a tiling iff the internal point y = STAR u lies in
its window Omega_a = c_a + sc R_a Z_0(S_a) (Z_0 from lakes_of_wada.spectre_window; the affine fit of the classes against a level-N patch is that of
brentharts/spectre boundary_certificate/window_boundary_closure.py).  A tile (b, u + w) overlapping the placed tile (a, u) is present iff
y in Omega_b - STAR w, i.e. in the frame z of Omega_a:   z_b = R_b^T R_a z + R_b^T (c_a - c_b + STAR w)/sc  must lie in Z_0(S_b).
Lattice constraint (the point that the first attempt, spectre_window_neighbours.py, missed): the tile vectors of one class lie in a coset of a sublattice
of Z^4 of index 36 (Smith form 1,1,6,6), so w must be congruent to (v_b - u_a) modulo that sublattice, v_b and u_a being tile vectors of the patch.
Outputs: pixel statistics in the z frame of Z_0(S_a): inside Z_0 the number of overlapping tiles present must be 0 (a tiling has no overlaps); outside,
the number is >= 1 (the plane is covered).  Multiplicity of distinct neighbour windows at the boundary is the Wada-type statistic.
Usage: SPECTRE_DIR=/path/to/spectre python3 figures/spectre_window_covering.py N class [npix]    -> figures/spectre_window_covering_N{N}_c{class}.csv
"""
import io, contextlib, os, sys, itertools, time
SP = os.environ.get("SPECTRE_DIR", "../spectre")
sys.path.insert(0, os.path.join(SP, "boundary_certificate")); sys.path.insert(0, SP)
import numpy as np
import mpmath as mp
t00 = time.time()
with contextlib.redirect_stdout(io.StringIO()):
    import spectre_window as sw
from spectre_exact import slot_transforms, SPECTRE, mul, apply, IDENT, SUPER_RULES, ttrans, trot, to_c
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
print("imports %.0fs" % (time.time() - t00), flush=True)
lam = float(sw.LAM); RAD = float(sw.RAD)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
A_CLASS = int(sys.argv[2]) if len(sys.argv) > 2 else 0
NPIX = int(sys.argv[3]) if len(sys.argv) > 3 else 700
quads = [[SPECTRE[3], SPECTRE[5], SPECTRE[7], SPECTRE[11]]]; TS = [None]
for j in range(N):
    Ts, nq = slot_transforms(quads[-1]); TS.append(Ts); quads.append(nq)
Bm = mp.matrix(4, 4)
for c, v in enumerate(sw.E[0] + sw.F[0]):
    for r in range(4): Bm[r, c] = v[r]
Bim = mp.inverse(Bm)
STARm = sw.T[0] * mp.matrix([[Bim[r, c] for c in range(4)] for r in (2, 3)])
STAR = np.array(STARm.tolist(), dtype=float)
groups = {}
VEC = {}
def rec(s, j, T, z):
    if j == 0:
        g = groups.setdefault((tuple(T[0].flatten()), 'G' if s == 'Gamma' else 'O'), [T[0], [], []])
        g[1].append(STARm * mp.matrix([mp.mpf(int(x)) for x in T[1]])); g[2].append(z)
        VEC.setdefault((tuple(T[0].flatten()), 'G' if s == 'Gamma' else 'O'), []).append(tuple(int(x) for x in T[1])); return
    p = (j - 1) % 2
    for i, c in enumerate(SUPER_RULES[s]):
        if c is not None:
            rec(c, j - 1, mul(T, TS[j][i]), sw.O[(i, p)].T * (sw.GA[(i, p)] + z / sw.LAM))
rec('Delta', N, IDENT, mp.matrix(2, 1))
FR, FRm = [], []
for (fk, S), (A, Y, Z) in groups.items():
    X = mp.matrix([[z[0], z[1], 1] for z in Z]); Ym = mp.matrix([[y[0], y[1]] for y in Y])
    NX, B2 = X.T * X, X.T * Ym
    cols = [mp.lu_solve(NX, B2[:, k]) for k in range(2)]
    sol = mp.matrix([[cols[k][i] for k in range(2)] for i in range(3)])
    R, c = sol[0:2, 0:2].T, sol[2, 0:2].T
    res = max(abs(x) for x in (X * sol - Ym)); G = R.T * R; s2 = G[0, 0]
    assert res < mp.mpf(10) ** -40 and max(abs(x) for x in (G - s2 * mp.eye(2))) < mp.mpf(10) ** -40, (res, G)
    sc = mp.sqrt(s2); Rm = R / sc
    FRm.append((Rm, c, sc))
    FR.append((S, A, np.array(Rm.tolist(), dtype=float), np.array(c.tolist(), dtype=float).ravel(), float(sc), len(Z)))
scale = FR[0][4]
print('window classes (shape, frame):', len(FR), ' tiles per class', [f[5] for f in FR], ' shapes', [f[0] for f in FR], ' similarity ratio', scale, ' %.0fs' % (time.time() - t00), flush=True)
LEAVES = {'O': [IDENT], 'G': [IDENT, mul(ttrans(SPECTRE[8]), trot(1))]}
def shape(S, A, t):
    return unary_union([Polygon([(lambda w: (w.real, w.imag))(to_c(apply(mul((A, t), L), v))) for v in SPECTRE])
                        for L in LEAVES[S]])

PH = np.array([[to_c(np.eye(4, dtype=object)[k]).real for k in range(4)], [to_c(np.eye(4, dtype=object)[k]).imag for k in range(4)]])
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lakes_of_wada import spectre_window as SWR
from scipy.spatial import cKDTree
import sympy
from sympy.matrices.normalforms import hermite_normal_form
keys = list(groups.keys())
# lattice of tile-vector differences (index 36 in Z^4), common to all classes
D = [tuple(a - b for a, b in zip(v, VEC[keys[0]][0])) for v in VEC[keys[0]][1:]]
H = hermite_normal_form(sympy.Matrix(D).T)
assert H.shape == (4, 4)
detH = abs(int(H.det())); ADJ = np.array((H.inv() * detH).tolist(), dtype=object).astype(np.int64)
print("difference lattice index in Z^4:", detH, flush=True)
def in_lattice(W):                                   # rows w; w in lattice iff H^-1 w integral
    return np.all((W @ ADJ.T) % detH == 0, axis=1)
for k in keys:                                       # every class has differences in the same lattice
    d = np.array([[a - b for a, b in zip(v, VEC[k][0])] for v in VEC[k]], dtype=np.int64)
    assert in_lattice(d).all()
trees = {S: cKDTree(SWR.points(S, 0, 7)[0]) for S in "OG"}
TOL = float(os.environ.get("TOL", "0.008"))
SUF = "" if NPIX == 700 else "_n%d" % NPIX
a = A_CLASS
S1, A1, R1, c1, _, _ = FR[a]; u_a = np.array(VEC[keys[a]][0], dtype=np.int64)
P1 = shape(S1, A1, np.zeros(4, dtype=object))
# candidate relative vectors: small physical part, bounded internal part, in the right coset for each class b
ymax = 2 * RAD * scale + 2 * max(np.linalg.norm(f[3]) for f in FR) + 1
M4 = np.vstack([PH, STAR]); Mi = np.linalg.inv(M4)
bnd = int(np.ceil(np.linalg.norm(Mi, axis=1).max() * np.hypot(12.0, ymax))) + 1
rng = np.arange(-bnd, bnd + 1)
cands = []
for x0 in rng:
    Vc = np.array(list(itertools.product([x0], rng, rng, rng)))
    ok = np.linalg.norm(Vc @ PH.T, axis=1) < 12.0
    Vc = Vc[ok]
    if len(Vc): cands.append(Vc[np.linalg.norm(Vc @ STAR.T, axis=1) < ymax])
cands = np.concatenate(cands); print("lattice candidates", len(cands), "box edge", bnd, flush=True)
nbrs = []                                            # (b, w, B, d)
for b in range(len(FR)):
    Sb, Ab, Rb, cb, _, _ = FR[b]
    v_b = np.array(VEC[keys[b]][0], dtype=np.int64)
    W = cands[in_lattice(cands - (v_b - u_a))]
    dd = ((c1 - cb)[None, :] + W @ STAR.T) @ Rb / scale          # R_b^T (c_a - c_b + STAR w)/sc as rows
    W, dd = W[np.linalg.norm(dd, axis=1) <= 2 * RAD + 1e-9], dd[np.linalg.norm(dd, axis=1) <= 2 * RAD + 1e-9]
    for w, d in zip(W, dd):
        if (w == 0).all() and b == a: continue
        if shape(Sb, Ab, np.array([int(x) for x in w], dtype=object)).intersection(P1).area < 1e-6: continue
        nbrs.append((b, Sb, Rb.T @ R1, d, w))
print("overlapping placements with window within 2 RAD, coset-restricted:", len(nbrs), flush=True)
from shapely.geometry import Point
q = P1.representative_point()
# the placements that contain the probe point q of T_a: their windows (with that of T_a) partition internal space, each point once
cover = [k for k, (b, Sb, B, d, w) in enumerate(nbrs) if shape(Sb, FR[b][1], np.array([int(x) for x in w], dtype=object)).contains(q)]
print("probe point (%.3f, %.3f); placements containing it: %d" % (q.x, q.y, len(cover)), flush=True)
ALL = os.environ.get("ALL", "0") == "1"
# pixel grid in the z frame of Z_0(S_a)
Zc = trees[S1].data
lo, hi = Zc.min(0) - 0.45, Zc.max(0) + 0.45
xs = np.linspace(lo[0], hi[0], NPIX); ys = np.linspace(lo[1], hi[1], int(NPIX * (hi[1] - lo[1]) / (hi[0] - lo[0])))
GX, GY = np.meshgrid(xs, ys); G = np.stack([GX.ravel(), GY.ravel()], 1)
own = trees[S1].query(G)[0] < TOL

cnt = np.zeros(len(G), int)
MASK = {}
for k, (b, Sb, B, d, w) in enumerate(nbrs):                  # z_b = B z + d with B = R_b^T R_a
    if not ALL and k not in cover: continue
    lo_b, hi_b = trees[Sb].data.min(0) - 2 * TOL, trees[Sb].data.max(0) + 2 * TOL
    Zb = G @ B.T + d
    cand_px = np.where(np.all((Zb >= lo_b) & (Zb <= hi_b), axis=1))[0]
    m = np.zeros(len(G), bool)
    if len(cand_px): m[cand_px] = trees[Sb].query(Zb[cand_px])[0] < TOL
    cnt += m
    if k in cover: MASK[k] = m
from scipy.ndimage import distance_transform_edt
img = own.reshape(GX.shape); px = xs[1] - xs[0]
d_in = distance_transform_edt(img) * px                         # depth inside Z_0
d_out = distance_transform_edt(~img) * px                       # distance outside
print("pixels:", len(G), " in Z_0:", own.sum(), " area %.3f" % (own.sum() * px * (ys[1] - ys[0])), " pixel %.4f, tolerance %.3f" % (px, TOL))
np.savez("figures/spectre_window_covering_N%d_c%d%s.npz" % (N, a, SUF), own=img, cnt=cnt.reshape(GX.shape), xs=xs, ys=ys, nnbrs=len(nbrs))
import csv
rows = []
def stat(lab, mask):
    c = cnt[mask]; rows.append((lab, int(mask.sum()), *(float(np.mean(c == k)) for k in range(3)), float(np.mean(c >= 3))))
    print("%-26s n=%7d  P(0)=%.3f P(1)=%.3f P(2)=%.3f P(>=3)=%.3f" % ((lab, mask.sum()) + rows[-1][2:]), flush=True)
din, dout = d_in.ravel(), d_out.ravel()
for lo_, hi_ in ((0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 9)):
    stat("inside, depth %.2f-%.2f" % (lo_, hi_ if hi_ < 9 else 9), own & (din >= lo_) & (din < hi_))
for lo_, hi_ in ((0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 9)):
    stat("outside, dist %.2f-%.2f" % (lo_, hi_ if hi_ < 9 else 9), ~own & (dout >= lo_) & (dout < hi_))
with open("figures/spectre_window_covering_N%d_c%d%s.csv" % (N, a, SUF), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["subset", "pixels", "P0", "P1", "P2", "P3plus"]); w.writerows(rows)

# ---- partition of internal space by the windows of the placements containing q (own window included) ----
from scipy.ndimage import binary_dilation
pieces = [own] + [MASK[k] for k in cover]
tot = sum(m.astype(int) for m in pieces)
far = ~own & (dout > 0.02) | own & (din > 0.02)
# the partition property is only meaningful away from the thin tolerance collar of each window's boundary
print("PARTITION: pixels with exactly one piece: %.4f, none: %.4f, two or more: %.4f (all pixels); outside the collar of the own window: %.4f / %.4f / %.4f"
      % (np.mean(tot == 1), np.mean(tot == 0), np.mean(tot >= 2), np.mean(tot[far] == 1), np.mean(tot[far] == 0), np.mean(tot[far] >= 2)), flush=True)
shp = GX.shape
band = (dout < 0.012) & ~own                      # pixels just outside Z_0(S_a)
band_all = (dout < 0.012)
rows2 = []
for r in (0.008, 0.016, 0.032, 0.064):
    it = max(1, int(round(r / px)))
    dil = [binary_dilation(m.reshape(shp), iterations=it).ravel() for m in pieces]
    npieces = sum(d.astype(int) for d in dil)
    # boundary of Z_0(S_a): pixels within one pixel-tolerance of it, counted over the r-dilated pieces
    bd = (din < 0.012) & own | (dout < 0.012) & ~own
    c = npieces[bd]
    rows2.append((r, int(bd.sum()), *(float(np.mean(c == k)) for k in (1, 2, 3)), float(np.mean(c >= 4)), float(np.mean(c >= 3))))
    print("r=%.3f: boundary pixels %d: pieces within r: 1: %.3f  2: %.3f  3: %.3f  >=4: %.3f   (>=3: %.3f)" % rows2[-1], flush=True)
with open("figures/spectre_window_partition_N%d_c%d%s.csv" % (N, a, SUF), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["r", "boundary_pixels", "P1", "P2", "P3", "P4plus", "P3plus"]); w.writerows(rows2)
    w.writerow(["probe_pieces", len(pieces)])
