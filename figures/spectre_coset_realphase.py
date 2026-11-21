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
PBOX = float(os.environ.get('PBOX','12'))
bnd = int(np.ceil(np.linalg.norm(Mi, axis=1).max() * np.hypot(PBOX, ymax))) + 1
rng = np.arange(-bnd, bnd + 1)
cands = []
for x0 in rng:
    Vc = np.array(list(itertools.product([x0], rng, rng, rng)))
    ok = np.linalg.norm(Vc @ PH.T, axis=1) < PBOX
    Vc = Vc[ok]
    if len(Vc): cands.append(Vc[np.linalg.norm(Vc @ STAR.T, axis=1) < ymax])
cands = np.concatenate(cands); print("lattice candidates", len(cands), "box edge", bnd, flush=True)

# ====================== relative coset of the tile covering a far point, as a function of the phase ======================
import json, collections
IN_T, OUT_T = float(os.environ.get("IN_T", "0.003")), float(os.environ.get("OUT_T", "0.010"))
from shapely.geometry import Point
QX, QY = float(os.environ.get("QX", "-6.4")), float(os.environ.get("QY", "0.3"))
q = Point(QX, QY)
NPH = int(os.environ.get("NPH", "3000"))
def resid(W):                                  # coset of Z^4 / H of each row
    return (np.atleast_2d(W) @ ADJ.T) % detH
Q = []                                         # all placements (any coset) whose footprint contains q
for b in range(len(FR)):
    Sb, Ab, Rb, cb, _, _ = FR[b]
    v_b = np.array(VEC[keys[b]][0], dtype=np.int64)
    base = resid(v_b - u_a)[0]
    for w in cands:
        if np.hypot(w @ PH[0] - QX, w @ PH[1] - QY) > float(os.environ.get('RADQ','7')): continue
        d = ((c1 - cb) + STAR @ w) @ Rb / scale
        if np.linalg.norm(d) > float(os.environ.get('DCUT',str(2*RAD+0.6))): continue
        poly = shape(Sb, Ab, np.array([int(t) for t in w], dtype=object))
        if not poly.contains(q): continue
        kappa = tuple(int(x) for x in (resid(w)[0] - base) % detH)
        Q.append((b, Sb, Rb.T @ R1, d, kappa, tuple(int(t) for t in w)))
print("placements containing q (any coset): %d, distinct relative cosets %d" % (len(Q), len(set(x[4] for x in Q))), flush=True)
cloud = trees[S1].data; lo, hi = cloud.min(0) - 0.3, cloud.max(0) + 0.3
rng_ = np.random.default_rng(21)
Zb = rng_.uniform(lo, hi, (NPH * 4, 2)); Zb = Zb[trees[S1].query(Zb)[0] < IN_T][:NPH]
# child label of Z_0(S1): nearest depth-7 point of the first-map index
pts7, kid7 = SWR.points(S1, 0, 7)
kids = kid7[trees[S1].query(Zb)[1]]
rows = []
for z, kd in zip(Zb, kids):
    present, amb = set(), False
    for (b, Sb, B, d, kappa, w) in Q:
        dist = trees[Sb].query(B @ z + d)[0]
        if dist < IN_T: present.add((b, kappa))
        elif dist <= OUT_T: amb = True
    if amb: continue
    allk = sorted(set(x[4] for x in Q)); pk = set(k for _, k in present)
    rows.append((float(z[0]), float(z[1]), int(kd), [list(k) for k in allk if k not in pk], len(present)))
print("unambiguous phases", len(rows))
kap0 = tuple([0] * 4)
cnt = collections.Counter()
by_kid = collections.defaultdict(collections.Counter)
for x, y, kd, miss, npres in rows:
    key = "none missing" if not miss else ("only (0,0,0,0) missing" if miss == [[0, 0, 0, 0]] else "missing %d cosets" % len(miss) + (" incl. 0" if [0, 0, 0, 0] in miss else ""))
    cnt[key] += 1; by_kid[kd][key] += 1
print("relative cosets with NO present covering tile at q, over phases:")
for k, v in cnt.most_common(10): print("  %-30s %d" % (k, v))
print("by child piece of Z_0(%s) (first-map index):" % S1)
for kd in sorted(by_kid): print("  child %2d" % kd, dict(by_kid[kd]))
json.dump(dict(q=[QX, QY], rows=rows, Qcosets=sorted(set(map(str, [x[4] for x in Q])))), open("figures/spectre_coset_realphase_q%.2f_%.2f_r%s_d%s_p%g.json" % (QX, QY, os.environ.get("RADQ","7"), os.environ.get("DCUT","def"), PBOX), "w"))
