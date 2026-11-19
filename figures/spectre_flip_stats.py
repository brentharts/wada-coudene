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
# ======================= phason flips with the contact graph =======================
# Placements (b, w) whose footprint lies within R2 of the footprint of the tile a at the origin; a (b = a, w = 0) is included.
# At internal point z (frame of Z_0(S_a)) the tile (b, w) is present iff z_b = B z + d lies in Z_0(S_b).
R2 = float(os.environ.get("R2", "1.0")); DIAM = float(os.environ.get("DIAM", "0")); DELTA = float(os.environ.get("DELTA", "0.03")); NS = int(os.environ.get("NS", "3000"))
IN_T, OUT_T = float(os.environ.get("IN_T", "0.008")), float(os.environ.get("OUT_T", "0.016"))                           # nearest-point distance: in if < IN_T, out if > OUT_T, else ambiguous
PL = []
for b in range(len(FR)):
    Sb, Ab, Rb, cb, _, _ = FR[b]
    v_b = np.array(VEC[keys[b]][0], dtype=np.int64)
    W = cands[in_lattice(cands - (v_b - u_a))]
    dd = ((c1 - cb)[None, :] + W @ STAR.T) @ Rb / scale
    keep = np.linalg.norm(dd, axis=1) <= 2 * RAD + 0.6      # z may lie up to 0.3 beyond the window radius
    for w, d in zip(W[keep], dd[keep]):
        poly = shape(Sb, Ab, np.array([int(x) for x in w], dtype=object))
        if poly.distance(P1) < R2:
            PL.append((b, Sb, Rb.T @ R1, d, w, poly))
print("tile bounds", P1.bounds, "max class diam", max(shape(f[0],f[1],np.zeros(4,dtype=object)).minimum_rotated_rectangle.length for f in FR));print("placements within %.2f of the tile: %d" % (R2, len(PL)), flush=True)
own_idx = [k for k, p in enumerate(PL) if p[0] == a and (p[4] == 0).all()][0]

import json, collections
from shapely.ops import unary_union as UU
rng_ = np.random.default_rng(11)
polys = [p[5] for p in PL]
isG = np.array([p[1] == "G" for p in PL])
cent = np.array([[q.centroid.x, q.centroid.y] for q in polys])
c0 = np.array([P1.centroid.x, P1.centroid.y])
def state(Zs):
    pres = np.zeros((len(Zs), len(PL)), bool); amb = np.zeros((len(Zs), len(PL)), bool)
    for k, (b, Sb, B, d, w, poly) in enumerate(PL):
        dist = trees[Sb].query(Zs @ B.T + d)[0]
        pres[:, k] = dist < IN_T; amb[:, k] = (dist >= IN_T) & (dist <= OUT_T)
    return pres, amb
rad = np.linalg.norm(cent - c0, axis=1)
RHOS = [1.5, 2, 3, 4, 5, 6, 7]
# ---------------- Part A: changed tiles in a disc (H8.5) and the Mystic-pair count (H8.4) ----------------
Zc = trees[S1].data; lo, hi = Zc.min(0) - 0.3, Zc.max(0) + 0.3
NSA = int(os.environ.get("NSA", "3500"))
Z0 = rng_.uniform(lo, hi, (NSA, 2)); ph = rng_.uniform(0, 2 * np.pi, NSA); Z1 = Z0 + DELTA * np.stack([np.cos(ph), np.sin(ph)], 1)
Pa, Aa = state(Z0); Pb, Ab = state(Z1)
ok = ~(Aa.any(1) | Ab.any(1)); print("part A: samples %d usable %d" % (NSA, ok.sum()), flush=True)
A = {}
for rho in RHOS:
    m = rad <= rho
    rem = Pa[ok][:, m] & ~Pb[ok][:, m]; add = ~Pa[ok][:, m] & Pb[ok][:, m]
    nch = rem.sum(1) + add.sum(1)
    gm = isG[m]
    dG = add[:, gm].sum(1) - rem[:, gm].sum(1); gch = add[:, gm].sum(1) + rem[:, gm].sum(1)
    dN = add.sum(1) - rem.sum(1)
    dA = dN + dG                                       # change of n_O + 2 n_G = area / 8.196
    ntiles = Pa[ok][:, m].sum(1)
    anych = nch > 0
    ac = anych & (dA == 0)
    A[rho] = dict(mean_present=float(ntiles.mean()), mean_changed=float(nch.mean()), frac_any=float(anych.mean()),
                  mean_G_changed=float(gch.mean()), P_dG_nonzero_given_change=float((dG[anych] != 0).mean()) if anych.any() else None,
                  P_dN_nonzero_given_change=float((dN[anych] != 0).mean()) if anych.any() else None, mean_abs_dG=float(np.abs(dG).mean()), mean_dG=float(dG.mean()), P_dArea_nonzero_given_change=float((dA[anych] != 0).mean()) if anych.any() else None, n_area_conserving=int(ac.sum()), P_dG_nonzero_given_area_conserved=float((dG[ac] != 0).mean()) if ac.any() else None, mean_dN_given_area_conserved=float(dN[ac].mean()) if ac.any() else None)
    print("rho %.1f  present %.1f changed %.3f any %.3f  G-changed %.3f  P(dG!=0|change) %.3f  P(dN!=0|change) %.3f  P(dArea!=0|change) %.3f  area-conserving events %d: P(dG!=0) %.3f mean dN %.2f" % (rho, A[rho]["mean_present"], A[rho]["mean_changed"], A[rho]["frac_any"], A[rho]["mean_G_changed"], A[rho]["P_dG_nonzero_given_change"] or -1, A[rho]["P_dN_nonzero_given_change"] or -1, A[rho]["P_dArea_nonzero_given_change"] or -1, A[rho]["n_area_conserving"], A[rho]["P_dG_nonzero_given_area_conserved"] if A[rho]["P_dG_nonzero_given_area_conserved"] is not None else -1, A[rho]["mean_dN_given_area_conserved"] if A[rho]["mean_dN_given_area_conserved"] is not None else 0), flush=True)
lr = np.log(RHOS); lm = np.log([max(A[r]["mean_changed"], 1e-9) for r in RHOS])
sl = np.polyfit(lr[1:], lm[1:], 1)[0]; print("log-log slope of changed tiles against rho (rho 2 to 7): %.3f (rho^2 -> 2)" % sl)
# ---------------- Part B: erasure of the tile at the origin, ambiguity against the multiplicity of outcomes (H8.6) ----------------
NB = int(os.environ.get("NB", "200")); M = int(os.environ.get("MB", "50")); EPS = float(os.environ.get("EPS", "0.05"))
inside = np.array([polys[k].intersection(P1).area > 1e-9 for k in range(len(PL))])
outside = ~inside
rows = []
zs = rng_.uniform(lo, hi, (NB, 2))
for z in zs:
    u = rng_.uniform(-1, 1, (4 * M, 2)); u = u[np.linalg.norm(u, axis=1) <= 1][:M]
    Zb = np.vstack([z[None, :], z[None, :] + EPS * u])
    P, Am = state(Zb)
    good = ~Am.any(1)
    if not good[0]: continue
    P, Zb = P[good], Zb[good]
    out0 = P[0][outside]
    cons = np.array([np.array_equal(p[outside], out0) for p in P])
    pats_all = set(tuple(np.where(p[inside])[0]) for p in P)
    pats_cons = set(tuple(np.where(p[inside])[0]) for p, c in zip(P, cons) if c)
    half = np.linalg.norm(Zb - z, axis=1) <= EPS / 2
    pats_all_h = set(tuple(np.where(p[inside])[0]) for p, h in zip(P, half) if h)
    pats_cons_h = set(tuple(np.where(p[inside])[0]) for p, c, h in zip(P, cons, half) if c and h)
    rows.append((len(pats_all), len(pats_cons), len(pats_all_h), len(pats_cons_h), int(cons.sum()), len(P)))
rows = np.array(rows) if rows else np.zeros((0, 6), int)
print("part B: z used %d" % len(rows))
for tag, (ja, ic) in {"eps": (0, 1), "eps/2": (2, 3)}.items():
    J, I = rows[:, ja], rows[:, ic]
    for lab, sel in (("multiplicity >= 3", J >= 3), ("multiplicity 2", J == 2), ("multiplicity 1", J == 1)):
        if sel.any(): print("  %-6s %-18s n=%3d  P(inside not determined by outside) = %.3f   mean consistent fraction %.3f" % (tag, lab, sel.sum(), (I[sel] >= 2).mean(), (rows[sel, 4] / rows[sel, 5]).mean()))
json.dump(dict(DELTA=DELTA, R2=R2, IN_T=IN_T, OUT_T=OUT_T, NSA=NSA, usable=int(ok.sum()), slope=float(sl), partA={str(k): v for k, v in A.items()},
               NB=NB, M=M, EPS=EPS, partB=rows.tolist()), open("figures/spectre_flip_stats_d%s.json" % DELTA, "w"))
