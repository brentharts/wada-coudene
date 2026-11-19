"""Neighbouring windows of a Spectre tile in a covering of internal space, from the placements that cover a probe point.

For a tile T0 (window class a: shape S1 and frame A1) and a probe point P inside T0, the tile placements (S2, A2, v) (v a lattice vector) that
cover P have windows that partition internal space (a tiling covers each point of the plane once).  This script finds them by the method of
brentharts/spectre boundary_certificate/window_boundary_closure.py (affine fit of the window classes against a level-N patch, lattice offsets v whose
placed tile overlaps T0 physically and whose window lies within 2 RAD) with one change: it keeps only the placements that contain P.
Each is stored as (class b, v, S2, rotation B = R1^T R2, offset d) in the frame of the window of T0, i.e. the neighbouring window is B Z_0(S2) + d.
Usage: SPECTRE_DIR=/path/to/spectre python3 figures/spectre_window_neighbours.py N K [class ...]     -> figures/spectre_window_neighbours_N{N}_K{K}.npz
(N: levels of the patch used to fit the classes; K: box of lattice offsets.)"""
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
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
K = int(sys.argv[2]) if len(sys.argv) > 2 else 10
only = [int(x) for x in sys.argv[3:]]
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
def rec(s, j, T, z):
    if j == 0:
        g = groups.setdefault((tuple(T[0].flatten()), 'G' if s == 'Gamma' else 'O'), [T[0], [], []])
        g[1].append(STARm * mp.matrix([mp.mpf(int(x)) for x in T[1]])); g[2].append(z); return
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
rng = np.arange(-K, K + 1)
V = np.array(list(itertools.product(rng, rng, rng, rng)))
ph, it = V @ PH.T, V @ STAR.T
DPH = 12.0
cand = V[(np.linalg.norm(ph, axis=1) < DPH) & (np.linalg.norm(it, axis=1) < 2 * RAD * scale + 1e-9 + 2 * max(np.linalg.norm(f[3]) for f in FR))]
print('lattice candidates', len(cand), ' (box edge %d; boundary hits %d)' % (K, int((np.abs(cand) == K).any(1).sum())), flush=True)
z4 = np.zeros(4, dtype=object)
out = {}
for a in (only or range(len(FR))):
    S1, A1, R1, c1, _, _ = FR[a]
    P1 = shape(S1, A1, z4)
    Pt = P1.representative_point()
    rows = []
    for b in range(len(FR)):
        S2, A2, R2, c2, _, _ = FR[b]
        for v in cand:
            d = R1.T @ ((c2 - c1 - STAR @ v) / scale)
            if np.linalg.norm(d) > 2 * RAD + 1e-9:
                continue
            vo = np.array([int(x) for x in v], dtype=object)
            sh = shape(S2, A2, vo)
            if not sh.contains(Pt):
                continue
            rows.append((b, *[int(x) for x in v], 0 if S2 == 'O' else 1, *(R1.T @ R2).ravel(), *d))
    out[f"a{a}"] = np.array(rows, dtype=float)
    out[f"a{a}_probe"] = np.array([Pt.x, Pt.y]); out[f"a{a}_shape"] = np.array([0 if S1 == 'O' else 1])
    print(f"class {a} ({S1}): probe ({Pt.x:.3f}, {Pt.y:.3f}); {len(rows)} placements cover it (including the tile itself)  {time.time() - t00:.0f}s", flush=True)
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"spectre_window_neighbours_N{N}_K{K}.npz"), **out)
