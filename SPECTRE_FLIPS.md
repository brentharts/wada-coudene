# Spectre phason flips: engine, checks and results

Companion to `YANGMILLS.md` sections 11 and 12 (the four window pieces and the covering by neighbouring windows). Here the windows are used the other way round: a tiling is a point `z` of internal space, a phason step `z -> z + delta*e` moves that point, and a **flip** is a step after which the set of tiles present has changed. The paper is `quantgrav.tex` in `spectre-monotile-py` (section "Phason flips"); the run folder is `research/toymodel_v8` there.

Labels: **Exact** (proved or by construction), **Counted** (measured, one class, one footprint), **Hypothesis**, **Negative**.

## 1. Setup

- A tile `(b, u+w)` is present at phase `z` iff `B_b z + d_{b,w}` lies in the window `Z_0(S_b)`. The windows are graph-directed fractals with ratio `1/lambda`, `lambda^2 = 4 + sqrt(15)`. **Exact.**
- All tile vectors lie in one coset of an index-36 sublattice of `Z^4`; the coset offset is the constant 0, so no finer coset label is needed. **Exact** (by the construction in `spectre_coset_label.py`).
- There are 12 window classes: 6 single tiles (`O`) and 6 Mystic pairs (`G`). Area ratio `|Z_0(G)| / |Z_0(O)| = 1/(3 + sqrt(15)) = 0.1455`; measured 0.1460 at tolerance 0.0015. **Exact / Counted.**
- Flip: ordered pair `(z, z' = z + delta*e)` with different present tile sets. `F(z, z') = T(z) triangle T(z')` is symmetric under reversal, so flips alone give no arrow of time. **Exact.**

## 2. Engine parameters that matter

| Parameter | Value | Why |
| :--- | :--- | :--- |
| `PBOX` (lattice box, `|PH w| < PBOX`) | 20 (default in `spectre_flip_*.py`) | A box of 12 dropped Mystic-pair tiles; 5 of 1934 real neighbour pairs were lost and 12% of `Z_0(O)` phases showed holes at the probe. Box 20 finds all 1934. |
| `IN_T / OUT_T` (membership tolerance) | 0.0015 / 0.0035 | A flip happens only where `z` is within `delta` of a boundary, so `IN_T` close to `delta` makes the old and the new tile both "present". With `IN_T = 0.003`, `delta = 0.01`: overlaps in 37 of 60 changed pairs, 0 of 60 unchanged. With `0.0015`, `delta = 0.03`: 5 and 1 of 60; 1180 of 1321 (89%) of pairs valid. |
| `delta` | 0.03 | Output file names carry it (`*_d0.03*`). |
| Cloud spacing at depth 7 | 0.0012 | Resolution of the window boundary. |

The covering script `spectre_window_covering.py` gives the same class-0 placement list at box 12 and box 20 (368 overlapping placements, 92 containing the probe), so section 12 of `YANGMILLS.md` is unaffected for class 0. Classes 2 and 5 have not been rerun at box 20.

## 3. Results (class 0, level 3, delta 0.03, box 20, tight tolerance; **Counted**)

- About 10 tiles change per phase within radius 14, uniform in distance. Changed tiles in a disc grow as `rho^2.25`; changed per present tile 0.056 (a lower bound).
- 818 clusters of changed tiles. The first twelve inspected have 8 to 17 tiles per side and area 80 to 172. None closed inside the validated region (radius 8.5; closure needs about 26), so no `(dV, dE, dF)`.
- Contact graph of the tiles meeting a fixed footprint: of 539 changed pairs, 468 valid; 224 (48%) change the tile set; 169 of those (75%) give a non-isomorphic graph, 55 an isomorphic graph with different classes; 157 of the 224 change the tile count. The reversed step gives the same counts.
- The tetrad flux of the lift is zero for every flip (the lift forms are exact), while the contact graph changes. This is **not** a Comisso-Asenjo reconnection: there is no field, no non-ideal term, no helicity on a planar tile. **Negative** for that reading.
- Parity: orientation-parity minority tiles are the Mystic tiles, frequency `(5 - sqrt(15))/10 = 0.1127`, ratio of the two frequencies `4 + sqrt(15) = lambda^2`. Tile counts at levels 0 to 5: 1, 9, 71, 559, 4401, 34649; minority counts at levels 1 to 5: 1, 8, 63, 496, 3905.
- In a closed cluster `n_O + 2 n_G` (area) is conserved; `n_G` changes in 0.32 of 163 area-conserving events at `rho = 7`, not smoothly in `rho`. **Counted.**
- Erasing one tile never leaves the inside ambiguous given the outside (P = 0 for multiplicity 1, 2, 3 or more; n = 143, 22, 19). **Negative** for a Landauer-type cost from this test.

## 4. What is not done

Closed clusters (needs radius about 26, long run); step scan (0.015, 0.0075) with the tight tolerance; the other 11 classes and level 4; flip rate against level (the inflation arrow); Pachner-move classification; a map to cohomology of the hull.

## 5. Reproduction

Scripts in `figures/` (they need `SPECTRE_DIR`, a clone of `brentharts/spectre`): `spectre_coset_label.py`, `spectre_coset_realphase.py` (`PBOX` environment variable), `spectre_realpairs.py`, `spectre_flip_contacts.py`, `spectre_flip_stats.py`, `spectre_flip_moves.py`, `spectre_window_covering_pbox.py`; plots `plot_coset_label.py` and, in `spectre-monotile-py`, `research/toymodel_v8/plot_rerun.py`. The stats and contact runs take well over ten minutes each; run them singly, two at once ran out of memory.
