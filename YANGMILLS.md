# Exit basins of homogeneous SU(2) Yang-Mills mechanics

Module `lakes_of_wada/yangmills.py`, tests in `tests/test_yangmills.py`, analysis `figures/analyse_yangmills.py` (about 8 minutes), figure `figures/yangmills.png`, numbers `figures/yangmills_analysis.csv`.
Labels: **Exact** (identity or test that can fail), **Counted** (a computed number with its conditions), **Not shown**.

## 1. The system and what an exit is

H = (p_x^2 + p_y^2)/2 + x^2 y^2 / 2, the spatially homogeneous SU(2) Yang-Mills field in a two-variable truncation (Matinyan, Savvidy and others; not re-derived here). The scaling x -> s x, t -> t/s, E -> s^4 E lets us set E = 1.

The potential vanishes on the four half axes, so there are four valleys, of width about 1/|x|. A particle in a valley is not free: the transverse oscillation has frequency |x|, its adiabatic invariant J = E_transverse/|x| acts as a potential J|x| for the motion along the valley, so a given orbit is
turned back (this is a heuristic argument, not shown here). The system has no escape to infinity in the sense of the Henon-Heiles scattering of this repository. `PHYSICS.md` 4.3 item 2 says the same (that homogeneous Yang-Mills is bound and has no exits). What can be defined is a **first passage**: the label of an orbit is the half
axis (1: +x, 2: +y, 3: -x, 4: -y) on which it first crosses the circle of radius R. For R above (8E)^(1/4) = 1.68 the circle lies outside the allowed region x^2 y^2 < 2E except in the valleys, so the label is unambiguous. **The label depends on R**, and so does everything below: this is a family of
problems, not one exit-basin problem.

Initial conditions: the section x = 0 (V = 0 there for every y), (y, p_y) in [-2, 2] x (-sqrt 2, sqrt 2), p_x = + sqrt(2E - p_y^2). Integration by velocity Verlet.

- **Exact (tests).** The map has the symmetry (y, p_y) -> (-y, -p_y) with labels 2 and 4 exchanged (0 disagreements out of 400); energy drift below 1e-3 (about 1e-4 in the runs below); invalid states (p_y^2 >= 2E) are marked.

## 2. Counted results (E = 1)

| R | undecided (t_max) | neighbour agreement of labels, n = 200 (shuffled-label control) | step floor, dt 0.001 against 0.0005 | uncertainty exponent alpha (D = 2 - alpha) |
|---|---|---|---|---|
| 2 | 0 (60) | 0.87 (0.27) | 0.000 | 0.43 (1.57) |
| 3 | 0 (100) | 0.62 (0.26) | 0.069 | 0.147 (1.85) |
| 4 | 0.001 (150) | 0.51 (0.25) | 0.177 | 0.088 (1.91) |

The uncertainty fraction f(eps) is for random perturbations of (y, p_y) of size eps in 7 values from 1e-4 to 0.1 (6000 states, p_x recomputed), with alpha the log-log slope; at R = 3 the local slopes are 0.12 to 0.16 over three decades, so this
one is a clean power law. At R = 4 the values of f are only just above the step floor (the floor is 0.18, the smallest f is 0.32) and the exponent is reliable only as an upper bound on how far D can be below 2. At R = 6 the floor is 0.37 and no exponent was fitted. **D grows with R**: 1.57, 1.85, 1.91. The longer an orbit is followed before it is
assigned a channel, the more fractal the boundary, which is what a first-passage label does when there is no true escape; the number is a property of the chosen R, and I do not read it as the dimension of anything physical.

The labels themselves are not converged in the step at R >= 3: between dt 0.002 and 0.001, 8 % of the pixels change label at R = 3 and 37 % at R = 6 (about a third of them at R = 4), falling only slowly as dt is halved (0.098, 0.084, 0.071 at R = 3 for 0.004 to 0.002, 0.002 to 0.001, 0.001 to 0.0005). The size of this change matches the uncertainty fraction at an equivalent
eps, so it is the same hypersensitivity, not an error that a better integrator removes; the labels at R = 2 are converged (floor 0).

## 3. The Wada test

Wada profile (fraction of boundary pixels that see all four basins within r pixels), n = 200:

| R | r = 1 | 2 | 3 | 4 | 8 |
|---|---|---|---|---|---|
| 2 | 0.078 | 0.50 | 0.85 | 0.94 | 0.99 |
| 3 | 0.149 | 0.71 | 0.96 | 0.99 | 1.0 |
| 4 | 0.168 | 0.78 | 0.98 | 1.0 | 1.0 |
| R = 3, labels shuffled over the pixels (control) | 0.217 | 0.87 | 1.0 | 1.0 | 1.0 |

Radius needed to see all four basins, in pixels, percentiles 50, 90, 99, at n = 200, 400, 800:
R = 2: 2.2, 3.6, 8.1 then 2.2, 3.6, 9.5 then 2.2, 4.1, 10.6; R = 3: 2.0, 2.8, 3.6 then 1.4, 2.2, 3.6 then 1.4, 2.2, 3.6.

- At **R = 3** the needed radius stays constant in pixels as the grid is refined from 200 to 800, so in physical units it shrinks in proportion to the resolution, which is the criterion this repository uses for a Wada boundary; the profile does not move when dt is halved three times (0.958 at r = 3 for dt from 0.004 to 0.0005).
- At **R = 2** the 99th percentile grows slowly in pixels (8.1, 9.5, 10.6), so in physical units it falls by a factor of about 0.6 instead of 0.5 per doubling: a small number of boundary points are far from at least one basin. Not Wada at all points at these resolutions, and not ruled out.
- **The test cannot separate Wada from noise at R >= 3.** The shuffled-label control has a profile above the real one at every radius. The labels do have spatial structure (neighbour agreement 0.62 against 0.26 shuffled, rising to 0.65 at n = 400), but at the pixel scale the map is so sensitive that the boundary pixels are
  surrounded by all four labels in any case. The honest reading: consistent with Wada, not demonstrated, and not a distinct property from the intermingling of the labels at the resolution of the grid.

## 4. What this says

- **Says.** An exit-basin problem for homogeneous Yang-Mills mechanics can be defined by first passage to radius R, it has a sharp symmetry, its basins are intermingled at every resolution tried, and the uncertainty exponent at R = 3 is 0.147 (D = 1.85), a number in the range of the Mixmaster exit basins of `MIXMASTER.md` that I did not expect and that I have not explained.
- **Does not say.** That the basins are Wada in the sense of the repository's diagnostics, that D is R independent (it is not), or anything about the quantum Yang-Mills theory, the mass gap, or `4D_Wilson_Loop_RG`. Item 2 of `PHYSICS.md` 4.3 proposed gauge-Higgs variants with several degenerate vacua; this is not that, it is the pure x^2 y^2 system with an operational exit.
- **Not done.** The same with a Higgs potential that has true final states; a section other than x = 0 (the full energy surface is three-dimensional); other values of E are equivalent by scaling, and so is the choice of window only through the section; a time-resolved (instead of radius-resolved) exit; the limit R -> infinity.

## 5. The gauge-Higgs toy with true final states (`higgs_basins`, `figures/analyse_higgs.py`, `figures/higgs_basins.png`)

To remove the arbitrary radius, the valleys are closed off by a Higgs-like term: V = x^2 y^2 / 2 + (lambda/4)(x^2 + y^2 - 1)^2. It has four degenerate minima, (+-1, 0) and (0, +-1), at the foot of the four valleys (**Exact**, a test), the symmetry of the square, and, with friction gamma
(q'' = -grad V - gamma q'), every orbit ends in one of them. **This is a toy chosen for this repository, not a derived truncation of Yang-Mills-Higgs theory**; its Higgs sector (a single radial term) and the friction are assumptions. Orbits are released from rest on a grid of positions in [-1.8, 1.8]^2 and
labelled by the well they settle in (1: +x, 2: +y, 3: -x, 4: -y); settled means within 0.03 of a minimum at speed below 0.03, t_max = max(25/gamma, 40), RK4 with dt 0.01.

- **Exact (tests).** The picture has the symmetry of the square: (x, y) -> (y, x) exchanges labels 1 and 2, 3 and 4, and x -> -x exchanges 1 and 3 (0 disagreements out of 300); at gamma = 1 the basins are smooth (neighbour agreement above 0.95).

**Counted, lambda = 1** (n = 200 grid, 4000 states for f(eps) with eps from 1e-4 to 3e-2, dt 0.01 against 0.005; except gamma = 0.02, run on n = 100 and 200 and 2000 states):

| friction | undecided | neighbour agreement (shuffled 0.25) | uncertainty exponent alpha, D = 2 - alpha | Wada profile at r = 2, 4 px | 99th percentile of the radius needed, px (n, then 2n) | step floor |
|---|---|---|---|---|---|---|
| 1.0 | 0.010 | 0.99 | alpha above 1: boundaries are smooth, D = 1 | 0.01, 0.04 | 57, 115 | 0 |
| 0.3 | 0.007 | 0.96 | 1.07 (D = 1) | 0.03, 0.15 | 33, 66 | 0 |
| 0.1 | 0.007 | 0.70 | 0.55 (1.45) | 0.39, 0.85 | 14.5, 25.8 | 0 |
| 0.05 | 0.007 | 0.51 | 0.21 (1.79) | 0.72, 0.93 | 10.3, 13.0 | 0.001 |
| 0.02 | 0.015 | 0.41 | 0.03 (1.97) | 0.80, 0.99 | 4.2, 5.1 | 0.225 |

(The shuffled-label control has profile 0.89 to 0.90 at r = 2 and 1.0 at r = 4 in every row.) At gamma >= 0.3 an alpha above 1 only means that the boundaries are curves; the estimate has no meaning there beyond 'D is 1'. At gamma = 0.02 the labels depend on the step (floor 0.23) and f(eps) is above 1.5 times the floor
only slightly (0.51 to 0.60), so D = 1.97 is a statement that the boundary is almost everywhere and not an exponent.

**Other lambda, gamma = 0.1:** lambda = 0.3 gives alpha 0.55 (D = 1.45, agreement 0.78, profile 0.71 at r = 4, 99th percentile 23 then 43 px); lambda = 3 gives alpha 0.27 (D = 1.73, agreement 0.55, profile 0.97 at r = 4, 99th percentile 7.1 then 8.5 px). D at lambda = 0.3 and 1 is the same to the digits shown and larger at lambda = 3.

**What the picture shows (post hoc).** The figure has a central region, |x| + |y| below about 1, in which each basin is a smooth diamond: orbits released below the barrier go to the nearest well, and the diagonals x = +-y are separatrices between two basins, which cannot be Wada.
Outside it, the high-energy region, is where the fractal structure grows as the friction falls. Splitting the boundary by that region (the threshold |x| + |y| = 1.0 and 1.5 was picked after seeing the figure):

| friction | outer boundary, profile at r = 2, 4 | inner boundary, profile at r = 2, 4 |
|---|---|---|
| 0.3 | 0.04, 0.17 | 0.01, 0.05 |
| 0.1 | 0.42, 0.90 | 0.04, 0.20 |
| 0.05 | 0.80, 1.00 | 0.18, 0.55 |
| 0.02 (n = 100) | 0.87, 1.00 | 0.68, 0.95 |

**Reading.** With weak friction the outer boundary sees all four wells within a few pixels and the radius needed grows only mildly with the resolution (99th percentile: 10.3 to 13.0 px at gamma = 0.05, 4.2 to 5.1 at 0.02, against a doubling for the smooth cases), which is the Wada signature of the Fernandez-Seoane-Sanjuan type of result
that weak dissipation helps (`PHYSICS.md`, section 1; I did not re-read that paper for this). It is not as sharp as at R = 3 of section 3, where the radius was constant in pixels; here the 99th percentile still creeps up. The same limits as section 3 apply: the shuffled-label control saturates before the real profile, so the test cannot
tell Wada from pixel noise at the resolutions used, and the labels at gamma = 0.02 are not converged in the step. The smooth inner region is not Wada at any friction in the table. The value at gamma = 0.1 (D = 1.45) is a partly fractal boundary, not an everywhere-fractal one.

**Not done.** More wells (the number of final states is fixed at four by the square symmetry of the model; the repository's n-channel results suggest nothing special about four); a derived Higgs sector or a time-dependent coupling; friction smaller than 0.02; larger grids at weak friction (the cost grows as 1/gamma); the same section in 3D phase space (released with momentum).
