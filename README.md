# wada-coudene

Python/Sagemath scripts to make pictures of "Lakes of Wada" sets in the plane and in the 2-torus, using the method described by Yves Coudene in [Pictures of Dynamical Systems, Notices of AMS 53, No 1, pp8-13](http://www.ams.org/notices/200601/fea-coudene.pdf).

## Author

* [David Dumas ](http://dumas.io/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Rewrite (numpy/scipy, no Sage)

`wada-plane.py` and `wada-torus.py` are kept unchanged as reference. The
maintained code is the `lakes_of_wada/` package plus the `wada.py` front end:

    python3 wada.py torus -n 400 -o torus.png
    python3 wada.py plane -n 400 -o plane.png --data plane.npz   # also save raw arrays

Package layout: `maps` (Anosov map + pluggable perturbations), `classify`
(orbit/basin classification, colouring rules), `views` (pixel -> torus point),
`elliptic` (inverse Weierstrass P via scipy), `pipeline` (glue). Run tests with
`python3 -m pytest tests`.

### Analysis tools and plane dynamics

    python3 wada.py entire --family newton --degree 3 -n 600 --wada   # basins of z^3-1 + Wada profile
    python3 wada.py entire --family exp --lam 0.3 --radius 5 -o exp.png
    python3 wada.py torus -n 400 --maverick mav.png                   # orbit-divergence map
    python3 wada.py canal -n 256 --wada                               # canal construction
    python3 wada.py scatter -n 500 --energy 0.25 --wada --merging     # Henon-Heiles exit basins
    python3 wada.py scatter -n 64 --dim 3 --data hh3d.npz --wada      # 3D energy surface

* `lakes_of_wada.wadatest`: numerical Wada tests on a label array of **any
  dimension** (2D images, 3D volumes, ...). Labels: `>0` basins, `0` undecided
  (counted as boundary), `<0` invalid/outside the domain (ignored).
  * *Wada profile* (`--wada`): fraction of boundary points that see all K basins
    within r pixels; `basin_profile` does it per basin, to expose partially Wada
    sets. Verified controls: Newton for z^3-1 reaches ~1.0 by r = 4 px, an
    ordinary 3-sector partition stays near 0.
  * *Merging method* (`--merging`, after the Wagemakers-Daza-Sanjuan review):
    merge all basins but one, compare the resulting slim boundaries by Hausdorff
    distance. Judge by `max_d_normalized` (small = Wada) as well as by the
    ratio: a plain 3-sector partition has nearly identical boundaries
    (ratio ~ 0) but they sit ~40% of the diagonal apart.
* `lakes_of_wada.scattering`: continuous-time open systems. A vectorized RK4
  integrator with an exit classifier gives exit basins; Henon-Heiles (`E > 1/6`,
  three exits) is included, both as the 2D section x = 0 in (y, p_y) and as the
  full 3D energy surface. Validated against scipy's DOP853 and for energy
  conservation. At E = 0.25 the exit basins come out Wada by both tests, in 2D
  (profile 0.995 at r = 4 px, merging distance 0.8% of the diagonal) and in 3D
  at 48^3 to 64^3 (profile 0.98 at r = 3 voxels).
* `lakes_of_wada.maverick`: per-pixel *orbit divergence*: the largest distance
  over a late window of iterates between a boundary pixel's orbit and the orbit
  of the nearest settled interior pixel of each basin (a finite-window proxy for
  the limsup in the maverick-point definition below). Treat it as a heat map:
  the maps are chaotic, so double-precision orbits are only pseudo-orbits.
* `lakes_of_wada.entire`: maps of the complex plane (Newton, polynomial, `exp`,
  `sin`; any vectorized callable works), classified into escape / attracting
  cycle / undecided, with the chordal metric for the analysis tools.
* `lakes_of_wada.construction`: Wada-style canal construction on a raster (sea
  plus ponds; stage k digs lake k mod K with shrinking tolerance and canal
  width). Resolution-limited: at n = 512 the Wada profile is 0.61 at r = 16 px,
  0.90 at r = 32, 1.00 at r = 64; larger n goes deeper.

## Related work

* D. Martí-Pete, L. Rempe, J. Waterman, *Eremenko's conjecture, wandering Lakes
  of Wada, and maverick points*, J. Amer. Math. Soc. (arXiv:
  [2108.10256](https://arxiv.org/abs/2108.10256)). Constructs transcendental
  entire functions with infinitely many Fatou components sharing a common
  boundary (Lakes of Wada continua), and introduces *maverick points*: boundary
  points of a wandering domain whose orbits accumulate somewhere the domain's
  own iterates do not. Their Figure 1 is the Plykin attractor, the perturbed
  Anosov torus map drawn here.
* A. Wagemakers, A. Daza, M. A. F. Sanjuán, *How to detect Wada basins*
  (arXiv: [2005.05923](https://arxiv.org/abs/2005.05923)). Review of four
  detection methods, each built on a different characterization of the Wada
  property: Nusse-Yorke (unstable manifolds of periodic orbits), the *grid*
  method (refine until a third colour appears between two others), the *merging*
  method (merge all other basins into one colour and compare the resulting
  boundaries by Hausdorff distance) and *saddle-straddle* (compare the chaotic
  saddles found for each basin versus the rest). Newton's method with 3 to 6
  roots is among its test systems.
* J. Aguirre, J. C. Vallejo, M. A. F. Sanjuán, *Wada basins and chaotic
  invariant sets in the Hénon-Heiles system*, Phys. Rev. E 64, 066208 (2001),
  [doi:10.1103/PhysRevE.64.066208](https://doi.org/10.1103/PhysRevE.64.066208).
  Shows by several computational methods that the exit basins of this open
  Hamiltonian are Wada, not just fractal. (Cited from its listing; the journal
  page was not readable when this entry was written.)
* D. S. Fernández, J. M. Seoane, M. A. F. Sanjuán, *Weak dissipation drives and
  enhances Wada basins in three-dimensional chaotic scattering*
  (arXiv: [2202.03069](https://arxiv.org/abs/2202.03069)). In a 3D open
  Hamiltonian system the basins stop being Wada above a critical energy, and
  weak dissipation turns non-Wada basins into partially Wada ones.
* Á. Daza, J. O. Shipley, S. R. Dolan, M. A. F. Sanjuán, *Wada structures in a
  binary black hole system*, Phys. Rev. D 98, 084050 (2018),
  [doi:10.1103/PhysRevD.98.084050](https://doi.org/10.1103/PhysRevD.98.084050)
  (arXiv: [1807.10741](https://arxiv.org/abs/1807.10741)). Treats the shadow of
  a binary black hole as exit basins of an open Hamiltonian system and uses a
  numerical Wada-detection algorithm to show parts of the shadow have the Wada
  property, separating fractal regions from regular ones.
* P. Pankka, J.-M. Wu, *On Lakes of Wada* (arXiv:
  [2601.02266](https://arxiv.org/abs/2601.02266)). Shows Lakes of Wada exist in
  the sphere S^n for every n >= 3, with lakes that are quasiconformally
  equivalent to a ball and are John domains.
* Y. Coudène, *Pictures of Dynamical Systems*, Notices AMS 53 (2006), the
  method used by the original scripts.

