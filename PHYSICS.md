# Physics connections


---

## 1. Where Wada basins already appear in physics (established)

The common setting is an *open* Hamiltonian system: a particle (or a light ray)
comes in, interacts, and leaves through one of several exits. The set of
initial conditions that leave through exit *k* is the exit basin of *k*. When the
motion is chaotic, the boundaries between exit basins are fractal, and for some
systems they are Wada: every boundary point of one basin is also a boundary point
of all the others.

* **Binary black holes.** Daza, Shipley, Dolan and Sanjuán treat the shadow of a
  binary black hole as the exit basins of an open Hamiltonian system. Using a
  numerical Wada-detection algorithm on the Majumdar-Papapetrou binary, they show
  that parts of the shadow have the Wada property, and that the algorithm
  separates the fractal parts of the shadow from the regular ones
  ([arXiv:1807.10741](https://arxiv.org/abs/1807.10741), Phys. Rev. D 98, 084050).
  This is the physical reading of the pictures in this repository: a "lake" is
  where a photon ends up (which black hole captures it, or where it escapes), and
  Wada-ness means that no boundary point of the shadow is next to only two of the
  outcomes.
* **Multi-black-hole scattering.** For three or more static extremal black holes,
  the scattering of light in the plane is chaotic, with a lower bound on the
  topological entropy of order log(N-1)
  ([Kluitenberg, Roest, Seri, arXiv:2211.06880](https://arxiv.org/abs/2211.06880)).
  Pantig studies the symbolic dynamics of multi-black-hole photon scattering
  (Physica D, December 2026,
  [doi:10.1016/j.physd.2026.135415](https://doi.org/10.1016/j.physd.2026.135415));
  I have only the title.
* **pp-wave spacetimes.** Geodesics of pp-waves with polynomial profiles are
  equivalent to a particle in a 2D harmonic polynomial potential, and the paper
  reports that the Wada property of the escape basins is robust as the degree
  (number of escape channels) increases (Rossetto et al., Chaos 36, 063136,
  [doi:10.1063/5.0336333](https://doi.org/10.1063/5.0336333); abstract only).
* **Hamiltonian chaotic scattering generally.** Henon-Heiles exit basins are Wada
  (Aguirre, Vallejo, Sanjuan, Phys. Rev. E 64, 066208); in a 3D system Wada-ness
  can disappear above a critical energy, and weak dissipation can enhance it
  (Fernandez, Seoane, Sanjuan, [arXiv:2202.03069](https://arxiv.org/abs/2202.03069)).

### What this repository adds (our own numerical results)

The `scattering` module integrates the Henon-Heiles system and its generalisation
to *n* exit channels, V = r^2/2 + Im(z^n)/n. The `wadatest` module measures the
Wada property on the resulting images. The results are in the README; briefly, at
3 x the energy at which the channels open, three to seven exit channels all
behave like Wada boundaries down to the pixel scale (the radius needed to see all
exits shrinks in proportion to 1/resolution), and at eight channels part of the
boundary stops doing so. This is a statement about this family of 2D potentials.

### Possible next experiment (not done)

The Majumdar-Papapetrou binary is a natural next system for the scattering
module: it is the setting of Daza et al., and a model of photon orbits there
would let the same diagnostics (Wada profile, needed-radius scaling) be applied
to a physical shadow. I have not derived the equations of motion for this
repository, and would do so from the original paper before using them.

---

## 2. Quarks, colour and the flux-tube picture

The idea under discussion: the three colour charges, and the fact that quarks
are never found alone, might be related to a three-way Wada-like junction of
spaces, with quarks as the places where three dimensions of space connect. Here
are the parts that touch established physics, and where the idea departs from it.

### 2.1 Y-junctions (established, with an analogy attached)

* Three static quarks must be joined by three colour-flux strings to form a
  gauge-invariant object. Two proposals for the baryon potential are the
  **Delta law** (a sum of two-body terms, as for the perimeter of the triangle
  of quarks) and the **Y law** (the length of three flux tubes joined at a
  junction, the Steiner point)
  ([de Forcrand and Jahn review, arXiv:hep-ph/0502039](https://arxiv.org/abs/hep-ph/0502039)).
* Lattice QCD results have shifted over time. Early quenched studies reported
  evidence for a Delta-type long-range potential and others for a genuine
  three-body force
  ([Baryonic flux after Abelian projection, arXiv:hep-lat/0401026](https://arxiv.org/abs/hep-lat/0401026)).
  A later high-accuracy calculation finds agreement with the Delta law at short
  distances and a departure from it for quark separations of roughly 0.8 fm and
  more (de Forcrand and Jahn, same review), in the direction of the Y law. Another
  group found the ground-state potential well described by Coulomb plus a
  Y-type linear term to within about 1%, and observed Y-shaped flux-tube
  formation in maximally-Abelian-projected lattice QCD
  ([Suganuma et al., arXiv:hep-lat/0412026](https://arxiv.org/abs/hep-lat/0412026)).
  So the Y-shaped junction is a real feature of the baryon in lattice QCD, at
  least at larger separations.
* **Where the idea reverses the standard picture.** In the standard picture the
  quarks are the endpoints of the strings and the junction sits inside ordinary 3D
  space. The proposal here puts the connection at the quark itself. The lattice
  results do not distinguish these, since they do not probe a model of space.
* **Analogy with this toolkit.** A baryon of SU(N) has N colours and N strings
  meeting at a junction. The harmonic-polynomial family in the `scattering` module
  has *n* exit channels around one central region. Both are "n one-dimensional
  channels meeting at a point". There is no dynamical derivation connecting them:
  flux tubes are not trajectories of a particle in a potential. If the analogy is
  worth anything, it is as a source of questions, for example whether the number of
  channels affects how intermingled the basins are. Our own result there is that
  three is not special for that (README, figure 10).

### 2.2 String breaking (established)

* Pulling heavy quarks apart stretches the flux tube. With light quarks present, the
  string is expected to break and a quark-antiquark pair to appear at the break;
  this is why free quarks are not seen. Observing it in lattice QCD was
  difficult: the plain Wilson-loop measurement does not show it, and one has to
  include the coupling to the meson-antimeson state
  ([Bernard, DeGrand, Sugar, arXiv:hep-lat/0103012](https://arxiv.org/abs/hep-lat/0103012)).
  The SESAM collaboration reported resolving the signature of string breaking in
  two-flavour QCD at zero temperature
  ([arXiv:hep-lat/0512018](https://arxiv.org/abs/hep-lat/0512018) and
  [hep-lat/0505012](https://arxiv.org/abs/hep-lat/0505012)).
* **Reading in the proposed picture (speculation).** If quarks are places where
  three dimensions connect, "pulling them apart" would try to disconnect the
  dimensions, and the creation of a pair would be the system re-establishing the
  connection. This is a restatement of the observation in different words, not
  yet a mechanism that predicts anything the standard picture does not.

### 2.3 Six 1D spaces (speculation)

The proposal is that 3D space is composed of separate half-lines, +X, -X, +Y, -Y,
+Z, -Z, stitched together at points, and that time is a single +T.

* **Constraint: rotational and Lorentz symmetry.** Experiments constrain departures
  from isotropy and Lorentz invariance very tightly. A model of stitched axes would
  need to recover full rotational symmetry in the continuum.
* **A relevant existing example.** Lattice gauge theory is built out of links along
  the directions of a hypercubic grid, with quarks on sites and gluons on links.
  The grid breaks rotational symmetry at finite spacing, and the symmetry is
  recovered as the spacing goes to zero. This shows that a discrete axis-aligned
  structure is compatible with a rotation-invariant continuum, but it is a
  property of the limit, not of any one lattice. It is also worth noting that a
  link in the -X direction from a site is, in the standard formalism, the
  Hermitian conjugate of the +X link from the neighbouring site, so +X and -X are
  not independent in that construction; the proposal would have to treat them as
  independent, which is a real difference.
* 

---

## 3. Negative temperature (established, and a caution)

* Negative absolute temperatures are states in which high-energy levels are more
  populated than low-energy ones. They are possible only in systems whose energy
  has an upper bound, so they are hotter than any positive temperature, not colder
  than absolute zero. They were first shown in localised spin systems, and a
  negative-temperature state for *motional* degrees of freedom was prepared with
  ultracold potassium atoms in an optical lattice, using a Hamiltonian made to
  have an upper bound on kinetic energy
  ([Braun et al., Science 339, 52 (2013), arXiv:1211.0545](https://arxiv.org/abs/1211.0545)).
* **Speculation under discussion:** that an ordered, locked-together state at very
  high energy might correspond to negative temperature, and that quark
  confinement might be the result of such fusing.
* **Evidence against, as I know it.** (i) Negative temperature needs an energy
  upper bound, and the energy of relativistic particles in quantum field theory is
  not known to have one. (ii) In QCD at positive temperature, lattice calculations
  find that matter heated past a crossover temperature of roughly 150 to 160 MeV
  becomes a quark-gluon plasma, in which quarks are deconfined; the HotQCD
  collaboration quotes T_pc = 156.5 +/- 1.5 MeV for the chiral crossover
  ([arXiv:2312.12857](https://arxiv.org/abs/2312.12857) quoting
  Phys. Lett. B 795 (2019) 15), and the transition is a smooth crossover
  ([Aarts, Physics 7, 86](https://physics.aps.org/articles/v7/86)). Hotter means less
  confined in this established picture, the opposite of what the proposal needs.
  (iii) The Wada property describes the boundaries of basins of a dynamical
  system, which is a statement about how initial conditions are sorted, not about a
  thermodynamic ordering.

---

## 4. The Yang-Mills mass gap (speculative section)

### 4.1 The problem (established)

The Clay Millennium problem (Jaffe and Witten) asks for a proof that, for any
compact simple gauge group, a non-trivial quantum Yang-Mills theory exists on
four-dimensional space (Euclidean R^4 in the standard statement) with a mass gap
Delta > 0, to the level of rigour of constructive quantum field theory
([problem statement summary](https://en.wikipedia.org/wiki/Yang%E2%80%93Mills_existence_and_mass_gap);
original: Jaffe and Witten, in
[The Millennium Prize Problems, Clay Mathematics Institute](https://www.claymath.org/library/monographs/MPPc.pdf)).
The existence part means axiomatic properties at least as strong as the Wightman
axioms. As of the pages I could read it is listed as unsolved, and the page notes
that a mathematically complete four-dimensional quantum gauge theory example is
not yet known. Unvetted claims of a solution circulate; this file treats the
problem as open.

### 4.2 What is, and is not, known (established)

* At the level of physics rather than proof, confinement and a gap are well
  supported: lattice computations find massive glueballs and a mass gap in pure
  Yang-Mills (the same page cites Teper, Lucini et al., and Morningstar and Chen
  et al.). That is evidence, not a proof of the continuum statement the prize
  asks for.
* Classical Yang-Mills theory with spatially homogeneous fields reduces to a
  Hamiltonian system with a few degrees of freedom, and it is chaotic. Matinyan,
  Savvidy and Ter-Arutyunyan-Savvidy (Sov. Phys. JETP 53, 421, 1981) showed a
  positive Lyapunov exponent; Savvidy later argued it is a Kolmogorov K-system. A
  recent paper summarises that this classical motion is chaotic at any energy
  density in the absence of other fields, and can be regularised only if other
  fields, for example a Higgs field, are included
  ([arXiv:0710.5691](https://arxiv.org/abs/0710.5691); the effect of matter fields
  on the instability is studied in [arXiv:0903.3990](https://arxiv.org/abs/0903.3990)).
  (The original papers are cited there; I have not read them.)

### 4.3 Speculative connections to explore

These are ideas to examine, ordered from the most concrete. None has support yet.

1. **Same family of mathematical objects, different statement.** The mass gap is a
   statement about the *spectrum* of a closed quantum system: there is a positive
   energy between the vacuum and the next state. Wada basin boundaries are a
   statement about *classical trajectories* in an *open* system. I know of no
   result that turns one into the other. Any link has to be built, not assumed,
   and a classical toy such as the ones in this repository cannot by itself say
   anything about the quantum gap.
2. **Chaotic classical Yang-Mills as a source of exit-basin problems.** Homogeneous
   Yang-Mills mechanics is chaotic but, as usually written, bound (the potential
   is confining), so it has no exits. Adding a Higgs potential with several
   degenerate vacua, or adding weak dissipation, would give a system with several
   *final states* whose basins could be mapped and tested for Wada-ness with this
   toolkit. That would be a concrete computation, and it would show only whether
   classical gauge-Higgs dynamics produces Wada basins. It would not address the
   gap.
3. **Confinement as a geometric statement (the proposal).** If colour-singlet
   composites are the only finite-energy states because quarks are junctions of
   three one-dimensional spaces, a gap could appear as the energy cost of
   disconnecting a junction. To be more than words this needs: a definition of the
   junction in a gauge-theory setting, a rule for its energy, and a reason that the
   cost is bounded below by a positive number independent of the size of the
   system. The lattice picture gives a sense of the last point (the flux tube
   energy grows with its length, and the string breaks), but only as a
   non-rigorous argument.
4. **The role of three.** Three is the number of colours in the physical theory,
   and the Millennium problem is posed for *any* compact simple group. A
   mechanism that depends on the number three being special would be a mechanism
   for SU(3) only, not for the problem as stated. This toolkit's n-channel results
   point the same way: nothing it measures singles out three.
5. **Centre vortices (unchecked).** I recall confinement pictures in terms of
   centre-vortex lines related to the centre of the gauge group, for SU(3) the
   group of three elements. I have not looked up the literature here, and list it as
   something to read before claiming any connection to the "three".

### 4.4 What would count as progress

* A precise mathematical statement of the proposal (what the spaces are, what the
  junction is, what the field is), in a form where its predictions can be compared
  with something established, such as the Y-shaped flux tube or string breaking.
* A concrete calculation, in a simple model, where the proposal predicts a number
  different from the standard picture.
* For the Millennium problem itself: existence of a theory on R^4 satisfying the
  stated axioms, together with the gap. Heuristic or numerical evidence, including
  anything computed in this repository, does not meet that standard.

---

## 5. Checklist of what could be done in this repository

* Majumdar-Papapetrou binary shadow in the `scattering` module (section 1).
* A weakly dissipative Henon-Heiles variant with a trapped fourth basin, as the
  closer analogue of Fernandez et al.
* Gauge-Higgs homogeneous Yang-Mills with degenerate vacua, mapped for Wada basins
  (section 4.3, item 2).
* Anything from the owner's own work on the mass gap, once available.

---

## References

* A. Daza, J. O. Shipley, S. R. Dolan, M. A. F. Sanjuan, *Wada structures in a
  binary black hole system*, Phys. Rev. D 98, 084050 (2018),
  [arXiv:1807.10741](https://arxiv.org/abs/1807.10741).
* M. Kluitenberg, D. Roest, M. Seri, *Chaotic light scattering around extremal
  black holes*, [arXiv:2211.06880](https://arxiv.org/abs/2211.06880).
* R. C. Pantig, *Symbolic dynamics of multi-black-hole photon scattering*,
  Physica D (2026), [doi:10.1016/j.physd.2026.135415](https://doi.org/10.1016/j.physd.2026.135415).
* P. H. B. Rossetto et al., *Robust Wada boundaries and entropy scaling in pp-wave
  spacetimes*, Chaos 36, 063136 (2026),
  [doi:10.1063/5.0336333](https://doi.org/10.1063/5.0336333).
* D. S. Fernandez, J. M. Seoane, M. A. F. Sanjuan,
  [arXiv:2202.03069](https://arxiv.org/abs/2202.03069); J. Aguirre, J. C. Vallejo,
  M. A. F. Sanjuan, Phys. Rev. E 64, 066208 (2001).
* Ph. de Forcrand, O. Jahn, *The baryon static potential from lattice QCD*,
  [arXiv:hep-ph/0502039](https://arxiv.org/abs/hep-ph/0502039).
* H. Suganuma, T. T. Takahashi, F. Okiharu, H. Ichie, *Lattice QCD study for the
  interquark force in three-quark and multi-quark systems*,
  [arXiv:hep-lat/0412026](https://arxiv.org/abs/hep-lat/0412026).
* *Baryonic flux in quenched and two-flavor dynamical QCD after Abelian projection*,
  [arXiv:hep-lat/0401026](https://arxiv.org/abs/hep-lat/0401026).
* C. Bernard, T. DeGrand, R. L. Sugar, *Zero temperature string breaking in lattice
  quantum chromodynamics*, [arXiv:hep-lat/0103012](https://arxiv.org/abs/hep-lat/0103012).
* *String breaking* (QCD with two flavours),
  [arXiv:hep-lat/0512018](https://arxiv.org/abs/hep-lat/0512018); SESAM,
  *Observation of string breaking in QCD*,
  [arXiv:hep-lat/0505012](https://arxiv.org/abs/hep-lat/0505012).
* S. Braun et al., *Negative absolute temperature for motional degrees of freedom*,
  Science 339, 52 (2013), [arXiv:1211.0545](https://arxiv.org/abs/1211.0545).
* G. Aarts, *Testing a realistic quark-gluon plasma*, Physics 7, 86,
  <https://physics.aps.org/articles/v7/86>; HotQCD crossover temperature as quoted in
  [arXiv:2312.12857](https://arxiv.org/abs/2312.12857).
* A. Jaffe, E. Witten, *Quantum Yang-Mills theory*, in *The Millennium Prize
  Problems*, Clay Mathematics Institute,
  <https://www.claymath.org/library/monographs/MPPc.pdf>.
* S. G. Matinyan, G. K. Savvidy, N. G. Ter-Arutyunyan-Savvidy, Sov. Phys. JETP 53,
  421 (1981); further chaos-in-Yang-Mills context in
  [arXiv:0710.5691](https://arxiv.org/abs/0710.5691) and
  [arXiv:0903.3990](https://arxiv.org/abs/0903.3990).
